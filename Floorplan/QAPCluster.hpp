//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <cmath>
#include <metis.h>

#include "QAPSolver.hpp"
#include "Config.hpp"
#include "Instance.hpp"

namespace qapc {

	static double cal_distance(Config::LevelDist method, double x1, double y1, double x2, double y2) {
		double distance = 0;
		switch (method) {
		case Config::LevelDist::EuclideanDist:
			distance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)); break;
		case Config::LevelDist::ManhattanDist:
			distance = abs(x1 - x2) + abs(y1 - y2); break;
		case Config::LevelDist::ChebyshevDist:
			distance = max(abs(x1 - x2), abs(y1 - y2)); break;
		default:
			assert(false); break;
		}
		return distance;
	}

	/// 将算例预处理成QAP算例
	class QAPCluster {

	public:

		QAPCluster(const Instance &ins, int dimension) : _ins(ins), _dimension(dimension) {
			build_graph(Config::LevelGraphConnect::Direct);
		}

		/// 计算QAP的flow matrix
		vector<vector<int>> cal_flow_matrix(Config::LevelFlow method) {
			vector<vector<int>> flow_matrix;
			switch (method) {
			case Config::LevelFlow::Kway:
				metis_cluster(flow_matrix, _dimension * _dimension, METIS_PartGraphKway);
				break;
			case Config::LevelFlow::Recursive:
				metis_cluster(flow_matrix, _dimension * _dimension, METIS_PartGraphRecursive);
				break;
			default:
				assert(false);
				break;
			}
			return flow_matrix;
		}

		/// 计算QAP的distance matrix
		vector<vector<int>> cal_distance_matrix(Config::LevelDist method) {
			int node_num = _dimension * _dimension;
			_distance_nodes.resize(node_num);
			for (int x = 0; x < _dimension; ++x) {
				for (int y = 0; y < _dimension; ++y) {
					_distance_nodes[x + y * _dimension].first = x;
					_distance_nodes[x + y * _dimension].second = y;
				}
			}
			vector<vector<int>> distance_matrix(node_num, vector<int>(node_num, 0));
			for (int i = 0; i < node_num; ++i) {
				for (int j = i + 1; j < node_num; ++j) {
					distance_matrix[i][j] = round(cal_distance(method,
						_distance_nodes[i].first, _distance_nodes[i].second,
						_distance_nodes[j].first, _distance_nodes[j].second));
					distance_matrix[j][i] = distance_matrix[i][j];
				}
			}
			return distance_matrix;
		}

		/// 计算qap分组
		void cal_qap_sol(const vector<vector<int>> &flow_matrix, const vector<vector<int>> &distance_matrix) {
			qap::run(flow_matrix, distance_matrix, qap_sol);
		}

		/// 计算分组邻居信息
		vector<vector<bool>> cal_group_neighbors() const {
			int group_num = _dimension * _dimension;
			vector<vector<bool>> group_neighbors(group_num, vector<bool>(group_num, false));
			for (int gi = 0; gi < group_num; ++gi) {
				for (int gj = gi + 1; gj < group_num; ++gj) {
					if (1 == cal_distance(Config::LevelDist::ChebyshevDist, // 切比雪夫距离为1定义为邻居
						_distance_nodes.at(gi).first, _distance_nodes.at(gi).second,
						_distance_nodes.at(gj).first, _distance_nodes.at(gj).second)) {
						group_neighbors[gi][gj] = true;
						group_neighbors[gj][gi] = true;
					}
				}
			}
			return group_neighbors;
		}

		/// 计算分组边界信息
		vector<Boundary> cal_group_boundaries(int bin_width, int bin_height) const {
			int group_num = _dimension * _dimension;
			vector<Boundary> group_boundaries(group_num);
			double unit_width = 1.0 * bin_width / _dimension;
			double unit_height = 1.0 * bin_height / _dimension;
			for (int gi = 0; gi < group_num; ++gi) {
				group_boundaries[gi].x = unit_width * _distance_nodes.at(gi).first;
				group_boundaries[gi].y = unit_height * _distance_nodes.at(gi).second;
				group_boundaries[gi].width = unit_width;
				group_boundaries[gi].height = unit_height;
			}
			return group_boundaries;
		}

	private:
		/// 结合net_list还原出图
		void build_graph(Config::LevelGraphConnect method) {
			_graph.resize(_ins.get_block_num(), vector<int>(_ins.get_block_num(), 0));
			for (auto &net : _ins.get_netlist()) {
				for (int i = 0; i < net.block_list.size(); ++i) {
					for (int j = i + 1; j < net.block_list.size(); ++j) {
						int a = net.block_list[i];
						int b = net.block_list[j];
						_graph[a][b] += 1;
						_graph[b][a] += 1;
					}
				}
			}

			/// [todo] 考虑非直连的边，需要结合最短路径算法求跳数，但metis不支持浮点型权重
			if (method == Config::LevelGraphConnect::Indirect) {}
		}

		/// 邻接矩阵转压缩图（CSR）
		void transform_graph_2_csr(vector<idx_t> &xadj, vector<idx_t> &adjncy, vector<idx_t> &vwgt, vector<idx_t> &adjwgt) const {
			xadj.reserve(_graph.size() + 1);
			vwgt.reserve(_graph.size());
			adjncy.reserve(_graph.size() * _graph.size());
			adjwgt.reserve(_graph.size() * _graph.size());
			for (int i = 0; i < _graph.size(); ++i) {
				xadj.push_back(adjncy.size());
				vwgt.push_back(_ins.get_blocks().at(i).area);
				for (int j = 0; j < _graph.size(); ++j) {
					if (i == j) { continue; }
					if (_graph[i][j] >= 1) {
						adjncy.push_back(j);
						adjwgt.push_back(_graph[i][j]);
					}
				}
			}
			xadj.push_back(adjncy.size());
		}

		/// metis求解分组问题
		void metis_cluster(vector<vector<int>> &flow_matrix, int group_num, decltype(METIS_PartGraphKway) *METIS_PartGraphFunc) {
			flow_matrix.clear();
			flow_matrix.resize(group_num, vector<int>(group_num, 0));

			vector<idx_t> xadj, adjncy, vwgt, adjwgt;
			transform_graph_2_csr(xadj, adjncy, vwgt, adjwgt);
			idx_t nvtxs = xadj.size() - 1;
			idx_t ncon = 1;
			idx_t nparts = group_num;
			idx_t objval;
			part.resize(nvtxs, 0);
			int ret = METIS_PartGraphFunc(&nvtxs, &ncon, xadj.data(), adjncy.data(),
				vwgt.data(), NULL, adjwgt.data(), &nparts, NULL,
				NULL, NULL, &objval, part.data());
			assert(ret == rstatus_et::METIS_OK);

			for (unsigned part_i = 0; part_i < part.size(); ++part_i) {
				for (unsigned part_j = part_i + 1; part_j < part.size(); ++part_j) {
					if (part[part_i] == part[part_j]) { continue; }
					flow_matrix[part[part_i]][part[part_j]] += _graph[part_i][part_j];
					flow_matrix[part[part_j]][part[part_i]] += _graph[part_j][part_i];
				}
			}
		}

	public:
		vector<vector<int>> _graph;
		vector<int> part;
		vector<int> qap_sol;
		// `qap_sol[part[rect.id]] = distance_node_id = group_id`

	private:
		const Instance &_ins;
		const int _dimension;
		vector<pair<int, int>> _distance_nodes;
	};

}