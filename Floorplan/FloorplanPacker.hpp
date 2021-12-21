//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <list>
#include <numeric>

#include "Config.hpp"
#include "Instance.hpp"

namespace fbp {

	using namespace std;

	class FloorplanPacker {

	public:
		FloorplanPacker() = delete;

		FloorplanPacker(const Instance& ins, const vector<Rect>& src, int bin_width, default_random_engine& gen) :
			_ins(ins), _src(src), _bin_width(bin_width), _bin_height(INF),
			_graph(ins.get_block_num(), vector<int>(ins.get_block_num(), 0)),
			_gen(gen), _dst(), _objective(numeric_limits<double>::max()),
			_obj_area(numeric_limits<int>::max()), _obj_wirelength(numeric_limits<double>::max()) {
			for (auto& net : _ins.get_netlist()) {
				for (int i = 0; i < net.block_list.size(); ++i) {
					for (int j = i + 1; j < net.block_list.size(); ++j) {
						int a = net.block_list[i];
						int b = net.block_list[j];
						_graph[a][b] += 1;
						_graph[b][a] += 1;
					}
				}
			}
		}

		const vector<Rect>& get_dst() const { return _dst; }

		double get_objective() const { return _objective; }

		int get_area() const { return _obj_area; }

		double get_wirelength() const { return _obj_wirelength; }

		int get_bin_height() const { return _bin_height; }

		void set_bin_height(int height) { _bin_height = height; }

		virtual void run(int, double, double, Config::LevelWireLength, Config::LevelObjDist) = 0;

	protected:
		/// 目标函数
		double cal_objective(int area, double dist, double alpha, double beta) {
			return alpha * area + beta * dist;
		}

		/// 计算线长，默认引脚在中心
		double cal_wirelength(const vector<Rect>& dst, const vector<bool>& is_packed, double& dist,
			Config::LevelWireLength level_wl, Config::LevelObjDist level_dist) {
			double total_wirelength = 0;
			dist = 0;

			vector<pair<double, double>> pins(dst.size());
			unordered_set<int> nets;
			for (int i = 0; i < dst.size(); ++i) {
				if (!is_packed[i]) { continue; } // 只计算当前已放置的块
				pins[i].first = dst[i].x + dst[i].width * 0.5;
				pins[i].second = dst[i].y + dst[i].height * 0.5;
				nets.insert(_ins.get_blocks().at(i).net_ids.begin(), _ins.get_blocks().at(i).net_ids.end());
			}
			assert(!nets.empty());

			for (int nid : nets) {
				double max_x = 0, min_x = numeric_limits<double>::max();
				double max_y = 0, min_y = numeric_limits<double>::max();
				for (int bid : _ins.get_netlist().at(nid).block_list) {
					if (!is_packed[bid]) { continue; }
					max_x = max(max_x, pins[bid].first);
					min_x = min(min_x, pins[bid].first);
					max_y = max(max_y, pins[bid].second);
					min_y = min(min_y, pins[bid].second);
				}
				if (level_wl == Config::LevelWireLength::BlockAndTerminal) {
					for (int tid : _ins.get_netlist().at(nid).terminal_list) {
						double pad_x = _ins.get_terminals().at(tid).x_coordinate;
						double pad_y = _ins.get_terminals().at(tid).y_coordinate;
						max_x = max(max_x, pad_x);
						min_x = min(min_x, pad_x);
						max_y = max(max_y, pad_y);
						min_y = min(min_y, pad_y);
					}
				}
				double hpwl = max_x - min_x + max_y - min_y;
				total_wirelength += hpwl;
			}

			switch (level_dist) {
			case Config::LevelObjDist::WireLengthDist:
				dist = total_wirelength; // ① 线长
				break;
			case Config::LevelObjDist::SqrEuclideanDist:
				for (int i = 0; i < _graph.size(); ++i) {
					for (int j = i + 1; j < _graph.size(); ++j) {
						if (is_packed[i] && is_packed[j] && _graph[i][j]) {
							double dx = pins[i].first - pins[j].first;
							double dy = pins[i].second - pins[j].second;
							dist += dx * dx + dy * dy; // ② 两两矩形之间的欧氏平方距离：dx^2+dy^2
						}
					}
				}
				break;
			case Config::LevelObjDist::SqrManhattanDist:
				for (int i = 0; i < _graph.size(); ++i) {
					for (int j = i + 1; j < _graph.size(); ++j) {
						if (is_packed[i] && is_packed[j] && _graph[i][j]) {
							double dx = abs(pins[i].first - pins[j].first);
							double dy = abs(pins[i].second - pins[j].second);
							dist += (dx + dy) * (dx + dy); // ③ 两两矩形之间的曼哈顿平方距离：(dx+dy)^2
						}
					}
				}
				break;
			default:
				break;
			}

			return total_wirelength;
		}

		/// 更新最优解
		void update_objective(double objective, int area, double wirelength, const vector<Rect>& dst) {
			if (_objective > objective) {
				_objective = objective;
				_obj_area = area;
				_obj_wirelength = wirelength;
				_dst = dst;
			}
		}

		/// 将`SkylineNode`转换成`SkylineSpace`
		static SkylineSpace skyline_nodo_to_space(const Skyline& skyline, int skyline_index) {
			int hl, hr;
			if (skyline.size() == 1) {
				hl = hr = INF - skyline[skyline_index].y;
			}
			else if (skyline_index == 0) {
				hl = INF - skyline[skyline_index].y;
				hr = skyline[skyline_index + 1].y - skyline[skyline_index].y;
			}
			else if (skyline_index == skyline.size() - 1) {
				hl = skyline[skyline_index - 1].y - skyline[skyline_index].y;
				hr = INF - skyline[skyline_index].y;
			}
			else {
				hl = skyline[skyline_index - 1].y - skyline[skyline_index].y;
				hr = skyline[skyline_index + 1].y - skyline[skyline_index].y;
			}
			return { skyline[skyline_index].x, skyline[skyline_index].y, skyline[skyline_index].width, hl, hr };
		}

		/// 合并同一level的skyline节点.
		static void merge_skylines(Skyline& skyline) {
			skyline.erase(
				remove_if(skyline.begin(), skyline.end(), [](auto& rhs) { return rhs.width <= 0; }),
				skyline.end()
			);
			for (int i = 0; i < skyline.size() - 1; ++i) {
				if (skyline[i].y == skyline[i + 1].y) {
					skyline[i].width += skyline[i + 1].width;
					skyline.erase(skyline.begin() + i + 1);
					--i;
				}
			}
		}

	protected:
		// 输入
		const Instance& _ins;
		const vector<Rect>& _src;
		const int _bin_width;
		int _bin_height; // 非const，超出_bin_height提前剪枝
		vector<vector<int>> _graph; // 结合net_list还原出图，评估块之间连接的紧密程度
		default_random_engine& _gen;

		// 优化目标
		vector<Rect> _dst;
		double _objective;
		int _obj_area;
		double _obj_wirelength;
	};

}
