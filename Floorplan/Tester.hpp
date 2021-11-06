//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <string>
#include <sstream>
#include <metis.h>

#include "QAPSolver.hpp"
#include "AdaptiveSelecter.hpp" 


namespace test {

	// GSRC算例有初始排版，记录利用率和线长
	void record_gsrc_init_sol() {
		for (auto& gsrc : ins_list) {
			if (gsrc.first == "MCNC") { continue; }
			Environment env(gsrc.first, "H", gsrc.second);
			Instance ins(env);

			int bin_width = 0, bin_height = 0;
			DisjointRects disjoint_rects;
			for (auto& r : ins.get_rects(false)) {
				assert(disjoint_rects.add(r));
				bin_width = max(bin_width, r.x + r.width);
				bin_height = max(bin_height, r.y + r.height);
			}

			double hpwl_block = 0, hpwl_terminal = 0;
			for (auto& net : ins.get_netlist()) {
				double max_x = 0, min_x = numeric_limits<double>::max();
				double max_y = 0, min_y = numeric_limits<double>::max();
				for (int b : net.block_list) {
					double pin_x = ins.get_blocks().at(b).x_coordinate + ins.get_blocks().at(b).width * 0.5;
					double pin_y = ins.get_blocks().at(b).y_coordinate + ins.get_blocks().at(b).height * 0.5;
					max_x = max(max_x, pin_x);
					min_x = min(min_x, pin_x);
					max_y = max(max_y, pin_y);
					min_y = min(min_y, pin_y);
				}
				hpwl_block += max_x - min_x + max_y - min_y;
				for (int t : net.terminal_list) {
					double pad_x = ins.get_terminals().at(t).x_coordinate;
					double pad_y = ins.get_terminals().at(t).y_coordinate;
					max_x = max(max_x, pad_x);
					min_x = min(min_x, pad_x);
					max_y = max(max_y, pad_y);
					min_y = min(min_y, pad_y);
				}
				hpwl_terminal += max_x - min_x + max_y - min_y;
			}

			ofstream log_file("Instance/GSRC/HARD/GSRC_init.csv", ios::app);
			log_file.seekp(0, ios::end);
			if (log_file.tellp() <= 0) {
				log_file << "Instance,Area,FillRatio,WireLengthBlock,WireLengthTerminal" << endl;
			}
			log_file << env._ins_name << "," << bin_width * bin_height << ","
				<< 1.0 * ins.get_total_area() / (bin_width * bin_height) << ","
				<< hpwl_block << "," << hpwl_terminal << endl;
		}
	}

	void test_floorplan_packer() {
		Environment env("GSRC", "H", "n300");
		Instance ins(env);
		vector<Rect> src = ins.get_rects();
		default_random_engine gen(cfg.random_seed);
		double dead_ratio = 1.05;
		int bin_width = ceil(sqrt(dead_ratio * ins.get_total_area()));

		printf("Perform the packing...\n");

		vector<shared_ptr<FloorplanPacker>> fbp_solvers;
		QAPCluster qap_cluster(ins, min(cfg.dimension, static_cast<int>(round(sqrt(ins.get_block_num())))));
		fbp_solvers.emplace_back(make_shared<RandomLocalSearcher>(ins, src, bin_width, qap_cluster._graph, gen));
		fbp_solvers.emplace_back(make_shared<BeamSearcher>(ins, src, bin_width, qap_cluster._graph, gen));
		for_each(fbp_solvers.begin(), fbp_solvers.end(), [&](auto& fbp_solver) {
			fbp_solver->run(1, cfg.alpha, cfg.beta, cfg.level_fbp_wl, cfg.level_fbp_dist);
			if (fbp_solver->get_dst().size() != ins.get_block_num()) { printf("Failed!\n"); }
			else { printf("Successful! Fill Ratio: %.2f%%\n", fbp_solver->get_fill_ratio()); }
		});
	}

}

namespace metis {

	using namespace std;

	vector<idx_t> func(vector<idx_t>& xadj, vector<idx_t>& adjncy, vector<idx_t>& adjwgt, decltype(METIS_PartGraphKway)* METIS_PartGraphFunc) {
		idx_t nVertices = xadj.size() - 1; // 节点数
		idx_t nEdges = adjncy.size() / 2;  // 边数
		idx_t nWeights = 1;                // 节点权重维数
		idx_t nParts = 2;                  // 子图个数≥2
		idx_t objval;                      // 目标函数值
		vector<idx_t> part(nVertices, 0);  // 划分结果

		int ret = METIS_PartGraphFunc(&nVertices, &nWeights, xadj.data(), adjncy.data(),
			NULL, NULL, adjwgt.data(), &nParts, NULL,
			NULL, NULL, &objval, part.data());

		if (ret != rstatus_et::METIS_OK) { cout << "METIS_ERROR" << endl; }
		cout << "METIS_OK" << endl;
		cout << "objval: " << objval << endl;
		for (unsigned part_i = 0; part_i < part.size(); part_i++) {
			cout << part_i + 1 << " " << part[part_i] << endl;
		}

		return part;
	}

	void test() {
		ifstream ingraph("Instance/METIS/graph_b.txt");

		int vexnum, edgenum;
		string line;
		getline(ingraph, line);
		istringstream tmp(line);
		tmp >> vexnum >> edgenum;
		vector<idx_t> xadj(0);
		vector<idx_t> adjncy(0); // 压缩图表示
		vector<idx_t> adjwgt(0); // 节点权重

		idx_t a, w;
		for (int i = 0; i < vexnum; i++) {
			xadj.push_back(adjncy.size());
			getline(ingraph, line);
			istringstream tmp(line);
			while (tmp >> a >> w) {
				adjncy.push_back(a - 1); // 节点id从0开始
				adjwgt.push_back(w);
			}
		}
		xadj.push_back(adjncy.size());
		ingraph.close();

		vector<idx_t> part = func(xadj, adjncy, adjwgt, METIS_PartGraphRecursive);
		//vector<idx_t> part = func(xadj, adjncy, adjwgt, METIS_PartGraphKway);

		ofstream outpartition("Instance/METIS/partition_b.txt");
		for (int i = 0; i < part.size(); i++) { outpartition << i + 1 << " " << part[i] << endl; }
		outpartition.close();
	}

}


namespace qap {

	void load_problem(int& n, type_matrix& a, type_matrix& b, long& best_objective) {
		cin >> best_objective >> n;
		a = new long* [n + 1];
		b = new long* [n + 1];
		for (int i = 1; i <= n; i = i + 1) {
			a[i] = new long[n + 1];
			b[i] = new long[n + 1];
		}

		for (int i = 1; i <= n; i = i + 1)
			for (int j = 1; j <= n; j = j + 1)
				cin >> a[i][j];
		for (int i = 1; i <= n; i = i + 1)
			for (int j = 1; j <= n; j = j + 1)
				cin >> b[i][j];
	}

	void load_problem_from_datfile(int& n, type_matrix& a, type_matrix& b, long& best_objective) {
		ifstream ifs("Instance/QAP/tai20a.dat");
		if (!ifs.is_open()) { return; }

		best_objective = 703482;
		ifs >> n;
		a = new long* [n + 1];
		b = new long* [n + 1];
		for (int i = 1; i <= n; i = i + 1) {
			a[i] = new long[n + 1];
			b[i] = new long[n + 1];
		}

		for (int i = 1; i <= n; i = i + 1)
			for (int j = 1; j <= n; j = j + 1)
				ifs >> a[i][j];
		for (int i = 1; i <= n; i = i + 1)
			for (int j = 1; j <= n; j = j + 1)
				ifs >> b[i][j];
	}

	void test() {
		//load_problem(n, a, b, best_objective);
		load_problem_from_datfile(n, a, b, best_objective);

		run_bls(n, a, b, best_objective);

		// check cost
		int check_cost = 0;
		for (int i = 1; i <= n; ++i)
			for (int j = 1; j <= n; ++j)
				check_cost += a[i][j] * b[solution[i]][solution[j]];
		if (check_cost != cost) { cout << "QAP_OK" << endl; }
		cout << "QAP_OK" << endl;
		cout << "objval: " << best_objective << endl;

		delete[] solution;
		for (int i = 1; i <= n; i = i + 1) {
			delete[] a[i];
			delete[] b[i];
		}
		delete[] a; delete[] b;
	}

}
