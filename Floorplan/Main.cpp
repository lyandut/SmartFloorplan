//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#include "Tester.hpp"
#include "AdaptiveSelecter.hpp" 

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

		// record init gsrc
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

void run_single_ins(const string& ins_bench, const string& ins_name) {
	Environment env(ins_bench, "H", ins_name);
	AdaptiveSelecter asa(env, cfg);
	asa.run();
	asa.draw_ins();
	asa.record_fp(env.fp_path());
	asa.record_fp(env.fp_path_with_time());
	asa.draw_fp(env.fp_html_path(), cfg.beta);
	asa.draw_fp(env.fp_html_path_with_time(), cfg.beta);
	asa.record_log();
}

void run_all_ins() { for_each(ins_list.begin(), ins_list.end(), [](auto& ins) { run_single_ins(ins.first, ins.second); }); }

int main(int argc, char** argv) {

	//qap::test_qap();

	//metis::test_metis();

	//record_gsrc_init_sol();

	//test_floorplan_packer();

	//run_single_ins("GSRC", "n10");

	run_all_ins();

	return 0;
}
