#include "AdaptiveSelection.hpp"

void test_qap() { qap::test_qap(); }

void test_floorplan_bin_pack() {
	int dimension = 5;
	double dead_ratio = 0.5;

	Environment env("GSRC", "H", "n300");
	Instance ins(env);
	QAPCluster cluster(ins, min(dimension, static_cast<int>(floor(sqrt(ins.get_block_num())))));

	int bin_width = ceil(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int lb_height = bin_width;

	printf("Calculating group info...\n");
	vector<rbp::Rect> src = ins.get_rects();
	cluster.cal_qap_sol(cluster.cal_flow_matrix(Config::LevelFlow::Kway), cluster.cal_distance_matrix(Config::LevelDistance::ManhattanDis));
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, lb_height);
	default_random_engine gen(random_device{}());
	fbp::FloorplanBinPack fbp_solver(ins, src, group_neighbors, group_boundaries, bin_width, gen);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinWasteFit);
	int envelope_area = fbp_solver.insert_bottom_left_score(dst, Config::LevelGroupSearch::NoGroup) * bin_width;

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Fill Ratio: %.2f%%\n", (float)ins.get_total_area() / envelope_area * 100.f); }
}

void run_single_ins() {
	Config cfg;
	cfg.random_seed = random_device{}();
	cfg.init_fill_ratio = 0.5;
	cfg.alpha = 1;
	cfg.beta = 1 - cfg.alpha;
	cfg.dimension = 5;
	cfg.ub_time = 60 * 1;
	cfg.ub_iter = 9999;
	cfg.level_asa_cw = Config::LevelCandidateWidth::Interval;
	cfg.level_qapc_gc = Config::LevelGraphConnection::Direct;
	cfg.level_qapc_flow = Config::LevelFlow::Kway;
	cfg.level_qapc_dis = Config::LevelDistance::ManhattanDis;
	cfg.level_fbp_gs = Config::LevelGroupSearch::NeighborAll;
	cfg.level_fbp_wl = Config::LevelWireLength::BlockAndTerminal;
	cfg.level_fbp_norm = Config::LevelObjNorm::Average;
	cfg.level_fbp_hs = Config::LevelHeuristicSearch::BottomLeftScore;

	Environment env("GSRC", "H", "n30");
	AdaptiveSelection asa(env, cfg);
	asa.run();
	asa.record_fp(env.fp_path());
	asa.record_fp(env.fp_path_with_time());
	asa.draw_html(env.html_path());
	asa.draw_html(env.html_path_with_time());
	asa.record_sol(env.solution_path());
}

static vector<pair<string, string>> ins_list{
	{"MCNC", "ami33"}, {"MCNC", "ami49"},{"MCNC", "apte"},{"MCNC", "hp"},
	//{"MCNC", "xerox"}, // xerox算例不适用fixed-outline
	{"GSRC", "n10"},{"GSRC", "n30"},{"GSRC", "n50"},{"GSRC", "n100"},{"GSRC", "n200"},{"GSRC", "n300"}
};

void run_all_ins() {
	Config cfg;
	cfg.random_seed = random_device{}();
	cfg.init_fill_ratio = 0.5;
	cfg.alpha = 1;
	cfg.beta = 1 - cfg.alpha;
	cfg.dimension = 5;
	cfg.ub_time = 60 * 30;
	cfg.ub_iter = 9999;
	cfg.level_asa_cw = Config::LevelCandidateWidth::Interval;
	cfg.level_qapc_gc = Config::LevelGraphConnection::Direct;
	cfg.level_qapc_flow = Config::LevelFlow::Kway;
	cfg.level_qapc_dis = Config::LevelDistance::ManhattanDis;
	cfg.level_fbp_gs = Config::LevelGroupSearch::NoGroup;
	cfg.level_fbp_wl = Config::LevelWireLength::BlockOnly;
	cfg.level_fbp_norm = Config::LevelObjNorm::NoNorm;
	cfg.level_fbp_hs = Config::LevelHeuristicSearch::BottomLeftScore;

	for (auto &ins : ins_list) {
		Environment env(ins.first, "H", ins.second);
		AdaptiveSelection asa(env, cfg);
		asa.run();
		asa.record_fp(env.fp_path());
		asa.record_fp(env.fp_path_with_time());
		asa.draw_html(env.html_path());
		asa.draw_html(env.html_path_with_time());
		asa.record_sol(env.solution_path());
	}
}

void record_gsrc_init_sol() {
	for (auto &gsrc : ins_list) {
		if (gsrc.first == "MCNC") { continue; }
		Environment env(gsrc.first, "H", gsrc.second);
		Instance ins(env);

		int min_bin_width = 0, min_bin_height = 0;
		double hpwl_block = 0, hpwl_terminal = 0;

		// draw init gsrc
		utils_visualize_drawer::Drawer html_drawer("Instance/GSRC/HARD/" + env._ins_name + ".html",
			ins.get_fixed_width(), ins.get_fixed_height());

		auto dst = ins.get_rects(false);
		utils_visualize_transform::dst_to_boxes(dst); // 检查是否重叠
		for (auto &b : ins.get_blocks()) {
			html_drawer.rect(b.x_coordinate, b.y_coordinate, b.width, b.height);
			min_bin_width = max(min_bin_width, b.x_coordinate + b.width);
			min_bin_height = max(min_bin_height, b.y_coordinate + b.height);
		}
		for (auto &t : ins.get_terminals()) {
			html_drawer.circle(t.x_coordinate, t.y_coordinate);
		}
		auto wire_block_boxes = utils_visualize_transform::wire_to_boxes(ins, dst, Config::LevelWireLength::BlockOnly);
		for (auto &w : wire_block_boxes) {
			//html_drawer.wire(w.min_corner().x(), w.min_corner().y(),
			//	w.max_corner().x() - w.min_corner().x(),
			//	w.max_corner().y() - w.min_corner().y());
			hpwl_block +=
				w.max_corner().x() - w.min_corner().x() +
				w.max_corner().y() - w.min_corner().y();
		}
		auto wire_terminal_boxes = utils_visualize_transform::wire_to_boxes(ins, dst, Config::LevelWireLength::BlockAndTerminal);
		for (auto &w : wire_terminal_boxes) {
			hpwl_terminal +=
				w.max_corner().x() - w.min_corner().x() +
				w.max_corner().y() - w.min_corner().y();
		}

		// record init gsrc
		ofstream log_file("Instance/GSRC/HARD/GSRC_init.csv", ios::app);
		log_file.seekp(0, ios::end);
		if (log_file.tellp() <= 0) {
			log_file << "Instance,Area,FillRatio,WireLengthBlock,WireLengthTerminal" << endl;
		}
		log_file << env._ins_name << ","
			<< min_bin_width * min_bin_height << ","
			<< 1.0*ins.get_total_area() / (min_bin_width * min_bin_height) << ","
			<< hpwl_block << "," << hpwl_terminal << endl;
	}
}


int main(int argc, char **argv) {

	//record_gsrc_init_sol();

	//test_floorplan_bin_pack();

	//run_single_ins();

	run_all_ins();

	system("pause");
	return 0;
}
