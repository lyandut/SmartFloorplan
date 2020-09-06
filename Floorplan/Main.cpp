#include "AdaptiveSelection.hpp"

void test_qap() { qap::test_qap(); }

// [todo] 修改为无参函数，用于测试贪心算法
void test_floorplan_bin_pack(const Instance &ins, const QAPCluster &cluster) {
	double dead_ratio = 0.5;
	int bin_width = ceil(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	printf("Calculating group info...\n");
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	default_random_engine gen(random_device{}());
	fbp::FloorplanBinPack fbp_solver(ins, src, group_neighbors, group_boundaries, bin_width, gen);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, Config::LevelGroupSearch::Selfishly);

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
}

static Config cfg{
	random_device{}(),
	0.5,
	1,
	0,
	5,
	60 * 10,
	9999,
	Config::LevelCandidateWidth::Interval,
	Config::LevelGraphConnection::Direct,
	Config::LevelFlow::Kway,
	Config::LevelDistance::ManhattanDis,
	Config::LevelGroupSearch::NoGroup,
	Config::LevelWireLength::BlockOnly,
	Config::LevelObjNorm::NoNorm,
	Config::LevelHeuristicSearch::BottomLeftScore
};

static vector<pair<string, string>> ins_list{
	{"MCNC", "ami33"}, {"MCNC", "ami49"},{"MCNC", "apte"},{"MCNC", "hp"},
	//{"MCNC", "xerox"}, // xerox算例不适用fixed-outline
	{"GSRC", "n10"},{"GSRC", "n30"},{"GSRC", "n50"},{"GSRC", "n100"},{"GSRC", "n200"},{"GSRC", "n300"}
};

void run_all_ins() {
	for (auto &ins : ins_list) {
		Environment env(ins.first, "H", ins.second);
		AdaptiveSelection asa(env, cfg);
		asa.run();
		asa.record_fp(env.fp_path());
		asa.record_fp(env.fp_path_with_time());
		asa.record_sol(env.solution_path());
		//asa.record_sol(env.solution_path_with_time());
	}
}

void run_single_ins() {
	Config cfg;
	cfg.random_seed = random_device{}();
	cfg.init_fill_ratio = 0.5;
	cfg.alpha = 1;
	cfg.beta = 1 - cfg.alpha;
	cfg.dimension = 5;
	cfg.ub_time = 60 * 10;
	cfg.ub_iter = 9999;
	cfg.level_asa_cw = Config::LevelCandidateWidth::Interval;
	cfg.level_qapc_gc = Config::LevelGraphConnection::Direct;
	cfg.level_qapc_flow = Config::LevelFlow::Kway;
	cfg.level_qapc_dis = Config::LevelDistance::ManhattanDis;
	cfg.level_fbp_gs = Config::LevelGroupSearch::NeighborAll;
	cfg.level_fbp_wl = Config::LevelWireLength::BlockAndTerminal;
	cfg.level_fbp_norm = Config::LevelObjNorm::Average;

	Environment env("GSRC", "H", "n300");
	AdaptiveSelection asa(env, cfg);
	asa.run();
	asa.record_fp(env.fp_path());
	asa.record_fp(env.fp_path_with_time());
	asa.record_sol(env.solution_path());
	//asa.record_sol(env.solution_path_with_time());
}

int main(int argc, char **argv) {

	run_single_ins();

	//run_all_ins();

	system("pause");
	return 0;
}
