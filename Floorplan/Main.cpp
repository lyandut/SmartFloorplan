#include "AdaptiveSelection.hpp"

void test_qap() { qap::test_qap(); }

void test_floorplan_bin_pack(const Instance &ins, const QAPCluster &cluster) {
	float dead_ratio = 0.5f;
	int bin_width = ceil(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	printf("Calculating group info...\n");
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	default_random_engine gen(random_device{}());
	fbp::FloorplanBinPack fbp_solver(src, group_neighbors, group_boundaries, bin_width, gen);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, Config::LevelHeuristicSearch::MinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, Config::LevelGroupSearch::Selfishly);

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
}

int main(int argc, char **argv) {

	Environment env("GSRC", "H", "n300");
	//Environment env("MCNC", "H", "ami49");

	Config cfg;
	cfg.dimension = 5;
	cfg.random_seed = random_device{}();
	cfg.init_fill_ratio = 0.5f;
	cfg.ub_iter = 9999;
	cfg.level_asa_cw = Config::LevelCandidateWidth::Interval;
	cfg.level_qapc_flow = Config::LevelFlow::Kway;
	cfg.level_qapc_dis = Config::LevelDistance::ManhattanDis;
	cfg.level_fbp_gs = Config::LevelGroupSearch::None;

	AdaptiveSelection asa(env, cfg);
	asa.run();

	system("pause");
	return 0;
}
