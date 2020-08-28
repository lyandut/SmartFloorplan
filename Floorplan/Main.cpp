#include "AdaptiveSelection.hpp"

using namespace qapc;

void test_qap() { qap::test_qap(); }

void test_floorplan_bin_pack(const Instance &ins, const QAPCluster &cluster) {
	default_random_engine gen(random_device{}());
	float dead_ratio = 0.5f;
	int bin_width = ceil(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	printf("Calculating group info...\n");
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(src, group_neighbors, group_boundaries, bin_width, gen);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
}

int main(int argc, char **argv) {

	Environment env("GSRC", "H", "n300");
	//Environment env("MCNC", "H", "ami49");

	Configuration cfg;
	cfg.dimension = 5;
	cfg.random_seed = random_device{}();
	cfg.init_fill_ratio = 0.5f;
	cfg.ub_iter = 9999;
	cfg.level_cw = Configuration::LevelCandidateWidth::Interval;
	cfg.level_flow = QAPCluster::LevelMetis::Kway;
	cfg.level_dis = QAPCluster::LevelDistance::ManhattanDis;
	cfg.level_gs = FloorplanBinPack::LevelGroupSearch::LevelNone;

	AdaptiveSelection asa(env, cfg);
	asa.run();

	system("pause");
	return 0;
}
