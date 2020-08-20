#include <cmath>

#include "Instance.hpp"
#include "FloorplanBinPack.hpp"
#include "QAPCluster.hpp"
#include "../QAPSolver/QAPSolver.hpp"

/*
void test_floorplan_bin_pack() {

	//Environment env("GSRC", "H", "n10");
	Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	float dead_ratio = 0.1452f;
	int bin_width = sqrt(ins.get_total_area() * (1 + dead_ratio));
	int bin_height = bin_width;

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(ins, bin_width, bin_height);

	printf("Perform the packing...\n");
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	//fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);
	fbp_solver.random_local_search(1);
	fbp_solver.random_local_search(2);

	// Test success or failure.
	//if (src.empty()) { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
	//else { printf("Failed!\n"); }
	//printf("Done. All rectangles packed.\n");
} */

// void test_qap() { qap::test_bls(); }

void test_cluster() {
	Environment env("GSRC", "H", "n300");
	//Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	float dead_ratio = 0.5;
	int bin_width = sqrt(ins.get_total_area() * (1 + dead_ratio));
	int bin_height = bin_width;

	QAPCluster cluster(ins);
	cluster.cal_flow_matrix(5);
}

int main(int argc, char **argv) {

	//test_floorplan_bin_pack();

	//test_qap();

	test_cluster();


	system("pause");
	return 0;
}
