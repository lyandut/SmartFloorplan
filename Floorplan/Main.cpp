#include "Instance.hpp"
#include "FloorplanBinPack.hpp"
#include "QAPCluster.hpp"

void test_cluster() {
	Environment env("GSRC", "H", "n300");
	//Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	QAPCluster cluster(ins);
	vector<vector<int>> flow_matrix = cluster.cal_flow_matrix(5);
	vector<vector<int>> distance_matrix = cluster.cal_distance_matrix(5, QAPCluster::EuclideanDis);
	cluster.cal_qap_sol(flow_matrix, distance_matrix);
}

void test_floorplan_bin_pack() {
	//Environment env("GSRC", "H", "n10");
	Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	QAPCluster cluster(ins);
	vector<vector<int>> flow_matrix = cluster.cal_flow_matrix(5);
	vector<vector<int>> distance_matrix = cluster.cal_distance_matrix(5, QAPCluster::EuclideanDis);
	cluster.cal_qap_sol(flow_matrix, distance_matrix);

	float dead_ratio = 0.5f;
	int bin_width = round(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(ins, cluster, bin_width, bin_height, 5 * 5);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);

	if (dst.size() == ins.get_block_num()) { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
	else { printf("Failed!\n"); }
}

void run() {
	//Environment env("GSRC", "H", "n10");
	Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	QAPCluster cluster(ins);
	vector<vector<int>> flow_matrix = cluster.cal_flow_matrix(5);
	vector<vector<int>> distance_matrix = cluster.cal_distance_matrix(5, QAPCluster::EuclideanDis);
	cluster.cal_qap_sol(flow_matrix, distance_matrix);

	float dead_ratio = 0.1452f;
	int bin_width = round(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(ins, cluster, bin_width, bin_height);

	printf("Perform the packing...\n");
	fbp_solver.random_local_search(1);

	// Test success or failure.
	//if (src.empty()) { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
	//else { printf("Failed!\n"); }
	//printf("Done. All rectangles packed.\n");
}

int main(int argc, char **argv) {
	//qap::test_qap();

	//test_cluster();

	test_floorplan_bin_pack();

	system("pause");
	return 0;
}
