#include <unordered_set>

#include "Instance.hpp"
#include "FloorplanBinPack.hpp"
#include "QAPCluster.hpp"

void test_floorplan_bin_pack(const Instance &ins, const QAPCluster &cluster) {
	float dead_ratio = 0.5f;
	int bin_width = ceil(sqrt(ins.get_total_area() * (1 + dead_ratio)));
	int bin_height = bin_width;

	// 计算分组信息
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)); }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(src, group_neighbors, group_boundaries, bin_width, bin_height);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
}

/// 候选宽度定义
struct CandidateWidth {
	int value;
	double fill_ratio;
	unique_ptr<fbp::FloorplanBinPack> fbp_solver;
};

void adaptive_selection(const Instance &ins, const QAPCluster &cluster) {
	// 计算分组信息
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)); }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();

	// 贪心算法或者人为设置一个初始填充率
	float best_fill_ratio = 0.9f;
	vector<rbp::Rect> best_dst;

	// [todo] 穷举所有宽度组合: 2c_n1 + 2²c_n2 + 2³c_n3 + ...
	unordered_set<int> visited_widths;
	vector<CandidateWidth> candidate_widths;
	for (auto &rect : src) {
		int tmp_width;
		if (visited_widths.find(tmp_width) != visited_widths.end()) { continue; }
		visited_widths.insert(tmp_width);
		candidate_widths.push_back({}); // [todo] 是否有必要使用智能指针？
	}

	for (auto &cw : candidate_widths) {
		int bin_width = cw.value;
		int bin_height = ceil(ins.get_total_area() / (bin_width * best_fill_ratio)); // 向上取整
		vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

		printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
		cw.fbp_solver->update_height(bin_height, group_boundaries);

		printf("Perform the packing...\n");
		cw.fbp_solver->random_local_search(1);
	}

	// Sort each width in W by the increasing filling ratio
	sort();

	while (true) {

	}
}

int main(int argc, char **argv) {
	//qap::test_qap();

	//Environment env("GSRC", "H", "n300");
	Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	int dimension = 5;
	QAPCluster cluster(ins, dimension);
	vector<vector<int>> flow_matrix = cluster.cal_flow_matrix(QAPCluster::LevelMetis::Kway);
	vector<vector<int>> distance_matrix = cluster.cal_distance_matrix(QAPCluster::LevelDistance::EuclideanDis);
	cluster.cal_qap_sol(flow_matrix, distance_matrix);

	//test_floorplan_bin_pack(ins, cluster);

	adaptive_selection(ins, cluster);

	system("pause");
	return 0;
}
