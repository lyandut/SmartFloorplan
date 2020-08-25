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
	int iter;
	float fill_ratio;
	unique_ptr<fbp::FloorplanBinPack> fbp_solver;
};

void adaptive_selection(const Instance &ins, const QAPCluster &cluster, fbp::FloorplanBinPack::LevelGroupSearch method) {
	// 计算分组信息
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)); }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();

	// 参数设置
	int iter_ub = 9999;           // 最大迭代次数
	float best_fill_ratio = 0.5f; // 贪心算法或者人为设置一个初始填充率
	vector<fbp::Rect> best_dst;

	// Calculate the set of candidate widths W: 2c_n1 + 2²c_n2 + 2³c_n3 + ...
	cout << "==== Calculate the set of candidate widths W ====" << endl;
	vector<int> comb_widths;
	vector<int> rects(src.size());
	iota(rects.begin(), rects.end(), 0);
	for (int k = 1; k <= rects.size(); ++k) {
		vector<int> comb_rects, ncomb_rects;
		while (utils::next_combination(rects, rects.size(), k, comb_rects, ncomb_rects)) {
			for (int kk = 1; kk <= comb_rects.size(); ++kk) {
				vector<int> rotated_rects, nrotated_rects;
				while (utils::next_combination(comb_rects, comb_rects.size(), kk, rotated_rects, nrotated_rects)) {
					comb_widths.resize(comb_widths.size() + 1, 0);
					for (int r : rotated_rects) { comb_widths.back() += src.at(r).height; }
					for (int nr : nrotated_rects) { comb_widths.back() += src.at(nr).width; }
				}
			}
		}
	}
	cout << comb_widths.size() << endl;
	cout << unordered_set<int>(comb_widths.begin(), comb_widths.end()).size() << endl;

	//unordered_set<int> visited_widths;
	//vector<CandidateWidth> candidate_widths;
	//for (int cw : comb_widths) {

	//	int bin_width = cw; // [todo]
	//	if (visited_widths.find(bin_width) != visited_widths.end()) { continue; }
	//	visited_widths.insert(bin_width);

	//	int bin_height = ceil(ins.get_total_area() / (bin_width * best_fill_ratio)); // 向上取整
	//	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);
	//	candidate_widths.push_back({ bin_width, 1, best_fill_ratio, unique_ptr<fbp::FloorplanBinPack>(
	//		new fbp::FloorplanBinPack(src, group_neighbors, group_boundaries, bin_width, bin_height)) });
	//	candidate_widths.back().fbp_solver->random_local_search(1, method);
	//	candidate_widths.back().fill_ratio = candidate_widths.back().fbp_solver->get_fill_ratio();
	//	if (candidate_widths.back().fill_ratio > best_fill_ratio) {
	//		best_fill_ratio = candidate_widths.back().fill_ratio;
	//		best_dst = candidate_widths.back().fbp_solver->get_dst();
	//	}
	//}
	//// Sort each width in W by the increasing filling ratio
	//sort(candidate_widths.begin(), candidate_widths.end(), [](auto &lhs, auto &rhs) { return lhs.fill_ratio < rhs.fill_ratio; });

	//while (true) { // [todo] time limit is not exceeded
	//	static default_random_engine gen(random_device{}());
	//	vector<int> probs;
	//	probs.reserve(candidate_widths.size());
	//	for (int i = 1; i <= candidate_widths.size(); ++i) { probs.push_back(2 * i); }
	//	discrete_distribution<> discrete_dist(probs.begin(), probs.end());

	//	CandidateWidth &picked_width = candidate_widths[discrete_dist(gen)];
	//	int new_bin_height = ceil(ins.get_total_area() / (picked_width.value * best_fill_ratio)); // 向上取整
	//	vector<rbp::Boundary> new_group_boundaries = cluster.cal_group_boundaries(picked_width.value, new_bin_height);
	//	picked_width.iter = min(2 * picked_width.iter, iter_ub);
	//	picked_width.fbp_solver->update_bin_height(new_bin_height, new_group_boundaries);
	//	picked_width.fbp_solver->random_local_search(picked_width.iter, method);
	//	picked_width.fill_ratio = picked_width.fbp_solver->get_fill_ratio();
	//	if (picked_width.fill_ratio > best_fill_ratio) {
	//		best_fill_ratio = picked_width.fill_ratio;
	//		best_dst = picked_width.fbp_solver->get_dst();
	//	}
	//	sort(candidate_widths.begin(), candidate_widths.end(), [](auto &lhs, auto &rhs) { return lhs.fill_ratio < rhs.fill_ratio; });
	//}

	// return best found solution
	// best_dst;
	// best_fill_ratio;
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

	adaptive_selection(ins, cluster, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);


	system("pause");
	return 0;
}
