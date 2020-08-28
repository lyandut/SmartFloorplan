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
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();
	vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	fbp::FloorplanBinPack fbp_solver(src, group_neighbors, group_boundaries, bin_width);

	printf("Perform the packing...\n");
	vector<fbp::Rect> dst;
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);

	if (dst.size() != ins.get_block_num()) { printf("Failed!\n"); }
	else { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
}

/// 基于排列组合生成候选宽度组合
/// Calculate the set of candidate widths W: 2c_n1 + 2²c_n2 + 2³c_n3 + ...
/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
vector<int> calculate_candidate_widths_on_combination(const vector<fbp::Rect> &src, int area, double alpha = 1.05, int miniterms = 2, int maxiterms = 6) {
	int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
	int max_cw = floor(sqrt(area) * alpha);
	unordered_set<int> candidate_widths;
	vector<int> rects(src.size());
	iota(rects.begin(), rects.end(), 0);
	for (int k = miniterms; k <= maxiterms; ++k) {
		utils::Combination kth_comb_rects(rects, k);
		vector<int> comb_rects, ncomb_rects;
		while (kth_comb_rects.next_combination(comb_rects, ncomb_rects)) {
			/// Ⅰ. 不考虑旋转，短边之和为一个候选宽度组合
			int cw = 0, ncw = 0;
			for (int r : comb_rects) { cw += src.at(r).width; }
			if (cw >= min_cw && cw <= max_cw) { candidate_widths.insert(cw); }
			for (int nr : ncomb_rects) { ncw += src.at(nr).width; }
			if (ncw >= min_cw && ncw <= max_cw) { candidate_widths.insert(ncw); }

			/// Ⅱ. 考虑旋转，复杂度太高，不可避免重复计算
			//for (int kk = 1; kk <= comb_rects.size() / 2; ++kk) { // 利用组合数性质计算一半即可.
			//	utils::Combination kth_rotated_rects(comb_rects, kk);
			//	vector<int> rotated_rects, nrotated_rects;
			//	while (kth_rotated_rects.next_combination(rotated_rects, nrotated_rects)) {
			//		int rcw = 0, nrcw = 0;
			//		for (int r : rotated_rects) { rcw += src.at(r).height; }
			//		for (int nr : nrotated_rects) { rcw += src.at(nr).width; }
			//		if (rcw >= min_cw && rcw <= max_cw) { candidate_widths.insert(rcw); }
			//		for (int nr : rotated_rects) { nrcw += src.at(nr).width; }
			//		for (int r : nrotated_rects) { nrcw += src.at(r).height; }
			//		if (nrcw >= min_cw && nrcw <= max_cw) { candidate_widths.insert(nrcw); }
			//	}
			//}
			//for (int kk = 1; kk <= ncomb_rects.size() / 2; ++kk) {
			//	utils::Combination kth_rotated_rects(ncomb_rects, kk);
			//	vector<int> rotated_rects, nrotated_rects;
			//	while (kth_rotated_rects.next_combination(rotated_rects, nrotated_rects)) {
			//		int rcw = 0, nrcw = 0;
			//		for (int r : rotated_rects) { rcw += src.at(r).height; }
			//		for (int nr : nrotated_rects) { rcw += src.at(nr).width; }
			//		if (rcw >= min_cw && rcw <= max_cw) { candidate_widths.insert(rcw); }
			//		for (int nr : rotated_rects) { nrcw += src.at(nr).width; }
			//		for (int r : nrotated_rects) { nrcw += src.at(r).height; }
			//		if (nrcw >= min_cw && nrcw <= max_cw) { candidate_widths.insert(nrcw); }
			//	}
			//}

		}
	}
	return vector<int>(candidate_widths.begin(), candidate_widths.end());
}

/// 在区间[W_min, W_max]内，等距地生成候选宽度
/// alpha取值范围：[1.0, 1.05, 1.1, 1.15, 1.2]
vector<int> calculate_candidate_widths_on_interval(const vector<rbp::Rect> &src, int area, double alpha = 1.05, int interval = 1) {
	int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
	int max_cw = floor(sqrt(area) * alpha);
	vector<int> candidate_widths;
	candidate_widths.reserve(max_cw - min_cw + 1);
	for (int cw = min_cw; cw <= max_cw; cw += interval) { candidate_widths.push_back(cw); }
	return candidate_widths;
}

/// 候选宽度定义
struct CandidateWidthObj {
	int value;
	int iter;
	unique_ptr<fbp::FloorplanBinPack> fbp_solver;
};

void adaptive_selection(const Instance &ins, const QAPCluster &cluster, fbp::FloorplanBinPack::LevelGroupSearch method) {
	// 计算分组信息
	vector<rbp::Rect> src = ins.get_rects();
	for (auto &rect : src) { rect.gid = cluster.qap_sol.at(cluster.part.at(rect.id)) - 1; }
	vector<vector<bool>> group_neighbors = cluster.cal_group_neighbors();

	// 参数设置
	int ub_iter = 9999;           // 最大迭代次数
	float best_fill_ratio = 0.5f; // 贪心算法或者人为设置一个初始填充率
	vector<fbp::Rect> best_dst;

	cout << "==== Calculate the set of candidate widths W ====" << endl;
	//vector<int> candidate_widths = calculate_candidate_widths_on_combination(src, ins.get_total_area());
	vector<int> candidate_widths = calculate_candidate_widths_on_interval(src, ins.get_total_area());
	cout << candidate_widths.size() << endl;

	vector<CandidateWidthObj> cw_objs;
	for (int bin_width : candidate_widths) {
		int bin_height = ceil(ins.get_total_area() / (bin_width * best_fill_ratio)); // 向上取整
		vector<rbp::Boundary> group_boundaries = cluster.cal_group_boundaries(bin_width, bin_height);
		cw_objs.push_back({ bin_width, 1, unique_ptr<fbp::FloorplanBinPack>(
			new fbp::FloorplanBinPack(src, group_neighbors, group_boundaries, bin_width)) });
		cw_objs.back().fbp_solver->random_local_search(1, method);
		if (cw_objs.back().fbp_solver->get_fill_ratio() > best_fill_ratio) {
			best_fill_ratio = cw_objs.back().fbp_solver->get_fill_ratio();
			best_dst = cw_objs.back().fbp_solver->get_dst();
			cout << best_fill_ratio << endl;
		}
	}
	// Sort each width in W by the increasing filling ratio
	sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
		return lhs.fbp_solver->get_fill_ratio() < rhs.fbp_solver->get_fill_ratio(); });

	while (true) { // [todo] time limit is not exceeded
		static default_random_engine gen(random_device{}());
		vector<int> probs;
		probs.reserve(cw_objs.size());
		for (int i = 1; i <= cw_objs.size(); ++i) { probs.push_back(2 * i); }
		discrete_distribution<> discrete_dist(probs.begin(), probs.end());

		CandidateWidthObj &picked_width = cw_objs[discrete_dist(gen)];
		int new_bin_height = ceil(ins.get_total_area() / (picked_width.value * best_fill_ratio)); // 向上取整
		vector<rbp::Boundary> new_group_boundaries = cluster.cal_group_boundaries(picked_width.value, new_bin_height);
		picked_width.iter = min(2 * picked_width.iter, ub_iter);
		picked_width.fbp_solver->update_group_boundaries(new_group_boundaries);
		picked_width.fbp_solver->random_local_search(picked_width.iter, method);
		if (picked_width.fbp_solver->get_fill_ratio() > best_fill_ratio) {
			best_fill_ratio = picked_width.fbp_solver->get_fill_ratio();
			best_dst = picked_width.fbp_solver->get_dst();
			cout << best_fill_ratio << endl;
		}
		sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
			return lhs.fbp_solver->get_fill_ratio() < rhs.fbp_solver->get_fill_ratio(); });
	}

	//return best found solution
	//best_dst;
	cout << best_fill_ratio << endl;
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

	adaptive_selection(ins, cluster, fbp::FloorplanBinPack::LevelGroupSearch::LevelNone);

	system("pause");
	return 0;
}
