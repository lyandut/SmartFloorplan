//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <unordered_set>

#include "Instance.hpp"
#include "QAPCluster.hpp"
#include "FloorplanBinPack.hpp"

using namespace fbp;
using namespace qapc;

class AdaptiveSelection {

public:

	AdaptiveSelection(const Environment &env, const Config &cfg) :
		_env(env), _cfg(cfg), _ins(_env), _cluster(_ins, _cfg.dimension),
		_gen(_cfg.random_seed), _best_fill_ratio(_cfg.init_fill_ratio) {}

	/// 候选宽度定义
	struct CandidateWidthObj {
		int value;
		int iter;
		unique_ptr<FloorplanBinPack> fbp_solver;
	};

	void run() {
		vector<Rect> src = _ins.get_rects();

		// 计算分组信息
		_cluster.cal_qap_sol(_cluster.cal_flow_matrix(_cfg.level_qapc_flow), _cluster.cal_distance_matrix(_cfg.level_qapc_dis));
		for (auto &rect : src) { rect.gid = _cluster.qap_sol.at(_cluster.part.at(rect.id)) - 1; }
		vector<vector<bool>> group_neighbors = _cluster.cal_group_neighbors();

		// Calculate the set of candidate widths W
		vector<int> candidate_widths;
		switch (_cfg.level_asa_cw) {
		case Config::LevelCandidateWidth::CombRotate:
			candidate_widths = cal_candidate_widths_on_combrotate(src);
			break;
		case Config::LevelCandidateWidth::CombShort:
			candidate_widths = cal_candidate_widths_on_combshort(src);
			break;
		case Config::LevelCandidateWidth::Interval:
			candidate_widths = cal_candidate_widths_on_interval(src);
			break;
		default:
			assert(false);
			break;
		}

		// 分支初始化iter=1
		vector<CandidateWidthObj> cw_objs;
		for (int bin_width : candidate_widths) {
			int bin_height = ceil(_ins.get_total_area() / (bin_width * _best_fill_ratio)); // 向上取整
			vector<Boundary> group_boundaries = _cluster.cal_group_boundaries(bin_width, bin_height);
			cw_objs.push_back({ bin_width, 1, unique_ptr<FloorplanBinPack>(
				new FloorplanBinPack(src, group_neighbors, group_boundaries, bin_width, _gen)) });
			cw_objs.back().fbp_solver->random_local_search(1, _cfg.level_fbp_gs);
			if (cw_objs.back().fbp_solver->get_fill_ratio() > _best_fill_ratio) {
				_best_fill_ratio = cw_objs.back().fbp_solver->get_fill_ratio();
				_best_dst = cw_objs.back().fbp_solver->get_dst();
				cout << _best_fill_ratio << endl;
			}
		}
		sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
			return lhs.fbp_solver->get_fill_ratio() < rhs.fbp_solver->get_fill_ratio(); });

		// 离散概率分布
		vector<int> probs; probs.reserve(cw_objs.size());
		for (int i = 1; i <= cw_objs.size(); ++i) { probs.push_back(2 * i); }
		discrete_distribution<> discrete_dist(probs.begin(), probs.end());

		// 迭代优化
		while (true) { // [todo] time limit is not exceeded
			CandidateWidthObj &picked_width = cw_objs[discrete_dist(_gen)];
			int new_bin_height = ceil(_ins.get_total_area() / (picked_width.value * _best_fill_ratio)); // 向上取整
			vector<Boundary> new_group_boundaries = _cluster.cal_group_boundaries(picked_width.value, new_bin_height);
			picked_width.iter = min(2 * picked_width.iter, _cfg.ub_iter);
			picked_width.fbp_solver->update_group_boundaries(new_group_boundaries);
			picked_width.fbp_solver->random_local_search(picked_width.iter, _cfg.level_fbp_gs);
			if (picked_width.fbp_solver->get_fill_ratio() > _best_fill_ratio) {
				_best_fill_ratio = picked_width.fbp_solver->get_fill_ratio();
				_best_dst = picked_width.fbp_solver->get_dst();
				cout << _best_fill_ratio << endl;
			}
			sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
				return lhs.fbp_solver->get_fill_ratio() < rhs.fbp_solver->get_fill_ratio(); });
		}
	}

private:
	/// 基于排列组合生成候选宽度组合，考虑所有组合及旋转有`2c_n1 + 2²c_n2 + 2³c_n3 + ...`中情况
	/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
	vector<int> cal_candidate_widths_on_combrotate(const vector<Rect> &src, double alpha = 1.05, int miniterms = 3, int maxiterms = 6) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		unordered_set<int> candidate_widths;
		vector<int> rects(src.size());
		iota(rects.begin(), rects.end(), 0);
		for (int k = miniterms; k <= maxiterms; ++k) {
			utils::Combination kth_comb_rects(rects, k);
			vector<int> comb_rects, ncomb_rects;
			while (kth_comb_rects.next_combination(comb_rects, ncomb_rects)) {
				for (int kk = 1; kk <= comb_rects.size() / 2; ++kk) { // 组合数性质计算一半即可.
					utils::Combination kth_rotated_rects(comb_rects, kk);
					vector<int> rotated_rects, nrotated_rects;
					while (kth_rotated_rects.next_combination(rotated_rects, nrotated_rects)) {
						int rcw = 0, nrcw = 0;
						for (int r : rotated_rects) { rcw += src.at(r).height; }
						for (int nr : nrotated_rects) { rcw += src.at(nr).width; }
						if (rcw >= min_cw && rcw <= max_cw) { candidate_widths.insert(rcw); }
						for (int nr : rotated_rects) { nrcw += src.at(nr).width; }
						for (int r : nrotated_rects) { nrcw += src.at(r).height; }
						if (nrcw >= min_cw && nrcw <= max_cw) { candidate_widths.insert(nrcw); }
					}
				}
			}
		}
		return vector<int>(candidate_widths.begin(), candidate_widths.end());
	}

	/// 基于排列组合生成候选宽度组合，仅考虑短边的组合
	/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
	vector<int> cal_candidate_widths_on_combshort(const vector<Rect> &src, double alpha = 1.05, int miniterms = 3, int maxiterms = 6) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		unordered_set<int> candidate_widths;
		vector<int> rects(src.size());
		iota(rects.begin(), rects.end(), 0);
		for (int k = miniterms; k <= maxiterms; ++k) {
			utils::Combination kth_comb_rects(rects, k);
			vector<int> comb_rects, ncomb_rects;
			while (kth_comb_rects.next_combination(comb_rects, ncomb_rects)) {
				int cw = 0, ncw = 0;
				for (int r : comb_rects) { cw += src.at(r).width; }
				if (cw >= min_cw && cw <= max_cw) { candidate_widths.insert(cw); }
				for (int nr : ncomb_rects) { ncw += src.at(nr).width; }
				if (ncw >= min_cw && ncw <= max_cw) { candidate_widths.insert(ncw); }
			}
		}
		return vector<int>(candidate_widths.begin(), candidate_widths.end());
	}

	/// 在区间[W_min, W_max]内，等距地生成候选宽度
	/// 论文alpha取值范围：[1.0, 1.05, 1.1, 1.15, 1.2]
	vector<int> cal_candidate_widths_on_interval(const vector<Rect> &src, double alpha = 1.05, int interval = 1) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		vector<int> candidate_widths;
		candidate_widths.reserve(max_cw - min_cw + 1);
		for (int cw = min_cw; cw <= max_cw; cw += interval) { candidate_widths.push_back(cw); }
		return candidate_widths;
	}

private:
	const Environment &_env;
	const Config &_cfg;

	Instance _ins;
	QAPCluster _cluster;
	default_random_engine _gen;

	float _best_fill_ratio;
	vector<Rect> _best_dst;
};