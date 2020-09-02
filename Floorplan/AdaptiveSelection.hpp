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
		_env(env), _cfg(cfg), _ins(_env), _gen(_cfg.random_seed), _start(clock()),
		_cluster(_ins, min(_cfg.dimension, static_cast<int>(floor(sqrt(_ins.get_block_num()))))), // block数目不足以分组则求上界
		_objective(numeric_limits<double>::max()), _best_fill_ratio(_cfg.init_fill_ratio) {}

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

		//_start = clock(); // 不计算qap调用时间

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
			int lb_height = ceil(_ins.get_total_area() / (bin_width * _best_fill_ratio)); // bin下界高度，仅用于确定分组
			vector<Boundary> group_boundaries = _cluster.cal_group_boundaries(bin_width, lb_height);
			cw_objs.push_back({ bin_width, 1, unique_ptr<FloorplanBinPack>(
				new FloorplanBinPack(_ins, src, group_neighbors, group_boundaries, bin_width, _gen)) });
			cw_objs.back().fbp_solver->random_local_search(1, _cfg.alpha, _cfg.beta, _cfg.level_fbp_gs, _cfg.level_fbp_wl, _cfg.level_fbp_norm);
			check_cwobj(cw_objs.back());
		}
		// 降序排列，越后面的选中概率越大
		sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
			return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });

		// 初始化离散概率分布
		vector<int> probs; probs.reserve(cw_objs.size());
		for (int i = 1; i <= cw_objs.size(); ++i) { probs.push_back(2 * i); }
		discrete_distribution<> discrete_dist(probs.begin(), probs.end());

		// 迭代优化 
		while ((clock() - _start) / static_cast<double>(CLOCKS_PER_SEC) < _cfg.ub_time) {
			CandidateWidthObj &picked_width = cw_objs[discrete_dist(_gen)];
			int new_lb_height = ceil(_ins.get_total_area() / (picked_width.value * _best_fill_ratio));
			vector<Boundary> new_group_boundaries = _cluster.cal_group_boundaries(picked_width.value, new_lb_height);
			picked_width.iter = min(2 * picked_width.iter, _cfg.ub_iter);
			picked_width.fbp_solver->update_group_boundaries(new_group_boundaries);
			picked_width.fbp_solver->random_local_search(picked_width.iter, _cfg.alpha, _cfg.beta, _cfg.level_fbp_gs, _cfg.level_fbp_wl, _cfg.level_fbp_norm);
			check_cwobj(picked_width);
			sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
				return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });
		}
	}

	void record_fp(string fp_path) const {
		ofstream fp_file(fp_path);
		for (auto &rect : _best_dst) {
			fp_file << _ins.get_blocks().at(rect.id).name << " "
				<< rect.x << " "
				<< rect.y << endl;
		}
		fp_file << endl;
		for (auto &terminal : _ins.get_terminals()) {
			fp_file << terminal.name << " "
				<< terminal.x_coordinate << " "
				<< terminal.y_coordinate << endl;
		}
		fp_file.close();
	}

	void record_sol(string sol_path) const {
		ofstream log_file(sol_path, ios::app);
		log_file.seekp(0, ios::end);
		if (log_file.tellp() <= 0) {
			log_file << "Instance,"
				"Alpha(AreaWeight),Area,FillRatio,(Beta)WireWeight,WireLength,"
				"Objective,Duration,RandomSeed,Dimension,"
				"LevelCandidateWidth,LevelGraphConnection,LevelFlow,LevelDistance,LevelGroupSearch,LevelWireLength,LevelObjNorm"
				<< endl;
		}
		log_file << _env._ins_name << ","
			<< _cfg.alpha << "," << _best_area << "," << _best_fill_ratio << "," << _cfg.beta << "," << _best_wirelength << ","
			<< _objective << "," << _best_duration << "," << _cfg.random_seed << "," << _cluster._dimension << ","
			<< _cfg << endl;
	}

private:
	/// 基于排列组合生成候选宽度组合，考虑所有组合及旋转有`2c_n1 + 2²c_n2 + 2³c_n3 + ...`中情况
	/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
	vector<int> cal_candidate_widths_on_combrotate(const vector<Rect> &src, int miniterms = 3, int maxiterms = 6, double alpha = 1.05) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		//int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		int max_cw = _ins.get_fixed_width();
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
						if (rcw >= min_cw && rcw <= max_cw && rcw * _ins.get_fixed_height() > _ins.get_total_area()) {
							candidate_widths.insert(rcw);
						}
						for (int nr : rotated_rects) { nrcw += src.at(nr).width; }
						for (int r : nrotated_rects) { nrcw += src.at(r).height; }
						if (nrcw >= min_cw && nrcw <= max_cw && nrcw * _ins.get_fixed_height() > _ins.get_total_area()) {
							candidate_widths.insert(nrcw);
						}
					}
				}
			}
		}
		return vector<int>(candidate_widths.begin(), candidate_widths.end());
	}

	/// 基于排列组合生成候选宽度组合，仅考虑短边的组合
	/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
	vector<int> cal_candidate_widths_on_combshort(const vector<Rect> &src, int miniterms = 3, int maxiterms = 6, double alpha = 1.05) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		//int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		int max_cw = _ins.get_fixed_width();
		unordered_set<int> candidate_widths;
		vector<int> rects(src.size());
		iota(rects.begin(), rects.end(), 0);
		for (int k = miniterms; k <= maxiterms; ++k) {
			utils::Combination kth_comb_rects(rects, k);
			vector<int> comb_rects, ncomb_rects;
			while (kth_comb_rects.next_combination(comb_rects, ncomb_rects)) {
				int cw = 0, ncw = 0;
				for (int r : comb_rects) { cw += src.at(r).width; }
				if (cw >= min_cw && cw <= max_cw && cw * _ins.get_fixed_height() > _ins.get_total_area()) { candidate_widths.insert(cw); }
				for (int nr : ncomb_rects) { ncw += src.at(nr).width; }
				if (ncw >= min_cw && ncw <= max_cw && _ins.get_fixed_height() > _ins.get_total_area()) { candidate_widths.insert(ncw); }
			}
		}
		return vector<int>(candidate_widths.begin(), candidate_widths.end());
	}

	/// 在区间[W_min, W_max]内，等距地生成候选宽度
	/// packing论文不考虑fixed-outline限制，alpha取值范围：[1.0, 1.05, 1.1, 1.15, 1.2]
	vector<int> cal_candidate_widths_on_interval(const vector<Rect> &src, int interval = 1, double alpha = 1.05) {
		int min_cw = max_element(src.begin(), src.end(), [](auto &lhs, auto &rhs) { return lhs.height < rhs.height; })->height;
		//int max_cw = floor(sqrt(_ins.get_total_area()) * alpha);
		int max_cw = _ins.get_fixed_width();
		vector<int> candidate_widths;
		candidate_widths.reserve(max_cw - min_cw + 1);
		for (int cw = min_cw; cw <= max_cw; cw += interval) {
			if (cw * _ins.get_fixed_height() > _ins.get_total_area()) { candidate_widths.push_back(cw); }
		}
		return candidate_widths;
	}

	/// 检查cw_obj的RLS结果
	void check_cwobj(const CandidateWidthObj &cw_obj) {
		if (cw_obj.fbp_solver->get_area() / cw_obj.value > _ins.get_fixed_height()) { // 当前解无效
			cw_obj.fbp_solver->reset_objective();
		}
		if (cw_obj.fbp_solver->get_objective() < _objective) {
			_objective = cw_obj.fbp_solver->get_objective();
			_best_area = cw_obj.fbp_solver->get_area();
			_best_fill_ratio = cw_obj.fbp_solver->get_fill_ratio();
			_best_wirelength = cw_obj.fbp_solver->get_wirelength();
			_best_dst = cw_obj.fbp_solver->get_dst();
			_best_duration = (clock() - _start) / static_cast<double>(CLOCKS_PER_SEC);
		}
	}

private:
	const Environment &_env;
	const Config &_cfg;

	Instance _ins;
	QAPCluster _cluster;
	default_random_engine _gen;
	clock_t _start;

	int _best_area;
	double _best_fill_ratio; // _objective对应的填充率，不一定是最优填充率
	double _best_wirelength;
	double _objective;
	vector<Rect> _best_dst;
	double _best_duration;
};
