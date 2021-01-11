//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <memory>
#include <unordered_set>

#include "QAPCluster.hpp"
#include "RandomLocalSearcher.hpp"
#include "BeamSearcher.hpp"

using namespace qapc;
using namespace fbp;

class AdaptiveSelecter {

	/// 候选宽度定义
	struct CandidateWidth {
		int value;
		int iter;
		shared_ptr<FloorplanPacker> fbp_solver;
	};

public:

	AdaptiveSelecter() = delete;

	AdaptiveSelecter(const Environment &env, const Config &cfg) :
		_env(env), _cfg(cfg), _ins(_env), _gen(_cfg.random_seed), _start(clock()),
		_cluster(_ins, min(_cfg.dimension, static_cast<int>(round(sqrt(_ins.get_block_num()))))),
		_objective(numeric_limits<double>::max()) {}

	void run() {
		vector<Rect> src = _ins.get_rects();

		// Calculate the set of candidate widths W
		vector<int> candidate_widths;
		switch (_cfg.level_asa_cw) {
		case Config::LevelCandidateWidth::Interval:
			candidate_widths = cal_candidate_widths_on_interval(src);
			break;
		case Config::LevelCandidateWidth::Sqrt:
			candidate_widths = cal_candidate_widths_on_sqrt(src);
			break;
		case Config::LevelCandidateWidth::CombRotate:
			candidate_widths = cal_candidate_widths_on_combrotate(src);
			break;
		case Config::LevelCandidateWidth::CombShort:
			candidate_widths = cal_candidate_widths_on_combshort(src);
			break;
		default: assert(false); break;
		}

		// 初始化离散概率分布
		vector<int> probs; probs.reserve(candidate_widths.size());
		for (int i = 1; i <= candidate_widths.size(); ++i) { probs.push_back(2 * i); }
		discrete_distribution<> discrete_dist(probs.begin(), probs.end());
		// 初始化均匀分布
		uniform_int_distribution<> uniform_dist(0, candidate_widths.size() - 1);

		if (_cfg.level_asa_fbp == Config::LevelFloorplanPacker::RandomLocalSearch)
			random_local_search(src, candidate_widths, discrete_dist, uniform_dist);
		else
			beam_search(src, candidate_widths, discrete_dist, uniform_dist);
	}

	void random_local_search(vector<Rect> &src, vector<int> &candidate_widths,
		discrete_distribution<> &discrete_dist, uniform_int_distribution<> &uniform_dist) {
		// 分组信息①：分配`gid`
		if (_cfg.level_rls_qapc == Config::LevelQAPCluster::On) {
			_cluster.cal_qap_sol(
				_cluster.cal_flow_matrix(_cfg.level_qapc_flow),
				_cluster.cal_distance_matrix(_cfg.level_qapc_dist)
			);
			for (auto &rect : src) { rect.gid = _cluster.qap_sol.at(_cluster.part.at(rect.id)) - 1; }
		}
		else { assert(_cfg.level_rls_gs == Config::LevelGroupSearch::NoGroup); }

		_start = clock(); // 不计算qap调用时间

		// 分支初始化iter=1
		vector<CandidateWidth> cw_objs; cw_objs.reserve(candidate_widths.size());
		for (int bin_width : candidate_widths) {
			cw_objs.push_back({ bin_width, 1, make_shared<RandomLocalSearcher>(_ins, src, bin_width, _cluster._graph, _gen) });
			auto rls_solver = dynamic_pointer_cast<RandomLocalSearcher>(cw_objs.back().fbp_solver);
			rls_solver->run(1, _cfg.alpha, _cfg.beta, _cfg.level_fbp_wl, _cfg.level_fbp_dist); // 首次迭代默认LevelGroupSearch::NoGroup
			check_cwobj(cw_objs.back());

			// 分组信息②：计算`boundaries`和`neighbors`
			if (_cfg.level_rls_qapc == Config::LevelQAPCluster::On) {
				vector<Boundary> boundaries = _cluster.cal_group_boundaries(bin_width, rls_solver->get_skyline_height());
				vector<vector<bool>> neighbors = _cluster.cal_group_neighbors();
				rls_solver->set_group_boundaries(boundaries);
				rls_solver->set_group_neighbors(neighbors);
			}
		}
		// 降序排列，越后面的选中概率越大
		sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
			return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });

		// 迭代优化 
		while (static_cast<double>(clock() - _start) / CLOCKS_PER_SEC < _cfg.ub_time) {
			// 疏散性：90%概率选择，10%随机选择
			CandidateWidth &picked_width = _gen() % 10 ? cw_objs[discrete_dist(_gen)] : cw_objs[uniform_dist(_gen)];
			picked_width.iter = min(2 * picked_width.iter, _cfg.ub_iter);
			auto picked_solver = dynamic_pointer_cast<RandomLocalSearcher>(picked_width.fbp_solver);
			picked_solver->run(picked_width.iter, _cfg.alpha, _cfg.beta, _cfg.level_fbp_wl, _cfg.level_fbp_dist, _cfg.level_rls_gs);
			check_cwobj(picked_width);
			sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
				return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });
		}
	}

	void beam_search(vector<Rect> &src, vector<int> &candidate_widths,
		discrete_distribution<> &discrete_dist, uniform_int_distribution<> &uniform_dist) {
		vector<CandidateWidth> cw_objs; cw_objs.reserve(candidate_widths.size());
		for (int bin_width : candidate_widths) {
			cw_objs.push_back({ bin_width, 1, make_shared<BeamSearcher>(_ins, src, bin_width, _cluster._graph, _gen) });
			auto bs_solver = dynamic_pointer_cast<BeamSearcher>(cw_objs.back().fbp_solver);
			bs_solver->run(1, _cfg.alpha, _cfg.beta, _cfg.level_fbp_wl, _cfg.level_fbp_dist);
			check_cwobj(cw_objs.back());
		}
		sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
			return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });

		while (static_cast<double>(clock() - _start) / CLOCKS_PER_SEC < _cfg.ub_time) {
			CandidateWidth &picked_width = _gen() % 10 ? cw_objs[discrete_dist(_gen)] : cw_objs[uniform_dist(_gen)];
			picked_width.iter = min(2 * picked_width.iter, _cfg.ub_iter);
			auto picked_solver = dynamic_pointer_cast<BeamSearcher>(picked_width.fbp_solver);
			picked_solver->run(picked_width.iter, _cfg.alpha, _cfg.beta, _cfg.level_fbp_wl, _cfg.level_fbp_dist);
			check_cwobj(picked_width);
			sort(cw_objs.begin(), cw_objs.end(), [](auto &lhs, auto &rhs) {
				return lhs.fbp_solver->get_objective() > rhs.fbp_solver->get_objective(); });
		}
	}

	void record_fp(const string &fp_path) const {
		ofstream fp_file(fp_path);
		for (auto &r : _dst) {
			fp_file << _ins.get_blocks().at(r.id).name << " " << r.x << " " << r.y << endl;
		}
		fp_file << endl;
		for (auto &t : _ins.get_terminals()) {
			fp_file << t.name << " " << t.x_coordinate << " " << t.y_coordinate << endl;
		}
	}

	void draw_fp(string html_path, bool draw_wire = false) const {
		utils_visualize_drawer::Drawer html_drawer(html_path, _ins.get_fixed_width() * 2, _ins.get_fixed_height() * 2);
		for (auto &r : _dst) { html_drawer.rect(r.x, r.y, r.width, r.height); }
		for (auto &t : _ins.get_terminals()) { html_drawer.circle(t.x_coordinate, t.y_coordinate); }
		if (draw_wire) {
			for (auto &net : _ins.get_netlist()) {
				html_drawer.rc.next();
				for (int i = 0; i < net.block_list.size(); ++i) {
					for (int j = i + 1; j < net.block_list.size(); ++j) {
						html_drawer.wire(
							_dst[net.block_list[i]].x + _dst[net.block_list[i]].width * 0.5,
							_dst[net.block_list[i]].y + _dst[net.block_list[i]].height * 0.5,
							_dst[net.block_list[j]].x + _dst[net.block_list[j]].width * 0.5,
							_dst[net.block_list[j]].y + _dst[net.block_list[j]].height * 0.5
						);
					}
				}
			}
		}
	}

	void draw_ins() const {
		ifstream ifs(_env.ins_html_path());
		if (ifs.good()) { return; }
		utils_visualize_drawer::Drawer html_drawer(_env.ins_html_path(), _ins.get_fixed_width(), _ins.get_fixed_height());
		for (auto &r : _ins.get_blocks()) { html_drawer.rect(r.x_coordinate, r.y_coordinate, r.width, r.height); }
		for (auto &t : _ins.get_terminals()) { html_drawer.circle(t.x_coordinate, t.y_coordinate); }
	}

	void record_log() const {
		ofstream log_file(_env.log_path(), ios::app);
		log_file.seekp(0, ios::end);
		if (log_file.tellp() <= 0) {
			log_file << "Instance,"
				"Alpha,Area,FillRatio,WHRatio,"
				"Beta,WireLength,Objective,CheckObj,"
				"Duration,Iteration,RandomSeed,"
				"LevelFloorplanPacker,LevelWireLength,LevelObjDist,LevelQAPCluster,"
				"LevelGroupSearch,LevelGraphConnect,LevelFlow,LevelDist" << endl;
		}
		log_file << _env._ins_name << ","
			<< _cfg.alpha << "," << _best_area << "," << _best_fillratio << "," << _best_whratio << ","
			<< _cfg.beta << "," << _best_wirelength << "," << _objective << "," << check_dst() << ","
			<< _duration << "," << _iteration << "," << _cfg.random_seed << "," << _cfg << endl;
	}

private:
	/// 基于排列组合生成候选宽度组合，考虑所有组合及旋转有`2c_n1 + 2²c_n2 + 2³c_n3 + ...`中情况
	/// 不需要从k=1开始计算组合数，通过[miniterms, maxiterms]参数控制；论文设置maxiterms=3,4,6
	vector<int> cal_candidate_widths_on_combrotate(const vector<Rect> &src, int miniterms = 3, int maxiterms = 6, double alpha = 1.05) {
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
				if (cw >= min_cw && cw <= max_cw && cw * _ins.get_fixed_height() > _ins.get_total_area()) { candidate_widths.insert(cw); }
				for (int nr : ncomb_rects) { ncw += src.at(nr).width; }
				if (ncw >= min_cw && ncw <= max_cw && _ins.get_fixed_height() > _ins.get_total_area()) { candidate_widths.insert(ncw); }
			}
		}
		return vector<int>(candidate_widths.begin(), candidate_widths.end());
	}

	/// 在区间[W_min, W_max]内，等距地生成候选宽度
	vector<int> cal_candidate_widths_on_interval(const vector<Rect> &src, int interval = 1) {
		int min_cw = 0, max_cw = 0;
		for_each(src.begin(), src.end(), [&](const auto &rect) {
			min_cw = max(min_cw, rect.height);
			max_cw += rect.height;
		});

		vector<int> candidate_widths; candidate_widths.reserve(max_cw - min_cw + 1);
		for (int cw = min_cw; cw <= max_cw; cw += interval) { candidate_widths.push_back(cw); }
		return candidate_widths;
	}

	/// 开平方限制长宽比，削减分支数目
	vector<int> cal_candidate_widths_on_sqrt(const vector<Rect> &src, int interval = 1) {
		int min_cw = floor(_cfg.lb_scale * sqrt(_ins.get_total_area()));
		int max_cw = ceil(_cfg.ub_scale * sqrt(_ins.get_total_area()));
		for_each(src.begin(), src.end(), [&](const auto &rect) { min_cw = max(min_cw, rect.height); });
		max_cw = max(max_cw, min_cw);

		vector<int> candidate_widths; candidate_widths.reserve(max_cw - min_cw + 1);
		for (int cw = min_cw; cw <= max_cw; cw += interval) { candidate_widths.push_back(cw); }
		return candidate_widths;
	}

	/// 更新历史最优解
	void check_cwobj(const CandidateWidth &cw_obj) {
		if (cw_obj.fbp_solver->get_objective() < _objective) {
			_duration = static_cast<double>(clock() - _start) / CLOCKS_PER_SEC;
			_iteration = cw_obj.iter;
			_objective = cw_obj.fbp_solver->get_objective();
			_best_area = cw_obj.fbp_solver->get_area();
			_best_fillratio = cw_obj.fbp_solver->get_fill_ratio();
			int cw_height = _best_area / cw_obj.value;
			_best_whratio = 1.0 * max(cw_obj.value, cw_height) / min(cw_obj.value, cw_height);
			_best_wirelength = cw_obj.fbp_solver->get_wirelength();
			_dst = cw_obj.fbp_solver->get_dst();
		}
	}

	/// 解的合法性检查
	bool check_dst() const {
		DisjointRects disjoint_rects;
		disjoint_rects.rects.reserve(_dst.size());
		for (int i = 0; i < _dst.size(); ++i) {
			if (min(_dst.at(i).width, _dst.at(i).height) != min(_ins.get_blocks().at(i).width, _ins.get_blocks().at(i).height) ||
				max(_dst.at(i).width, _dst.at(i).height) != max(_ins.get_blocks().at(i).width, _ins.get_blocks().at(i).height)) {
				fprintf(stderr, "id:%d, name:%s, has wrong width/height.\n", i, _ins.get_blocks().at(i).name.c_str());
				return false;
			}
			if (!disjoint_rects.add(_dst.at(i))) { // 重叠检测
				fprintf(stderr, "id:%d, name:%s, has overlap error.\n", i, _ins.get_blocks().at(i).name.c_str());
				return false;
			}
		}
		fprintf(stdout, "success!\n");
		return true;
	}

private:
	const Environment &_env;
	const Config &_cfg;

	Instance _ins;
	QAPCluster _cluster;
	default_random_engine _gen;
	clock_t _start;
	double _duration;
	int _iteration;

	double _objective;
	int _best_area;
	double _best_fillratio;
	double _best_whratio;
	double _best_wirelength;
	vector<Rect> _dst;
};
