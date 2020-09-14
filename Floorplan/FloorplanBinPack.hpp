//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <list>
#include <random>
#include <numeric>

#include "Instance.hpp"
#include "SkylineBinPack.hpp"
#include "Visualizer.hpp"

namespace fbp {

	using namespace rbp;

	class FloorplanBinPack final : public SkylineBinPack {

	public:

		FloorplanBinPack() = delete;

		FloorplanBinPack(const Instance &ins, const vector<Rect> &src, const vector<vector<bool>> &group_neighbors,
			vector<Boundary> &group_boundaries, int bin_width, default_random_engine &gen) :
			_ins(ins), _src(src), _objective(numeric_limits<double>::max()),
			_group_neighbors(group_neighbors), _group_boundaries(group_boundaries),
			SkylineBinPack(bin_width, INF, false), // binHeight = INF; 转为完成SP
			_gen(gen), _uniform_dist(0, _src.size() - 1) {
			init_sort_rules();
		}

		const vector<Rect>& get_dst() const { return _dst; }

		int get_area() const { return _best_area; }

		double get_fill_ratio() const { return _best_fill_ratio; }

		double get_wirelength() const { return _best_wirelength; }

		double get_objective() const { return _objective; }

		void reset_objective() { _objective = numeric_limits<double>::max(); }

		/// 高度上界变化导致分组边界变化
		void update_group_boundaries(const vector<Boundary> &new_boundaries) { _group_boundaries = new_boundaries; }

		/// 基于binWidth进行随机局部搜索
		void random_local_search(int iter, double alpha, double beta,
			Config::LevelGroupSearch level_gs, Config::LevelWireLength level_wl, Config::LevelObjNorm level_norm) {
			// the first time to call RLS on W_k
			if (iter == 1) {
				for (auto &rule : _sort_rules) {
					_rects.assign(rule.sequence.begin(), rule.sequence.end());
					vector<Rect> target_dst;
					int target_area = insert_bottom_left_score(target_dst, level_gs) * binWidth;
					double target_wirelength = cal_wirelength(target_dst, level_wl);
					add_his_sol(target_area, target_wirelength);
					rule.target_objective = cal_objective(target_area, target_wirelength, alpha, beta, level_norm);
					check_rule_sol(rule, target_area, target_wirelength, target_dst);
				}
				// 降序排列，越后面的目标函数值越小选中概率越大
				sort(_sort_rules.begin(), _sort_rules.end(), [](auto &lhs, auto &rhs) { return lhs.target_objective > rhs.target_objective; });
			}

			// 迭代优化
			SortRule &picked_rule = _sort_rules[_discrete_dist(_gen)];
			for (int i = 1; i <= iter; ++i) {
				SortRule new_rule = picked_rule;
				int a = _uniform_dist(_gen);
				int b = _uniform_dist(_gen);
				while (a == b) { b = _uniform_dist(_gen); }
				swap(new_rule.sequence[a], new_rule.sequence[b]);
				_rects.assign(new_rule.sequence.begin(), new_rule.sequence.end());
				vector<Rect> target_dst;
				int target_area = insert_bottom_left_score(target_dst, level_gs) * binWidth;
				double target_wirelength = cal_wirelength(target_dst, level_wl);
				add_his_sol(target_area, target_wirelength);
				new_rule.target_objective = min(new_rule.target_objective, cal_objective(target_area, target_wirelength, alpha, beta, level_norm));
				if (new_rule.target_objective < picked_rule.target_objective) {
					picked_rule = new_rule;
					check_rule_sol(picked_rule, target_area, target_wirelength, target_dst);
				}
			}

			// 更新排序规则列表
			sort(_sort_rules.begin(), _sort_rules.end(), [](auto &lhs, auto &rhs) { return lhs.target_objective > rhs.target_objective; });
		}

		int insert_bottom_left_score(vector<Rect> &dst, Config::LevelGroupSearch method) {
			reset();
			int min_bin_height = 0;
			dst.clear();
			dst.reserve(_rects.size());
			while (!_rects.empty()) {
				auto bottom_skyline_iter = min_element(skyLine.begin(), skyLine.end(),
					[](auto &lhs, auto &rhs) { return lhs.y < rhs.y; });
				int best_skyline_index = distance(skyLine.begin(), bottom_skyline_iter);
				list<int> candidate_rects = get_candidate_rects(*bottom_skyline_iter, method);
				auto min_width_iter = min_element(candidate_rects.begin(), candidate_rects.end(),
					[this](int lhs, int rhs) { return _src.at(lhs).width < _src.at(rhs).width; });
				int min_width = _src.at(*min_width_iter).width;

				if (bottom_skyline_iter->width < min_width) { // 最小宽度不能满足，需要填坑
					if (best_skyline_index == 0) { skyLine[best_skyline_index].y = skyLine[best_skyline_index + 1].y; }
					else if (best_skyline_index == skyLine.size() - 1) { skyLine[best_skyline_index].y = skyLine[best_skyline_index - 1].y; }
					else { skyLine[best_skyline_index].y = min(skyLine[best_skyline_index - 1].y, skyLine[best_skyline_index + 1].y); }
					MergeSkylines();
					continue;
				}

				int best_rect_index = -1;
				Rect best_node = find_rect_for_skyline_bottom_left(best_skyline_index, candidate_rects, best_rect_index);

				// 所有矩形必能放下
				assert(best_rect_index != -1);

				// 执行放置
				debug_assert(disjointRects.Disjoint(best_node));
				debug_run(disjointRects.Add(best_node));
				if (best_node.x == skyLine[best_skyline_index].x) { // 靠左
					AddSkylineLevel(best_skyline_index, best_node);
				}
				else { // 靠右
					SkylineNode new_skyline{ best_node.x, best_node.y + best_node.height, best_node.width };
					assert(new_skyline.x + new_skyline.width <= binWidth);
					assert(new_skyline.y <= binHeight);
					skyLine.insert(skyLine.begin() + best_skyline_index + 1, new_skyline);
					skyLine[best_skyline_index].width -= best_node.width;
					MergeSkylines();
				}
				usedSurfaceArea += _src.at(best_rect_index).width * _src.at(best_rect_index).height;
				_rects.remove(best_rect_index);
				dst.push_back(best_node);
				min_bin_height = max(min_bin_height, best_node.y + best_node.height);
			}

			assert(usedSurfaceArea == _ins.get_total_area());
			debug_run(utils_visualize_transform::dst_to_boxes(dst)); // 可视化packing过程
			return min_bin_height;
		}

		/// [deprecated]
		int insert_greedy_fit(vector<Rect> &dst, Config::LevelHeuristicSearch method) {
			reset();
			int min_bin_height = 0;
			dst.clear();
			dst.reserve(_rects.size());
			while (!_rects.empty()) {
				Rect best_node;
				int best_score_1 = numeric_limits<int>::max();
				int best_score_2 = numeric_limits<int>::max();
				int best_skyline_index = -1;
				int best_rect_index = -1;

				for (int i = 0; i < skyLine.size(); ++i) {
					Rect new_node;
					int score_1, score_2, rect_index;
					switch (method) {
					case Config::LevelHeuristicSearch::MinHeightFit:
						new_node = find_rect_for_skyline_min_height(i, _rects, score_1, score_2, rect_index);
						debug_assert(disjointRects.Disjoint(new_node));
						break;
					case Config::LevelHeuristicSearch::MinWasteFit:
						new_node = find_rect_for_skyline_min_waste(i, _rects, score_1, score_2, rect_index);
						debug_assert(disjointRects.Disjoint(new_node));
						break;
					default:
						assert(false);
						break;
					}
					if (new_node.height != 0) {
						if (score_1 < best_score_1 || (score_1 == best_score_1 && score_2 < best_score_2)) {
							best_node = new_node;
							best_score_1 = score_1;
							best_score_2 = score_2;
							best_skyline_index = i;
							best_rect_index = rect_index;
						}
					}
				}

				// 所有矩形必能放下
				assert(best_rect_index != -1);

				// 执行放置
				debug_assert(disjointRects.Disjoint(best_node));
				debug_run(disjointRects.Add(best_node));
				AddSkylineLevel(best_skyline_index, best_node);
				usedSurfaceArea += _src.at(best_rect_index).width * _src.at(best_rect_index).height;
				_rects.remove(best_rect_index);
				dst.push_back(best_node);
				min_bin_height = max(min_bin_height, best_node.y + best_node.height);
			}

			assert(usedSurfaceArea == _ins.get_total_area());
			debug_run(utils_visualize_transform::dst_to_boxes(dst)); // 可视化packing过程
			return min_bin_height;
		}

		/// 最优解可视化（debug_run）
		void visualize_best_sol(Config::LevelWireLength method) {
			debug_run(_terminal_points = utils_visualize_transform::terminal_to_points(_ins));
			debug_run(_dst_boxes = utils_visualize_transform::dst_to_boxes(_dst));
			debug_run(_wire_boxes = utils_visualize_transform::wire_to_boxes(_ins, _dst, method));
		}

	private:
		/// 每次迭代重置usedSurfaceArea和skyLine
		void reset() { Init(binWidth, binHeight, useWasteMap); }

		/// 排序规则定义
		struct SortRule {
			vector<int> sequence;
			double target_objective;
		};

		/// 初始化排序规则列表
		void init_sort_rules() {
			vector<int> seq(_src.size());
			// 0_输入顺序
			iota(seq.begin(), seq.end(), 0);
			for (int i = 0; i < 5; ++i) { _sort_rules.push_back({ seq, numeric_limits<double>::max() }); }
			// 1_面积递减
			sort(_sort_rules[1].sequence.begin(), _sort_rules[1].sequence.end(), [this](int lhs, int rhs) {
				return (_src.at(lhs).width * _src.at(lhs).height) > (_src.at(rhs).width * _src.at(rhs).height); });
			// 2_高度递减
			sort(_sort_rules[2].sequence.begin(), _sort_rules[2].sequence.end(), [this](int lhs, int rhs) {
				return _src.at(lhs).height > _src.at(rhs).height; });
			// 3_宽度递减
			sort(_sort_rules[3].sequence.begin(), _sort_rules[3].sequence.end(), [this](int lhs, int rhs) {
				return _src.at(lhs).width > _src.at(rhs).width; });
			// 4_随机排序
			shuffle(_sort_rules[4].sequence.begin(), _sort_rules[4].sequence.end(), _gen);

			// 默认赋给_rect输入顺序，可用于测试贪心算法
			_rects.assign(_sort_rules[0].sequence.begin(), _sort_rules[0].sequence.end());

			// 离散概率分布初始化
			vector<int> probs; probs.reserve(_sort_rules.size());
			for (int i = 1; i <= _sort_rules.size(); ++i) { probs.push_back(2 * i); }
			_discrete_dist = discrete_distribution<>(probs.begin(), probs.end());
		}

		/// 基于分组策略挑选候选矩形，减小搜索规模
		list<int> get_candidate_rects(const SkylineNode &skyline, Config::LevelGroupSearch method) {
			int gid = 0;
			for (; gid < _group_boundaries.size(); ++gid) {
				if (skyline.x >= _group_boundaries[gid].x && skyline.y >= _group_boundaries[gid].y
					&& skyline.x < _group_boundaries[gid].x + _group_boundaries[gid].width
					&& skyline.y < _group_boundaries[gid].y + _group_boundaries[gid].height) {
					break;
				}
			}
			list<int> candidate_rects;
			switch (method) {
			case Config::LevelGroupSearch::Selfishly:
				for (int r : _rects) { // 必须顺序遍历_rects才能使得sort_rule生效
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
				}
				break;
			case Config::LevelGroupSearch::NeighborAll:
				for (int r : _rects) {
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
					else if (_group_neighbors[_src.at(r).gid][gid]) {
						candidate_rects.push_back(r);
					}
				}
				break;
			case Config::LevelGroupSearch::NeighborPartial:
				for (int r : _rects) {
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
					else if (_group_neighbors[_src.at(r).gid][gid] && _src.at(r).gid < gid) {
						candidate_rects.push_back(r);
					}
				}
				break;
			case Config::LevelGroupSearch::NoGroup: // 仅在纯优化面积时使用
				candidate_rects = _rects;
				break;
			default:
				assert(false);
				break;
			}
			// 如果没有符合策略的矩形，则返回所有矩形
			if (candidate_rects.empty()) { candidate_rects = _rects; }
			return candidate_rects;
		}

		/// 论文idea，基于最下/最左和打分策略
		Rect find_rect_for_skyline_bottom_left(int skyline_index, const list<int> &rects, int &best_rect) {
			int best_score = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node)); // 确保无解时返回高度为0
			for (int r : rects) {
				int x, score;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					int width = _src.at(r).width, height = _src.at(r).height;
					if (rotate) { swap(width, height); }
					if (score_rect_for_skyline_bottom_left(skyline_index, width, height, x, score)) {
						if (best_score < score) {
							best_score = score;
							best_rect = r;
							new_node = { _src.at(r).id, _src.at(r).gid, x, skyLine[skyline_index].y, width, height };
							debug_assert(disjointRects.Disjoint(new_node));
						}
					}
				}
			}
			// [todo] 未实现(d)(f)->(h)的退化情况
			if (best_score == 4) {}
			if (best_score == 2) {}
			if (best_score == 0) {}

			return new_node;
		}

		/// 基于最小高度，为当前skyline选择最佳放置的矩形
		Rect find_rect_for_skyline_min_height(int skyline_index, const list<int> &rects,
			int &best_height, int &best_width, int &best_rect) {
			best_height = numeric_limits<int>::max();
			best_width = numeric_limits<int>::max(); // 高度相同选择宽度较小的，[todo]可能选择宽度较大的比较好？
			best_rect = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node));
			for (int r : rects) {
				int y;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					int width = _src.at(r).width, height = _src.at(r).height;
					if (rotate) { swap(width, height); }
					if (RectangleFits(skyline_index, width, height, y)) {
						if (y + height < best_height || (y + height == best_height && width < best_width)) {
							best_height = y + height;
							best_width = width;
							best_rect = r;
							new_node = { _src.at(r).id, _src.at(r).gid, skyLine[skyline_index].x, y, width, height };
							debug_assert(disjointRects.Disjoint(new_node));
						}
					}
				}
			}
			return new_node;
		}

		/// 基于最小浪费
		Rect find_rect_for_skyline_min_waste(int skyline_index, const list<int> &rects,
			int &best_wasted_area, int &best_height, int &best_rect) {
			best_wasted_area = numeric_limits<int>::max();
			best_height = numeric_limits<int>::max(); // 废料面积相同选择高度小的
			best_rect = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node));
			for (int r : rects) {
				int y, wasted_area;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					int width = _src.at(r).width, height = _src.at(r).height;
					if (rotate) { swap(width, height); }
					if (RectangleFits(skyline_index, width, height, y, wasted_area)) {
						if (wasted_area < best_wasted_area || (wasted_area == best_wasted_area && y + height < best_height)) {
							best_wasted_area = wasted_area;
							best_height = y + height;
							best_rect = r;
							new_node = { _src.at(r).id, _src.at(r).gid, skyLine[skyline_index].x, y, width, height };
							debug_assert(disjointRects.Disjoint(new_node));
						}
					}
				}
			}
			return new_node;
		}

		/// space定义
		struct SkylineSpace {
			int x;
			int y;
			int width;
			int hl;
			int hr;
		};

		/// 打分策略
		bool score_rect_for_skyline_bottom_left(int skyline_index, int width, int height, int &x, int &score) {
			if (width > skyLine[skyline_index].width) { return false; }
			if (skyLine[skyline_index].y + height > binHeight) { return false; }

			SkylineSpace space = skyline_nodo_to_space(skyline_index);
			if (space.hl >= space.hr) {
				if (width == space.width && height == space.hl) { score = 7; }
				else if (width == space.width && height == space.hr) { score = 6; }
				else if (width == space.width && height > space.hl) { score = 5; }
				else if (width < space.width && height == space.hl) { score = 4; }
				else if (width == space.width && height < space.hl && height > space.hr) { score = 3; }
				else if (width < space.width && height == space.hr) { score = 2; } // 靠右
				else if (width == space.width && height < space.hr) { score = 1; }
				else if (width < space.width && height != space.hl) { score = 0; }
				else { return false; }

				if (score == 2) { x = skyLine[skyline_index].x + skyLine[skyline_index].width - width; }
				else { x = skyLine[skyline_index].x; }
			}
			else { // hl < hr
				if (width == space.width && height == space.hr) { score = 7; }
				else if (width == space.width && height == space.hl) { score = 6; }
				else if (width == space.width && height > space.hr) { score = 5; }
				else if (width < space.width && height == space.hr) { score = 4; } // 靠右
				else if (width == space.width && height < space.hr && height > space.hl) { score = 3; }
				else if (width < space.width && height == space.hl) { score = 2; }
				else if (width == space.width && height < space.hl) { score = 1; }
				else if (width < space.width && height != space.hr) { score = 0; } // 靠右
				else { return false; }

				if (score == 4 || score == 0) { x = skyLine[skyline_index].x + skyLine[skyline_index].width - width; }
				else { x = skyLine[skyline_index].x; }
			}
			if (x + width > binWidth) { return false; }

			return true;
		}

		SkylineSpace skyline_nodo_to_space(int skyline_index) {
			int hl, hr;
			if (skyLine.size() == 1) {
				hl = hr = binHeight - skyLine[skyline_index].y;
			}
			else if (skyline_index == 0) {
				hl = binHeight - skyLine[skyline_index].y;
				hr = skyLine[skyline_index + 1].y - skyLine[skyline_index].y;
			}
			else if (skyline_index == skyLine.size() - 1) {
				hl = skyLine[skyline_index - 1].y - skyLine[skyline_index].y;
				hr = binHeight - skyLine[skyline_index].y;
			}
			else {
				hl = skyLine[skyline_index - 1].y - skyLine[skyline_index].y;
				hr = skyLine[skyline_index + 1].y - skyLine[skyline_index].y;
			}
			return { skyLine[skyline_index].x, skyLine[skyline_index].y, skyLine[skyline_index].width, hl, hr };
		}

		/// 计算线长并可视化
		/// [todo] 默认引脚在中心，MCNC算例引脚在block边上
		double cal_wirelength(vector<Rect> &dst, Config::LevelWireLength method) {
			int total_wirelength = 0;
			sort(dst.begin(), dst.end(), [](auto &lhs, auto &rhs) { return lhs.id < rhs.id; }); // 这里改变了插入顺序
			for (auto &net : _ins.get_net_list()) {
				double max_x = 0, min_x = numeric_limits<double>::max();
				double max_y = 0, min_y = numeric_limits<double>::max();
				for (int b : net.block_list) {
					double pin_x = dst.at(b).x + dst.at(b).width / 2.0;
					double pin_y = dst.at(b).y + dst.at(b).height / 2.0;
					max_x = max(max_x, pin_x);
					min_x = min(min_x, pin_x);
					max_y = max(max_y, pin_y);
					min_y = min(min_y, pin_y);
				}
				if (method == Config::LevelWireLength::BlockAndTerminal) {
					for (int t : net.terminal_list) {
						double pad_x = _ins.get_terminals().at(t).x_coordinate;
						double pad_y = _ins.get_terminals().at(t).y_coordinate;
						max_x = max(max_x, pad_x);
						min_x = min(min_x, pad_x);
						max_y = max(max_y, pad_y);
						min_y = min(min_y, pad_y);
					}
				}
				double hpwl = max_x - min_x + max_y - min_y;
				total_wirelength += hpwl;
			}
			return total_wirelength;
		}

		void add_his_sol(int area, double wirelength) {
			// 控制size，防止内存溢出
			if (_his_area.size() >= 9999) { _his_area.pop_front(); }
			if (_his_wirelength.size() >= 9999) { _his_wirelength.pop_front(); }
			_his_area.push_back(area);
			_his_wirelength.push_back(wirelength);
		}

		/// 尝试不同的归一化方法：1.最小面积/线长；2.平均面积/线长
		double cal_objective(int target_area, double target_wirelength, double alpha, double beta, Config::LevelObjNorm method) {
			double norm_area, norm_wirelength;
			switch (method) {
			case Config::LevelObjNorm::Average:
				norm_area = accumulate(_his_area.begin(), _his_area.end(), 0.0) / _his_area.size();
				norm_wirelength = accumulate(_his_wirelength.begin(), _his_wirelength.end(), 0.0) / _his_wirelength.size();
				break;
			case Config::LevelObjNorm::Minimum:
				norm_area = *min_element(_his_area.begin(), _his_area.end());
				norm_wirelength = *min_element(_his_wirelength.begin(), _his_wirelength.end());
				break;
			case Config::LevelObjNorm::NoNorm: // 仅在纯优化面积时使用
				norm_area = 1;
				norm_wirelength = 1;
				break;
			default:
				assert(false);
				break;
			}
			return alpha * target_area / norm_area + beta * target_wirelength / norm_wirelength;
		}

		/// 检查并更新最优解
		void check_rule_sol(const SortRule &rule, int target_area, double target_wirelength, const vector<Rect> &target_dst) {
			if (rule.target_objective < _objective) {
				_objective = rule.target_objective;
				_best_area = target_area;
				_best_fill_ratio = 1.0*usedSurfaceArea / _best_area;
				_best_wirelength = target_wirelength;
				_dst = target_dst;
			}
		}

	private:
		const Instance &_ins;
		const vector<Rect> &_src;
		vector<Rect> _dst;

		/// 优化目标相关
		int _best_area;
		double _best_fill_ratio;
		double _best_wirelength;
		double _objective;
		/// 记录最近的n个解，用于归一化
		deque<int> _his_area;
		deque<double> _his_wirelength;

		vector<SortRule> _sort_rules; // 排序规则列表，用于随机局部搜索
		list<int> _rects;			  // SortRule的sequence，相当于指针，使用list快速删除，放置完毕为空
		discrete_distribution<> _discrete_dist;   // 离散概率分布，用于挑选规则(即挑选sequence赋给_rects)
		uniform_int_distribution<> _uniform_dist; // 均匀分布，用于交换矩形顺序

		/// 分组信息，用于挑选候选矩形 `get_candidate_rects()`
		const vector<vector<bool>> &_group_neighbors; // `_group_neighbors[i][j]` 分组_i和_j是否为邻居
		vector<Boundary> _group_boundaries;           // `_group_boundaries[i]`   分组_i的边界

		default_random_engine &_gen;

		// 布局可视化，仅debug
		debug_run(vector<utils_visualize::point_t> _terminal_points;);
		debug_run(vector<utils_visualize::box_t> _dst_boxes;);
		debug_run(vector<utils_visualize::box_t> _wire_boxes;);
	};
}
