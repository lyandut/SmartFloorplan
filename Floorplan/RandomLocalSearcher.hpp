//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include "FloorplanPacker.hpp"

namespace fbp {

	using namespace std;

	class RandomLocalSearcher : public FloorplanPacker {

		/// 排序规则定义
		struct SortRule {
			vector<int> sequence;
			double target_objective;
		};

	public:

		RandomLocalSearcher() = delete;

		RandomLocalSearcher(const Instance& ins, const vector<Rect>& src, int bin_width, default_random_engine& gen) :
			FloorplanPacker(ins, src, bin_width, gen) {
			reset();
			init_sort_rules();
		}

		/// 基于_bin_width进行随机局部搜索
		void run(int iter, double alpha, double beta, Config::LevelWireLength level_wl, Config::LevelObjDist level_dist) {
			// the first time to call RLS on W_k
			if (iter == 1) {
				for (auto& rule : _sort_rules) {
					_rects.assign(rule.sequence.begin(), rule.sequence.end());
					vector<Rect> target_dst;
					vector<bool> is_packed(_src.size(), true);
					int target_area = insert_bottom_left_score(target_dst) * _bin_width;
					double target_dist;
					double target_wirelength = cal_wirelength(target_dst, is_packed, target_dist, level_wl, level_dist);
					rule.target_objective = cal_objective(target_area, target_dist, alpha, beta);
					update_objective(rule.target_objective, target_area, target_wirelength, target_dst);
				}
				// 降序排列，越后面的目标函数值越小选中概率越大
				sort(_sort_rules.begin(), _sort_rules.end(), [](auto& lhs, auto& rhs) {
					return lhs.target_objective > rhs.target_objective; });
			}

			// 迭代优化
			SortRule& picked_rule = _sort_rules[_discrete_dist(_gen)];
			bool is_resort_needed = false;
			for (int i = 1; i <= iter; ++i) {
				SortRule new_rule = picked_rule;
				if (iter % 4) { swap_sort_rule(new_rule); }
				else { rotate_sort_rule(new_rule); }
				_rects.assign(new_rule.sequence.begin(), new_rule.sequence.end());
				vector<Rect> target_dst;
				vector<bool> is_packed(_src.size(), true);
				int target_area = insert_bottom_left_score(target_dst) * _bin_width;
				double target_dist;
				double target_wirelength = cal_wirelength(target_dst, is_packed, target_dist, level_wl, level_dist);
				new_rule.target_objective = cal_objective(target_area, target_dist, alpha, beta);
				if (new_rule.target_objective <= picked_rule.target_objective) {
					picked_rule = new_rule;
					is_resort_needed = true;
					update_objective(picked_rule.target_objective, target_area, target_wirelength, target_dst);
				}
			}
			// 更新排序规则列表
			if (is_resort_needed) {
				sort(_sort_rules.begin(), _sort_rules.end(), [](auto& lhs, auto& rhs) {
					return lhs.target_objective > rhs.target_objective; });
			}
		}

		/// 基于最下最左和打分策略，贪心构造一个完整解
		int insert_bottom_left_score(vector<Rect>& dst) {
			int skyline_height = 0;
			reset();
			dst = _src;

			while (!_rects.empty()) {
				auto bottom_skyline_iter = min_element(_skyline.begin(), _skyline.end(),
					[](auto& lhs, auto& rhs) { return lhs.y < rhs.y; });
				int best_skyline_index = distance(_skyline.begin(), bottom_skyline_iter);
				int min_rect_width = _src.at(*min_element(_rects.begin(), _rects.end(),
					[this](int lhs, int rhs) { return _src.at(lhs).width < _src.at(rhs).width; })).width;

				if (_skyline[best_skyline_index].width < min_rect_width) { // 最小宽度矩形放不进去，需要填坑
					if (best_skyline_index == 0) { _skyline[best_skyline_index].y = _skyline[best_skyline_index + 1].y; }
					else if (best_skyline_index == _skyline.size() - 1) { _skyline[best_skyline_index].y = _skyline[best_skyline_index - 1].y; }
					else { _skyline[best_skyline_index].y = min(_skyline[best_skyline_index - 1].y, _skyline[best_skyline_index + 1].y); }
					merge_skylines(_skyline);
					continue;
				}

				int best_rect_index = find_rect_for_skyline_bottom_left(best_skyline_index, dst);
				assert(best_rect_index != -1);

				// 从未放置列表中删除
				_rects.remove(best_rect_index);

				// 更新skyline
				SkylineNode new_skyline_node = { dst[best_rect_index].x, dst[best_rect_index].y + dst[best_rect_index].height, dst[best_rect_index].width };
				if (new_skyline_node.x == _skyline[best_skyline_index].x) { // 靠左
					_skyline.insert(_skyline.begin() + best_skyline_index, new_skyline_node);
					_skyline[best_skyline_index + 1].x += new_skyline_node.width;
					_skyline[best_skyline_index + 1].width -= new_skyline_node.width;
					merge_skylines(_skyline);
				}
				else { // 靠右
					_skyline.insert(_skyline.begin() + best_skyline_index + 1, new_skyline_node);
					_skyline[best_skyline_index].width -= new_skyline_node.width;
					merge_skylines(_skyline);
				}
				skyline_height = max(skyline_height, new_skyline_node.y);
			}

			return skyline_height;
		}

	private:
		/// 每次迭代重置_skyLine
		void reset() {
			_skyline.clear();
			_skyline.push_back({ 0,0,_bin_width });
		}

		/// 初始化排序规则列表
		void init_sort_rules() {
			vector<int> seq(_src.size());
			// 0_输入顺序
			iota(seq.begin(), seq.end(), 0);
			_sort_rules.reserve(5);
			for (int i = 0; i < 5; ++i) { _sort_rules.push_back({ seq, numeric_limits<double>::max() }); }
			// 1_面积递减
			sort(_sort_rules[1].sequence.begin(), _sort_rules[1].sequence.end(), [this](int lhs, int rhs) {
				return _ins.get_blocks().at(lhs).area > _ins.get_blocks().at(rhs).area; });
			// 2_高度递减
			sort(_sort_rules[2].sequence.begin(), _sort_rules[2].sequence.end(), [this](int lhs, int rhs) {
				return _src.at(lhs).height > _src.at(rhs).height; });
			// 3_宽度递减
			sort(_sort_rules[3].sequence.begin(), _sort_rules[3].sequence.end(), [this](int lhs, int rhs) {
				return _src.at(lhs).width > _src.at(rhs).width; });
			// 4_随机排序
			shuffle(_sort_rules[4].sequence.begin(), _sort_rules[4].sequence.end(), _gen);

			// `_rects`默认输入顺序
			_rects.assign(_sort_rules[0].sequence.begin(), _sort_rules[0].sequence.end());

			// 离散概率分布初始化
			vector<int> probs; probs.reserve(_sort_rules.size());
			for (int i = 1; i <= _sort_rules.size(); ++i) { probs.push_back(2 * i); }
			_discrete_dist = discrete_distribution<>(probs.begin(), probs.end());
			// 均匀分布初始化
			_uniform_dist = uniform_int_distribution<>(0, _src.size() - 1);
		}

		/// 邻域动作1：交换两个块的顺序
		void swap_sort_rule(SortRule& rule) {
			int a = _uniform_dist(_gen);
			int b = _uniform_dist(_gen);
			while (a == b) { b = _uniform_dist(_gen); }
			swap(rule.sequence[a], rule.sequence[b]);
		}

		/// 邻域动作2：连续多个块移动
		void rotate_sort_rule(SortRule& rule) {
			int a = _uniform_dist(_gen);
			rotate(rule.sequence.begin(), rule.sequence.begin() + a, rule.sequence.end());
		}

		/// 基于打分策略为左下角选一个块
		int find_rect_for_skyline_bottom_left(int skyline_index, vector<Rect>& dst) {
			int best_rect = -1, best_score = -1;
			for (int r : _rects) {
				int width = _src.at(r).width, height = _src.at(r).height;
				int x, score;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					if (rotate) { swap(width, height); }
					if (score_rect_for_skyline_bottom_left(skyline_index, width, height, x, score)) {
						if (best_score < score) {
							best_score = score;
							best_rect = r;
							dst[r].x = x;
							dst[r].y = _skyline[skyline_index].y;
							dst[r].width = width;
							dst[r].height = height;
						}
					}
				}
			}
			// (d)(f)(h)的退化情况
			if ((best_score == 4 || best_score == 2 || best_score == 0) && _rects.size() > 1) {
				int min_unpacked_width = numeric_limits<int>::max();
				for (int r : _rects) {
					if (r == best_rect) { continue; }
					min_unpacked_width = min(min_unpacked_width, _src.at(r).width);
				}
				// 未放置的最小宽度放不下，导致浪费
				if (min_unpacked_width > _skyline[skyline_index].width - dst[best_rect].width) {
					SkylineSpace space = skyline_nodo_to_space(_skyline, skyline_index);
					int min_space_height = min(space.hl, space.hr);
					int new_best_rect = -1;
					int new_best_width = 0, new_best_height;
					for (int r : _rects) {
						for (int rotate = 0; rotate <= 1; ++rotate) {
							int width = _src.at(r).width, height = _src.at(r).height;
							if (rotate) { swap(width, height); }
							if (height >= min_space_height // 高不小于min_space_height
								&& width <= space.width // 宽能放下
								&& width > new_best_width) { // 宽最长
								new_best_rect = r;
								new_best_width = width;
								new_best_height = height;
							}
						}
					}
					if (new_best_rect != -1) {
						best_rect = new_best_rect;
						dst[best_rect].width = new_best_width;
						dst[best_rect].height = new_best_height;
						dst[best_rect].y = _skyline[skyline_index].y;
						dst[best_rect].x = space.hl >= space.hr ? // 必须靠高的一侧放
							_skyline[skyline_index].x : // 靠左
							_skyline[skyline_index].x + _skyline[skyline_index].width - new_best_width; // 靠右
					}
				}
			}

			return best_rect;
		}

		/// 打分策略
		bool score_rect_for_skyline_bottom_left(int skyline_index, int width, int height, int& x, int& score) {
			if (width > _skyline[skyline_index].width) { return false; }

			SkylineSpace space = skyline_nodo_to_space(_skyline, skyline_index);
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

				if (score == 2) { x = _skyline[skyline_index].x + _skyline[skyline_index].width - width; }
				else { x = _skyline[skyline_index].x; }
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

				if (score == 4 || score == 0) { x = _skyline[skyline_index].x + _skyline[skyline_index].width - width; }
				else { x = _skyline[skyline_index].x; }
			}
			if (x + width > _bin_width) { return false; }

			return true;
		}

	private:
		Skyline _skyline;

		// 排序规则列表，用于随机局部搜索  
		vector<SortRule> _sort_rules;
		list<int> _rects; // SortRule的sequence，相当于指针，使用list快速删除，放置完毕为空
		discrete_distribution<> _discrete_dist;   // 离散概率分布，用于挑选规则(即挑选sequence赋给_rects)
		uniform_int_distribution<> _uniform_dist; // 均匀分布，用于交换矩形顺序
	};

}
