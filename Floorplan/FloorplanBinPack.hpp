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

namespace fbp {

	using namespace rbp;

	class FloorplanBinPack final : public SkylineBinPack {

	public:
		FloorplanBinPack() = delete;

		FloorplanBinPack(const vector<Rect> &src, const vector<vector<bool>> &group_neighbors,
			vector<Boundary> &group_boundaries, int bin_width, int bin_height) :
			_src(src), _group_neighbors(group_neighbors), _group_boundaries(group_boundaries),
			SkylineBinPack(bin_width, bin_height, false) {

			_min_bin_height = numeric_limits<int>::max();

			_seed = random_device{}(); // set random seed

			init_sort_rules();
		}

		const vector<Rect>& get_dst() const { return _dst; }

		float get_fill_ratio() const { return _best_fill_ratio; }

		/// 固定宽度，更新高度
		void update_bin_height(int new_height, const vector<Boundary> &new_boundaries) {
			binHeight = new_height;
			_group_boundaries = new_boundaries;
		}

		/// 分组搜索策略
		enum LevelGroupSearch {
			LevelSelfishly,       // 仅当前分组
			LevelNeighborAll,     // 当前分组和全部邻居分组
			LevelNeighborPartial  // [todo] 当前分组和左/下邻居分组，或可考虑按百分比选取一部分矩形
		};

		/// 放置策略
		enum LevelHeuristicSearch {
			LevelMinHeightFit,    // 最小高度rectangle，不考虑靠skyline右侧放置
			LevelMinWasteFit,     // 最小浪费rectangle，不考虑靠skyline右侧放置
			LevelBottomLeftScore  // 最下最左skyline，考虑靠skyline右侧放置
		};

		/// 基于binWidth进行随机局部搜索，返回最小高度
		void random_local_search(int iter, LevelGroupSearch method) {
			// the first time to call RLS on W_k
			if (iter == 1) {
				for (auto &rule : _sort_rules) {
					_rects.assign(rule.sequence.begin(), rule.sequence.end());
					vector<Rect> target_dst;
					rule.target_height = insert_bottom_left_score(target_dst, method);
					if (rule.target_height < _min_bin_height) {
						_min_bin_height = rule.target_height;
						_best_fill_ratio = (float)usedSurfaceArea / (binWidth * _min_bin_height);
						_dst = target_dst;
					}
				}
				sort(_sort_rules.begin(), _sort_rules.end(), [](auto &lhs, auto &rhs) { return lhs.target_height > rhs.target_height; });
			}

			// 构造初始解
			static default_random_engine gen(_seed);
			vector<int> probs;
			probs.reserve(_sort_rules.size());
			for (int i = 1; i <= _sort_rules.size(); ++i) { probs.push_back(2 * i); }
			discrete_distribution<> discrete_dist(probs.begin(), probs.end());

			SortRule &picked_rule = _sort_rules[discrete_dist(gen)];
			_rects.assign(picked_rule.sequence.begin(), picked_rule.sequence.end());
			vector<Rect> target_dst;
			picked_rule.target_height = min(picked_rule.target_height, insert_bottom_left_score(target_dst, method));
			if (picked_rule.target_height < _min_bin_height) {
				_min_bin_height = picked_rule.target_height;
				_best_fill_ratio = (float)usedSurfaceArea / (binWidth * _min_bin_height);
				_dst = target_dst;
			}

			// 迭代优化
			for (int i = 1; i <= iter; ++i) {
				SortRule new_rule = picked_rule;
				uniform_int_distribution<> uniform_dist(0, new_rule.sequence.size() - 1);
				int a = uniform_dist(gen);
				int b = uniform_dist(gen);
				while (a == b) { b = uniform_dist(gen); }
				swap(new_rule.sequence[a], new_rule.sequence[b]);
				_rects.assign(new_rule.sequence.begin(), new_rule.sequence.end());
				vector<Rect> target_dst;
				new_rule.target_height = min(new_rule.target_height, insert_bottom_left_score(target_dst, method));
				if (new_rule.target_height < picked_rule.target_height) {
					picked_rule = new_rule;
				}
				if (picked_rule.target_height < _min_bin_height) {
					_min_bin_height = picked_rule.target_height;
					_best_fill_ratio = (float)usedSurfaceArea / (binWidth * _min_bin_height);
					_dst = target_dst;
				}
			}

			// 更新排序规则列表
			sort(_sort_rules.begin(), _sort_rules.end(), [](auto &lhs, auto &rhs) { return lhs.target_height > rhs.target_height; });
		}

		int insert_bottom_left_score(vector<Rect> &dst, LevelGroupSearch method) {
			Init(binWidth, binHeight, usedSurfaceArea); // 每次调用重置usedSurfaceArea和skyLine
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

				// 没有矩形能放下
				if (best_rect_index == -1) { return binHeight + 1; }

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

			return min_bin_height;
		}

		/// [deprecated]
		int insert_greedy_fit(vector<Rect> &dst, LevelHeuristicSearch method) {
			Init(binWidth, binHeight, usedSurfaceArea); // 每次调用重置usedSurfaceArea和skyLine
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
					case LevelHeuristicSearch::LevelMinHeightFit:
						new_node = find_rect_for_skyline_min_height(i, _rects, score_1, score_2, rect_index);
						debug_assert(disjointRects.Disjoint(new_node));
						break;
					case LevelHeuristicSearch::LevelMinWasteFit:
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

				// 没有矩形能放下
				if (best_rect_index == -1) { return binHeight + 1; }

				// 执行放置
				debug_assert(disjointRects.Disjoint(best_node));
				debug_run(disjointRects.Add(best_node));
				AddSkylineLevel(best_skyline_index, best_node);
				usedSurfaceArea += _src.at(best_rect_index).width * _src.at(best_rect_index).height;
				_rects.remove(best_rect_index);
				dst.push_back(best_node);
				min_bin_height = max(min_bin_height, best_node.y + best_node.height);
			}

			return min_bin_height;
		}

	private:
		/// 初始化排序规则列表
		void init_sort_rules() {
			vector<int> seq(_src.size());
			// 0_输入顺序
			iota(seq.begin(), seq.end(), 0);
			for (int i = 0; i < 5; ++i) { _sort_rules.push_back({ seq, numeric_limits<int>::max() }); }
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
			shuffle(_sort_rules[4].sequence.begin(), _sort_rules[4].sequence.end(), mt19937{ _seed });

			// 默认赋给_rect输入顺序，用于测试贪心算法
			_rects.assign(_sort_rules[0].sequence.begin(), _sort_rules[0].sequence.end());
		}

		/// 基于分组策略挑选候选矩形，减小搜索规模
		list<int> get_candidate_rects(const SkylineNode &skyline, LevelGroupSearch method) {
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
			case LevelGroupSearch::LevelSelfishly:
				for (int r : _rects) { // 必须顺序遍历_rects才能使得sort_rule生效
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
				}
				break;
			case LevelGroupSearch::LevelNeighborAll:
				for (int r : _rects) {
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
					else if (_group_neighbors[_src.at(r).gid][gid]) {
						candidate_rects.push_back(r);
					}
				}
				break;
			case LevelGroupSearch::LevelNeighborPartial:
				for (int r : _rects) {
					if (_src.at(r).gid == gid) { candidate_rects.push_back(r); }
					else if (_group_neighbors[_src.at(r).gid][gid] && _src.at(r).gid < gid) {
						candidate_rects.push_back(r);
					}
				}
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

	private:
		const vector<Rect> &_src;
		vector<Rect> _dst;
		int _min_bin_height;
		float _best_fill_ratio;

		/// 排序规则定义
		struct SortRule {
			vector<int> sequence;
			int target_height;
		};
		vector<SortRule> _sort_rules; // 排序规则列表，用于随机局部搜索
		list<int> _rects;			  // SortRule的sequence，相当于指针，使用list快速删除，放置完毕为空

		/// 分组信息
		const vector<vector<bool>> &_group_neighbors; // `_group_neighbors[i][j]` 分组_i和_j是否为邻居	public:
		vector<Boundary> &_group_boundaries;          // `_group_boundaries[i]`   分组_i的边界

		unsigned int _seed;
	};
}
