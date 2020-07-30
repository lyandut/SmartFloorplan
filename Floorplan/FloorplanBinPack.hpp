//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include "SkylineBinPack.hpp"

namespace fbp {

	using namespace rbp;

	class FloorplanBinPack final : public SkylineBinPack {

	public:
		FloorplanBinPack() = delete;

		FloorplanBinPack(vector<RectSize> &rects, int bin_width, int bin_height, bool use_waste_map = false, int group_num = 1) :
			_rects(rects), SkylineBinPack(bin_width, bin_height, use_waste_map) {
			init_rects_groups(group_num);
		}

		/// 分组搜索规则
		enum LevelGroupSearch {
			LevelSelfishly,      // 仅当前分组
			LevelNeighborAll,    // 当前分组和全部邻居分组
			LevelNeighborPartial // [todo] 当前分组和部分邻居分组，或按百分比选取邻居的部分矩形
		};

		enum LevelHeuristicSearch {
			LevelBottomLeftScore, // 最下最左skyline，考虑靠skyline右侧放置
			LevelMinHeightFit,    // 最小高度rectangle，不考虑右侧放置
			LevelMinWasteFit      // 最小浪费rectangle，不考虑右侧放置
		};

		void insert(vector<Rect> &dst, LevelGroupSearch method_1, LevelHeuristicSearch method_2) {
			dst.clear();
			while (!_rects.empty()) {
				Rect best_node;
				int best_score_1 = numeric_limits<int>::max();
				int best_score_2 = numeric_limits<int>::max();
				int best_skyline_index = -1;
				int best_rect_index = -1;

				if (method_2 == LevelHeuristicSearch::LevelBottomLeftScore) {
					// [todo] 将skyline按照y坐标排序，每次调用min_element复杂度较高
					auto bottom_iter = min_element(skyLine.begin(), skyLine.end(), [](auto &lhs, auto &rhs) { return lhs.y < rhs.y; });
					vector<RectSize> candidate_rects = get_candidate_rects(*bottom_iter, method_1);
					best_node = find_rect_for_skyline_bottom_left(distance(skyLine.begin(), bottom_iter), candidate_rects, best_rect_index);
					best_skyline_index = 
				}
				else { // 搜索所有skyline
					for (int i = 0; i < skyLine.size(); ++i) {
						Rect new_node;
						int score_1, score_2, rect_index;
						vector<RectSize> candidate_rects = get_candidate_rects(skyLine[i], method_1);
						switch (method_2) {
						case LevelHeuristicSearch::LevelMinHeightFit:
							new_node = find_rect_for_skyline_min_height(i, candidate_rects, score_1, score_2, rect_index);
							debug_assert(disjointRects.Disjoint(new_node));
							break;
						case LevelHeuristicSearch::LevelMinWasteFit:
							new_node = find_rect_for_skyline_min_waste(i, candidate_rects, score_1, score_2, rect_index);
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
				}

				if (best_node.height == 0 || best_rect_index == -1) { return; } // 无解

				// 执行放置
				debug_assert(disjointRects.Disjoint(best_node));
				debug_run(disjointRects.Add(best_node));
				AddSkylineLevel(best_skyline_index, best_node);
				usedSurfaceArea += _rects[best_rect_index].width * _rects[best_rect_index].height;
				_rects.erase(_rects.begin() + best_rect_index);
				_group_rects[_rects[best_rect_index].group_id].erase(); // [todo] 删除分组中的矩形，改为指针？
				dst.push_back(best_node);
			}
		}

	private:
		/// 调用QAP将矩形分组，初始化分组信息
		void init_rects_groups(int group_num) {
			_group_rects.resize(group_num);
			// [todo] 
			_group_rects[0] = _rects;
		}

		/// 基于分组策略生成候选矩形，减小搜索规模
		vector<RectSize> get_candidate_rects(const SkylineNode &skyline, LevelGroupSearch method) {
			int group_id = 0;
			for (; group_id < _group_boundaries.size(); ++group_id) {
				if (skyline.x >= _group_boundaries[group_id].x && skyline.y >= _group_boundaries[group_id].y
					&& skyline.x < _group_boundaries[group_id].x + _group_boundaries[group_id].width
					&& skyline.y < _group_boundaries[group_id].y + _group_boundaries[group_id].height) {
					break;
				}
			}
			vector<RectSize> candidate_rects(_group_rects[group_id]);
			switch (method) {
			case LevelGroupSearch::LevelNeighborAll:
				for (int id : _group_neighbors[group_id]) {
					candidate_rects.insert(candidate_rects.end(), _group_rects[id].begin(), _group_rects[id].end());
				}
				break;
			case LevelGroupSearch::LevelNeighborPartial:
				for (int id : _group_neighbors[group_id]) {
					if (id < group_id) { // [todo] 仅考虑左/下的分组
						candidate_rects.insert(candidate_rects.end(), _group_rects[id].begin(), _group_rects[id].end());
					}
				}
				break;
			case LevelGroupSearch::LevelSelfishly:
				break;
			default:
				assert(false);
				break;
			}
			return candidate_rects;
		}

		/// 基于最小高度，为当前skyline选择最佳放置的矩形
		Rect find_rect_for_skyline_min_height(int skyline_index, const vector<RectSize> &rects,
			int &best_height, int &best_width, int &best_index) {
			best_height = numeric_limits<int>::max();
			best_width = numeric_limits<int>::max(); // 高度相同选择宽度较小的,[todo]可能选择宽度较大的比较好？
			best_index = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node)); // 确保无解时返回高度为0
			for (auto &rect : rects) {
				int y;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					int width = rect.width, height = rect.height;
					if (rotate) { swap(width, height); }
					if (RectangleFits(skyline_index, width, height, y)) {
						if (y + height < best_height || (y + height == best_height && width < best_width)) {
							best_height = y + height;
							best_width = width;
							best_index = rect.id;
							new_node = { skyLine[skyline_index].x, y, width, height };
							debug_assert(disjointRects.Disjoint(new_node));
						}
					}
				}
			}
			return new_node;
		}

		/// 基于最小浪费
		Rect find_rect_for_skyline_min_waste(int skyline_index, const vector<RectSize> &rects,
			int &best_wasted_area, int &best_height, int &best_index) {
			best_wasted_area = numeric_limits<int>::max();
			best_height = numeric_limits<int>::max();
			best_index = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node));
			for (auto &rect : rects) {
				int y, wasted_area;
				for (int rotate = 0; rotate <= 1; ++rotate) {
					int width = rect.width, height = rect.height;
					if (rotate) { swap(width, height); }
					if (RectangleFits(skyline_index, width, height, y, wasted_area)) {
						if (wasted_area < best_wasted_area || (wasted_area == best_wasted_area && y + height < best_height)) {
							best_wasted_area = wasted_area;
							best_height = y + height;
							best_index = rect.id;
							new_node = { skyLine[skyline_index].x, y, width, height };
							debug_assert(disjointRects.Disjoint(new_node));
						}
					}
				}
			}
			return new_node;
		}

		/// 论文idea，基于最左/最下和打分策略
		/// [todo] 未实现(d)(f)->(h)的退化情况
		// [todo] 这里实现错误，skyline可能变化，应该将填坑部分代码向上提出一层
		Rect find_rect_for_skyline_bottom_left(int skyline_index, const vector<RectSize> &rects, int &best_index) {
			best_index = -1;
			int best_socre = -1;
			Rect new_node;
			memset(&new_node, 0, sizeof(new_node));

			while (best_socre == -1) {
				for (auto &rect : rects) {
					int x, score;
					for (int rotate = 0; rotate <= 1; ++rotate) {
						int width = rect.width, height = rect.height;
						if (rotate) { swap(width, height); }
						if (score_rect_for_skyline_bottom_left(skyline_index, width, height, x, score)) {
							if (best_socre < score) {
								best_socre = score;
								best_index = rect.id;
								new_node = { x, skyLine[skyline_index].y, width, height };
								debug_assert(disjointRects.Disjoint(new_node));
							}
						}
					}
				}
				if (best_socre == -1) { // 填坑
					if (skyline_index == 0) { skyLine[skyline_index].y = skyLine[skyline_index + 1].y; }
					else if (skyline_index == skyLine.size() - 1) { skyLine[skyline_index].y = skyLine[skyline_index - 1].y; }
					else { skyLine[skyline_index].y = min(skyLine[skyline_index - 1].y, skyLine[skyline_index + 1].y); }
					MergeSkylines();
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

		bool score_rect_for_skyline_bottom_left(int skyline_index, int width, int height, int &x, int &score) {
			if (width > skyLine.front().width) { return false; }

			SkylineSpace space = skyline_nodo_to_space(skyline_index);
			if (space.hl >= space.hr) {
				if (width == space.width && height == space.hl) { score = 7; }
				else if (width == space.width && height == space.hr) { score = 6; }
				else if (width == space.width && height > space.hl) { score = 5; }
				else if (width < space.width && height == space.hl) { score = 4; }
				else if (width == space.width && height < space.hl && height > space.hr) { score = 3; }
				else if (width < space.width && height == space.hr) { score = 2; } // 靠右
				else if (width == space.width && height < space.hr) { score = 1; }
				else if (width < space.width && height > space.hl) { score = 0; }
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
				else if (width < space.width && height > space.hr) { score = 0; } // 靠右
				else { return false; }

				if (score == 4 || score == 0) { x = skyLine[skyline_index].x + skyLine[skyline_index].width - width; }
				else { x = skyLine[skyline_index].x; }
			}
			return true;
		}

		SkylineSpace skyline_nodo_to_space(int skyline_index) {
			int hl, hr;
			if (skyline_index == 0) {
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
		/// 全部待放置矩形，放置完毕为空
		vector<RectSize> _rects;

		/// 分组信息
		vector<vector<RectSize>> _group_rects;
		vector<vector<int>> _group_neighbors;
		vector<Rect> _group_boundaries;
	};
}




