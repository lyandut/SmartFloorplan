/** @file SkylineBinPack.h
	@author Jukka Jyl鋘ki

	@brief Implements different bin packer algorithms that use the SKYLINE data structure.

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <limits>
#include <cassert>
#include <cstring>
#include <cmath>

#include "Rect.hpp"
#include "GuillotineBinPack.hpp"

namespace rbp {

	using namespace std;

	/// Represents a single level (a horizontal line) of the skyline/horizon/envelope.
	struct SkylineNode {
		/// The starting x-coordinate (leftmost).
		int x;
		/// The y-coordinate of the skyline level line.
		int y;
		/// The line width. The ending coordinate (inclusive) will be x+width-1.
		int width;
	};

	/// Implements bin packing algorithms that use the SKYLINE data structure to store the bin contents.
	///	Uses GuillotineBinPack as the waste map.
	class SkylineBinPack {
	public:
		/// Instantiates a bin of size (0,0). Call Init to create a new bin.
		SkylineBinPack() : binWidth(0), binHeight(0) {}

		/// Instantiates a bin of the given size.
		SkylineBinPack(int width, int height, bool useWasteMap_) { Init(width, height, useWasteMap_); }

		/// (重新)初始化packer
		void Init(int width, int height, bool useWasteMap_) {
			binWidth = width;
			binHeight = height;
			useWasteMap = useWasteMap_;
			if (useWasteMap) {
				wasteMap.Init(binWidth, binHeight);
				wasteMap.GetFreeRectangles().clear();
			}
			usedSurfaceArea = 0;
			skyLine.clear();
			skyLine.push_back({ 0, 0, binWidth });
			debug_run(disjointRects.Clear());
		}

		/// 定义启发式打包规则.
		enum LevelChoiceHeuristic {
			LevelBottomLeft,   // 最小高度
			LevelMinWasteFit   // 最小浪费
		};

		/// 将一批矩形放置到bin中，矩形可能旋转
		/// @param rects 源矩形列表
		/// @param dst 目的矩形列表
		/// @param method 打包规则.
		void Insert(std::vector<RectSize> &rects, std::vector<Rect> &dst, LevelChoiceHeuristic method) {
			dst.clear();
			while (!rects.empty()) { // 每次迭代放置放置一个矩形
				Rect bestNode;
				int bestScore1 = std::numeric_limits<int>::max();
				int bestScore2 = std::numeric_limits<int>::max();
				int bestSkylineIndex = -1;
				int bestRectIndex = -1;

				for (int i = 0; i < rects.size(); ++i) {
					Rect newNode;
					int score1, score2, skylineIndex;
					switch (method) {
					case LevelChoiceHeuristic::LevelBottomLeft:
						newNode = FindPositionForNewNodeBottomLeft(rects[i].width, rects[i].height, score1, score2, skylineIndex);
						debug_assert(disjointRects.Disjoint(newNode));
						break;
					case LevelChoiceHeuristic::LevelMinWasteFit:
						newNode = FindPositionForNewNodeMinWaste(rects[i].width, rects[i].height, score1, score2, skylineIndex);
						debug_assert(disjointRects.Disjoint(newNode));
						break;
					default:
						assert(false);
						break;
					}
					if (newNode.height != 0) {
						if (score1 < bestScore1 || (score1 == bestScore1 && score2 < bestScore2)) {
							bestNode = newNode;
							bestScore1 = score1;
							bestScore2 = score2;
							bestSkylineIndex = skylineIndex;
							bestRectIndex = i;
						}
					}
				}

				// 没有矩形能放下
				if (bestRectIndex == -1) return;

				// 执行放置
				debug_assert(disjointRects.Disjoint(bestNode));
				debug_run(disjointRects.Add(bestNode));
				AddSkylineLevel(bestSkylineIndex, bestNode); // 更新skyline
				usedSurfaceArea += rects[bestRectIndex].width * rects[bestRectIndex].height;
				rects.erase(rects.begin() + bestRectIndex);
				dst.push_back(bestNode);
			}
		}

		/// 将单个矩形放置到bin中，矩形可能旋转
		Rect Insert(int width, int height, LevelChoiceHeuristic method) {
			// 首先在waste map中尝试放置
			Rect node = wasteMap.Insert(width, height, true, GuillotineBinPack::RectBestShortSideFit, GuillotineBinPack::SplitMaximizeArea);

			// waste map中可以放置
			if (node.height != 0) {
				debug_assert(disjointRects.Disjoint(node));
				debug_run(disjointRects.Add(node));
				usedSurfaceArea += width * height; // 无需更新skyline
				return node;
			}

			// waste map中不能放置
			switch (method) {
			case LevelChoiceHeuristic::LevelBottomLeft:
				return InsertMinHeight(width, height);
			case LevelChoiceHeuristic::LevelMinWasteFit:
				return InsertMinWaste(width, height);
			default:
				assert(false);
				return node;
			}
		}

		/// 计算利用率.
		float Occupancy() const { return (float)usedSurfaceArea / (binWidth * binHeight); }


	protected:
		/// 基于最小高度，遍历skyline选择最佳位置
		Rect FindPositionForNewNodeBottomLeft(int width, int height, int &bestHeight, int &bestWidth, int &bestIndex) const {
			bestHeight = std::numeric_limits<int>::max();
			bestWidth = std::numeric_limits<int>::max(); // 高度相同选择skyline宽度较小的
			bestIndex = -1;
			Rect newNode;
			memset(&newNode, 0, sizeof(newNode));
			for (int i = 0; i < skyLine.size(); ++i) {
				int y; // 矩形放置纵坐标，引用传递
				for (int rotate = 0; rotate <= 1; ++rotate) {
					if (rotate) { swap(width, height); } // 旋转
					if (RectangleFits(i, width, height, y)) {
						if (y + height < bestHeight || (y + height == bestHeight && skyLine[i].width < bestWidth)) {
							bestHeight = y + height;
							bestWidth = skyLine[i].width;
							bestIndex = i;
							newNode.x = skyLine[i].x;
							newNode.y = y;
							newNode.width = width;
							newNode.height = height;
							debug_assert(disjointRects.Disjoint(newNode));
						}
					}
				}
			}
			return newNode;
		}

		/// 基于“最小浪费”，遍历skyline选择最佳位置
		Rect FindPositionForNewNodeMinWaste(int width, int height, int &bestWastedArea, int &bestHeight, int &bestIndex) const {
			bestWastedArea = std::numeric_limits<int>::max();
			bestHeight = std::numeric_limits<int>::max(); // 浪费面积相同选择高度较小的
			bestIndex = -1;
			Rect newNode;
			memset(&newNode, 0, sizeof(newNode));
			for (int i = 0; i < skyLine.size(); ++i) {
				int y, wastedArea; // 引用传递
				for (int rotate = 0; rotate <= 1; ++rotate) {
					if (rotate) { swap(width, height); } // 旋转
					if (RectangleFits(i, width, height, y, wastedArea)) {
						if (wastedArea < bestWastedArea || (wastedArea == bestWastedArea && y + height < bestHeight)) {
							bestWastedArea = wastedArea;
							bestHeight = y + height;
							bestIndex = i;
							newNode.x = skyLine[i].x;
							newNode.y = y;
							newNode.width = width;
							newNode.height = height;
							debug_assert(disjointRects.Disjoint(newNode));
						}
					}
				}
			}
			return newNode;
		}

		/// 判断矩形能否在当前skylineIndex处放置，返回矩形放置点纵坐标y
		bool RectangleFits(int skylineIndex, int width, int height, int &y) const {
			int i = skylineIndex;
			int x = skyLine[i].x;
			if (x + width > binWidth) { return false; }
			y = skyLine[i].y;
			int widthLeft = width;
			while (widthLeft > 0) {
				y = max(y, skyLine[i].y); // 当rect.width > skyline.width，需要检查y坐标
				if (y + height > binHeight) { return false; }
				widthLeft -= skyLine[i].width;
				++i;
				assert(i < skyLine.size() || widthLeft <= 0);
			}
			return true;
		}

		/// 同上，同时返回废料面积（局部）
		bool RectangleFits(int skylineIndex, int width, int height, int &y, int &wastedArea) const {
			bool fits = RectangleFits(skylineIndex, width, height, y);
			if (fits) { wastedArea = ComputeWastedArea(skylineIndex, width, height, y); }
			return fits;
		}

		/// 计算废料面积（局部）
		int ComputeWastedArea(int skylineIndex, int width, int height, int y) const {
			int wastedArea = 0;
			const int rectLeft = skyLine[skylineIndex].x;
			const int rectRight = rectLeft + width;
			for (int i = skylineIndex; i < skyLine.size() && skyLine[i].x < rectRight; ++i) {
				if (skyLine[i].x >= rectRight || skyLine[i].x + skyLine[i].width <= rectLeft) { break; }

				int leftSide = skyLine[i].x;
				int rightSide = min(rectRight, skyLine[i].x + skyLine[i].width);
				assert(y >= skyLine[i].y);
				wastedArea += (rightSide - leftSide) * (y - skyLine[i].y);
			}
			return wastedArea;
		}

		/// 更新skyline，仅支持靠skyline左侧放置
		void AddSkylineLevel(int skylineIndex, const Rect &rect) {
			if (useWasteMap) { AddWasteMapArea(skylineIndex, rect.width, rect.height, rect.y); }

			SkylineNode newNode{ rect.x, rect.y + rect.height, rect.width };
			assert(newNode.x + newNode.width <= binWidth);
			assert(newNode.y <= binHeight);
			skyLine.insert(skyLine.begin() + skylineIndex, newNode);

			for (int i = skylineIndex + 1; i < skyLine.size(); ++i) {
				assert(skyLine[i - 1].x <= skyLine[i].x);
				if (skyLine[i].x < skyLine[i - 1].x + skyLine[i - 1].width) {
					int shrink = skyLine[i - 1].x + skyLine[i - 1].width - skyLine[i].x;
					skyLine[i].x += shrink;
					skyLine[i].width -= shrink;

					if (skyLine[i].width <= 0) {
						skyLine.erase(skyLine.begin() + i);
						--i;
					}
					else { break; }
				}
				else { break; }
			}
			MergeSkylines();
		}

		/// 更新waste map
		void AddWasteMapArea(int skylineIndex, int width, int height, int y) {
			const int rectLeft = skyLine[skylineIndex].x;
			const int rectRight = rectLeft + width;
			for (int i = skylineIndex; i < skyLine.size() && skyLine[i].x < rectRight; ++i) {
				if (skyLine[i].x >= rectRight || skyLine[i].x + skyLine[i].width <= rectLeft) { break; }

				int leftSide = skyLine[i].x;
				int rightSide = min(rectRight, leftSide + skyLine[i].width);
				assert(y >= skyLine[i].y);
				Rect waste{ leftSide, skyLine[i].y, rightSide - leftSide, y - skyLine[i].y };
				debug_assert(disjointRects.Disjoint(waste));
				wasteMap.GetFreeRectangles().push_back(waste);
			}
		}

		/// 合并同一level的skyline节点.
		void MergeSkylines() {
			for (int i = 0; i < skyLine.size() - 1; ++i) {
				if (skyLine[i].y == skyLine[i + 1].y) {
					skyLine[i].width += skyLine[i + 1].width;
					skyLine.erase(skyLine.begin() + i + 1);
					--i;
				}
			}
		}

		Rect InsertMinHeight(int width, int height) {
			int bestHeight, bestWidth, bestIndex;
			Rect newNode = FindPositionForNewNodeBottomLeft(width, height, bestHeight, bestWidth, bestIndex);
			if (bestIndex != -1) {
				debug_assert(disjointRects.Disjoint(newNode));
				debug_run(disjointRects.Add(newNode));
				AddSkylineLevel(bestIndex, newNode);
				usedSurfaceArea += width * height;
			}
			else {
				memset(&newNode, 0, sizeof(newNode));
			}
			return newNode;
		}

		Rect InsertMinWaste(int width, int height) {
			int bestWastedArea, bestHeight, bestIndex;
			Rect newNode = FindPositionForNewNodeMinWaste(width, height, bestWastedArea, bestHeight, bestIndex);
			if (bestIndex != -1) {
				debug_assert(disjointRects.Disjoint(newNode));
				debug_run(disjointRects.Add(newNode));
				AddSkylineLevel(bestIndex, newNode);
				usedSurfaceArea += width * height;
			}
			else {
				memset(&newNode, 0, sizeof(newNode));
			}
			return newNode;
		}


	protected:
		int binWidth;
		int binHeight;

		std::vector<SkylineNode> skyLine;

		/// 已使用面积
		unsigned long usedSurfaceArea;

		/// If true, we use the GuillotineBinPack structure to recover wasted areas into a waste map.
		bool useWasteMap;
		GuillotineBinPack wasteMap;

		debug_run(DisjointRectCollection disjointRects;);
	};

}