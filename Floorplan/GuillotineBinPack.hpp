/** @file GuillotineBinPack.h
	@author Jukka Jylänki

	@brief Implements different bin packer algorithms that use the GUILLOTINE data structure.

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

namespace rbp {

	using namespace std;

	/** GuillotineBinPack implements different variants of bin packer algorithms that use the GUILLOTINE data structure
		to keep track of the free space of the bin where rectangles may be placed. */
	class GuillotineBinPack {
	public:
		/// The initial bin size will be (0,0). Call Init to set the bin size.
		GuillotineBinPack() : binWidth(0), binHeight(0) {}

		/// Initializes a new bin of the given size.
		GuillotineBinPack(int width, int height) { Init(width, height); }

		/// (Re)initializes the packer to an empty bin of width x height units. Call whenever
		/// you need to restart with a new bin.
		void Init(int width, int height) {
			binWidth = width;
			binHeight = height;
			usedRectangles.clear();
			freeRectangles.clear();
			freeRectangles.push_back({ 0, 0, binWidth, binHeight });
			debug_run(disjointRects.Clear());
		}

		/// Specifies the different choice heuristics that can be used when deciding which of the free subrectangles
		/// to place the to-be-packed rectangle into.
		enum FreeRectChoiceHeuristic
		{
			RectBestAreaFit, ///< -BAF
			RectBestShortSideFit, ///< -BSSF
			RectBestLongSideFit, ///< -BLSF
			RectWorstAreaFit, ///< -WAF
			RectWorstShortSideFit, ///< -WSSF
			RectWorstLongSideFit ///< -WLSF
		};

		/// Specifies the different choice heuristics that can be used when the packer needs to decide whether to
		/// subdivide the remaining free space in horizontal or vertical direction.
		enum GuillotineSplitHeuristic
		{
			SplitShorterLeftoverAxis, ///< -SLAS
			SplitLongerLeftoverAxis, ///< -LLAS
			SplitMinimizeArea, ///< -MINAS, Try to make a single big rectangle at the expense of making the other small.
			SplitMaximizeArea, ///< -MAXAS, Try to make both remaining rectangles as even-sized as possible.
			SplitShorterAxis, ///< -SAS
			SplitLongerAxis ///< -LAS
		};

		/// Inserts a single rectangle into the bin. The packer might rotate the rectangle, in which case the returned
		/// struct will have the width and height values swapped.
		/// @param merge If true, performs free Rectangle Merge procedure after packing the new rectangle. This procedure
		///		tries to defragment the list of disjoint free rectangles to improve packing performance, but also takes up 
		///		some extra time.
		/// @param rectChoice The free rectangle choice heuristic rule to use.
		/// @param splitMethod The free rectangle split heuristic rule to use.
		Rect Insert(int width, int height, bool merge, FreeRectChoiceHeuristic rectChoice, GuillotineSplitHeuristic splitMethod) {
			// Find where to put the new rectangle.
			int freeNodeIndex = 0;
			Rect newRect = FindPositionForNewNode(width, height, rectChoice, &freeNodeIndex);

			// Abort if we didn't have enough space in the bin.
			if (newRect.height == 0)
				return newRect;

			// Remove the space that was just consumed by the new rectangle.
			SplitFreeRectByHeuristic(freeRectangles[freeNodeIndex], newRect, splitMethod);
			freeRectangles.erase(freeRectangles.begin() + freeNodeIndex);

			// Perform a Rectangle Merge step if desired.
			if (merge)
				MergeFreeList();

			// Remember the new used rectangle.
			usedRectangles.push_back(newRect);

			// Check that we're really producing correct packings here.
			debug_assert(disjointRects.Add(newRect) == true);

			return newRect;
		}

		/// Computes the ratio of used/total surface area. 0.00 means no space is yet used, 1.00 means the whole bin is used.
		float Occupancy() const {
			///\todo The occupancy rate could be cached/tracked incrementally instead
			///      of looping through the list of packed rectangles here.
			unsigned long usedSurfaceArea = 0;
			for (size_t i = 0; i < usedRectangles.size(); ++i)
				usedSurfaceArea += usedRectangles[i].width * usedRectangles[i].height;

			return (float)usedSurfaceArea / (binWidth * binHeight);
		}

		/// Returns the internal list of disjoint rectangles that track the free area of the bin.
		/// You may alter this vector any way desired, as long as the end result still is a list of disjoint rectangles.
		std::vector<Rect> &GetFreeRectangles() { return freeRectangles; }

		/// Returns the list of packed rectangles. You may alter this vector at will, for example, you can move a Rect from
		/// this list to the Free Rectangles list to free up space on-the-fly, but notice that this causes fragmentation.
		std::vector<Rect> &GetUsedRectangles() { return usedRectangles; }

		/// Performs a Rectangle Merge operation. This procedure looks for adjacent free rectangles and merges them if they
		/// can be represented with a single rectangle. Takes up Theta(|freeRectangles|^2) time.
		void MergeFreeList() {
			debug_run(DisjointRectCollection test);
			for (size_t i = 0; i < freeRectangles.size(); ++i) debug_assert(test.Add(freeRectangles[i]) == true);

			// Do a Theta(n^2) loop to see if any pair of free rectangles could me merged into one.
			// Note that we miss any opportunities to merge three rectangles into one. (should call this function again to detect that)
			for (size_t i = 0; i < freeRectangles.size(); ++i)
				for (size_t j = i + 1; j < freeRectangles.size(); ++j)
				{
					if (freeRectangles[i].width == freeRectangles[j].width && freeRectangles[i].x == freeRectangles[j].x)
					{
						if (freeRectangles[i].y == freeRectangles[j].y + freeRectangles[j].height)
						{
							freeRectangles[i].y -= freeRectangles[j].height;
							freeRectangles[i].height += freeRectangles[j].height;
							freeRectangles.erase(freeRectangles.begin() + j);
							--j;
						}
						else if (freeRectangles[i].y + freeRectangles[i].height == freeRectangles[j].y)
						{
							freeRectangles[i].height += freeRectangles[j].height;
							freeRectangles.erase(freeRectangles.begin() + j);
							--j;
						}
					}
					else if (freeRectangles[i].height == freeRectangles[j].height && freeRectangles[i].y == freeRectangles[j].y)
					{
						if (freeRectangles[i].x == freeRectangles[j].x + freeRectangles[j].width)
						{
							freeRectangles[i].x -= freeRectangles[j].width;
							freeRectangles[i].width += freeRectangles[j].width;
							freeRectangles.erase(freeRectangles.begin() + j);
							--j;
						}
						else if (freeRectangles[i].x + freeRectangles[i].width == freeRectangles[j].x)
						{
							freeRectangles[i].width += freeRectangles[j].width;
							freeRectangles.erase(freeRectangles.begin() + j);
							--j;
						}
					}
				}

			debug_run(test.Clear());
			for (size_t i = 0; i < freeRectangles.size(); ++i) debug_assert(test.Add(freeRectangles[i]) == true);
		}

	private:
		int binWidth;
		int binHeight;

		/// Stores a list of all the rectangles that we have packed so far. This is used only to compute the Occupancy ratio,
		/// so if you want to have the packer consume less memory, this can be removed.
		std::vector<Rect> usedRectangles;

		/// Stores a list of rectangles that represents the free area of the bin. This rectangles in this list are disjoint.
		std::vector<Rect> freeRectangles;

		/// Used to track that the packer produces proper packings.
		debug_run(DisjointRectCollection disjointRects;);

		/// Goes through the list of free rectangles and finds the best one to place a rectangle of given size into.
		/// Running time is Theta(|freeRectangles|).
		/// @param nodeIndex [out] The index of the free rectangle in the freeRectangles array into which the new
		///		rect was placed.
		/// @return A Rect structure that represents the placement of the new rect into the best free rectangle.
		Rect FindPositionForNewNode(int width, int height, FreeRectChoiceHeuristic rectChoice, int *nodeIndex) {
			Rect bestNode;
			memset(&bestNode, 0, sizeof(Rect));

			int bestScore = std::numeric_limits<int>::max();

			/// Try each free rectangle to find the best one for placement.
			for (size_t i = 0; i < freeRectangles.size(); ++i)
			{
				// If this is a perfect fit upright, choose it immediately.
				if (width == freeRectangles[i].width && height == freeRectangles[i].height)
				{
					bestNode.x = freeRectangles[i].x;
					bestNode.y = freeRectangles[i].y;
					bestNode.width = width;
					bestNode.height = height;
					bestScore = std::numeric_limits<int>::min();
					*nodeIndex = i;
					debug_assert(disjointRects.Disjoint(bestNode));
					break;
				}
				// If this is a perfect fit sideways, choose it.
				else if (height == freeRectangles[i].width && width == freeRectangles[i].height)
				{
					bestNode.x = freeRectangles[i].x;
					bestNode.y = freeRectangles[i].y;
					bestNode.width = height;
					bestNode.height = width;
					bestScore = std::numeric_limits<int>::min();
					*nodeIndex = i;
					debug_assert(disjointRects.Disjoint(bestNode));
					break;
				}
				// Does the rectangle fit upright?
				else if (width <= freeRectangles[i].width && height <= freeRectangles[i].height)
				{
					int score = ScoreByHeuristic(width, height, freeRectangles[i], rectChoice);

					if (score < bestScore)
					{
						bestNode.x = freeRectangles[i].x;
						bestNode.y = freeRectangles[i].y;
						bestNode.width = width;
						bestNode.height = height;
						bestScore = score;
						*nodeIndex = i;
						debug_assert(disjointRects.Disjoint(bestNode));
					}
				}
				// Does the rectangle fit sideways?
				else if (height <= freeRectangles[i].width && width <= freeRectangles[i].height)
				{
					int score = ScoreByHeuristic(height, width, freeRectangles[i], rectChoice);

					if (score < bestScore)
					{
						bestNode.x = freeRectangles[i].x;
						bestNode.y = freeRectangles[i].y;
						bestNode.width = height;
						bestNode.height = width;
						bestScore = score;
						*nodeIndex = i;
						debug_assert(disjointRects.Disjoint(bestNode));
					}
				}
			}
			return bestNode;
		}

		/// Returns the heuristic score value for placing a rectangle of size width*height into freeRect. Does not try to rotate.
		static int ScoreByHeuristic(int width, int height, const Rect &freeRect, FreeRectChoiceHeuristic rectChoice) {
			switch (rectChoice)
			{
			case RectBestAreaFit: return ScoreBestAreaFit(width, height, freeRect);
			case RectBestShortSideFit: return ScoreBestShortSideFit(width, height, freeRect);
			case RectBestLongSideFit: return ScoreBestLongSideFit(width, height, freeRect);
			case RectWorstAreaFit: return ScoreWorstAreaFit(width, height, freeRect);
			case RectWorstShortSideFit: return ScoreWorstShortSideFit(width, height, freeRect);
			case RectWorstLongSideFit: return ScoreWorstLongSideFit(width, height, freeRect);
			default: assert(false); return std::numeric_limits<int>::max();
			}
		}
		// The following functions compute (penalty) score values if a rect of the given size was placed into the 
		// given free rectangle. In these score values, smaller is better.
		static int ScoreBestAreaFit(int width, int height, const Rect &freeRect) {
			return freeRect.width * freeRect.height - width * height;
		}
		static int ScoreBestShortSideFit(int width, int height, const Rect &freeRect) {
			int leftoverHoriz = abs(freeRect.width - width);
			int leftoverVert = abs(freeRect.height - height);
			int leftover = min(leftoverHoriz, leftoverVert);
			return leftover;
		}
		static int ScoreBestLongSideFit(int width, int height, const Rect &freeRect) {
			int leftoverHoriz = abs(freeRect.width - width);
			int leftoverVert = abs(freeRect.height - height);
			int leftover = max(leftoverHoriz, leftoverVert);
			return leftover;
		}
		static int ScoreWorstAreaFit(int width, int height, const Rect &freeRect) {
			return -ScoreBestAreaFit(width, height, freeRect);
		}
		static int ScoreWorstShortSideFit(int width, int height, const Rect &freeRect) {
			return -ScoreBestShortSideFit(width, height, freeRect);
		}
		static int ScoreWorstLongSideFit(int width, int height, const Rect &freeRect) {
			return -ScoreBestLongSideFit(width, height, freeRect);
		}

		/// Splits the given L-shaped free rectangle into two new free rectangles after placedRect has been placed into it.
		/// Determines the split axis by using the given heuristic.
		void SplitFreeRectByHeuristic(const Rect &freeRect, const Rect &placedRect, GuillotineSplitHeuristic method) {
			// Compute the lengths of the leftover area.
			const int w = freeRect.width - placedRect.width;
			const int h = freeRect.height - placedRect.height;

			// Placing placedRect into freeRect results in an L-shaped free area, which must be split into
			// two disjoint rectangles. This can be achieved with by splitting the L-shape using a single line.
			// We have two choices: horizontal or vertical.	

			// Use the given heuristic to decide which choice to make.

			bool splitHorizontal;
			switch (method)
			{
			case SplitShorterLeftoverAxis:
				// Split along the shorter leftover axis.
				splitHorizontal = (w <= h);
				break;
			case SplitLongerLeftoverAxis:
				// Split along the longer leftover axis.
				splitHorizontal = (w > h);
				break;
			case SplitMinimizeArea:
				// Maximize the larger area == minimize the smaller area.
				// Tries to make the single bigger rectangle.
				splitHorizontal = (placedRect.width * h > w * placedRect.height);
				break;
			case SplitMaximizeArea:
				// Maximize the smaller area == minimize the larger area.
				// Tries to make the rectangles more even-sized.
				splitHorizontal = (placedRect.width * h <= w * placedRect.height);
				break;
			case SplitShorterAxis:
				// Split along the shorter total axis.
				splitHorizontal = (freeRect.width <= freeRect.height);
				break;
			case SplitLongerAxis:
				// Split along the longer total axis.
				splitHorizontal = (freeRect.width > freeRect.height);
				break;
			default:
				splitHorizontal = true;
				assert(false);
			}

			// Perform the actual split.
			SplitFreeRectAlongAxis(freeRect, placedRect, splitHorizontal);
		}

		/// Splits the given L-shaped free rectangle into two new free rectangles along the given fixed split axis.
		/// This function will add the two generated rectangles into the freeRectangles array. 
		/// The caller is expected to remove the original rectangle from the freeRectangles array after that.
		void SplitFreeRectAlongAxis(const Rect &freeRect, const Rect &placedRect, bool splitHorizontal) {
			// Form the two new rectangles.
			Rect bottom;
			bottom.x = freeRect.x;
			bottom.y = freeRect.y + placedRect.height;
			bottom.height = freeRect.height - placedRect.height;

			Rect right;
			right.x = freeRect.x + placedRect.width;
			right.y = freeRect.y;
			right.width = freeRect.width - placedRect.width;

			if (splitHorizontal)
			{
				bottom.width = freeRect.width;
				right.height = placedRect.height;
			}
			else // Split vertically
			{
				bottom.width = placedRect.width;
				right.height = freeRect.height;
			}

			// Add the new rectangles into the free rectangle pool if they weren't degenerate.
			if (bottom.width > 0 && bottom.height > 0)
				freeRectangles.push_back(bottom);
			if (right.width > 0 && right.height > 0)
				freeRectangles.push_back(right);

			debug_assert(disjointRects.Disjoint(bottom));
			debug_assert(disjointRects.Disjoint(right));
		}
	};

}
