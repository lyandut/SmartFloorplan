/** @file Rect.h
	@author Jukka Jyl鋘ki

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include <vector>
#include <cassert>
#include <cstdlib>
#include <utility>

#ifdef _DEBUG
/// debug_assert is an assert that also requires debug mode to be defined.
#define debug_assert(x) assert(x)
#define debug_run(x) x
#else
#define debug_assert(x)
#define debug_run(x)
#endif

namespace rbp {

	struct Rect
	{
		int id;
		int gid;

		int x;
		int y;
		int width;
		int height;
	};

	struct Boundary
	{
		double x;
		double y;
		double width;
		double height;
	};

	/// Returns true if a is contained in b.
	bool IsContainedIn(const Rect &a, const Rect &b)
	{
		return a.x >= b.x && a.y >= b.y
			&& a.x + a.width <= b.x + b.width
			&& a.y + a.height <= b.y + b.height;
	}

	/// 不重叠矩形集合
	class DisjointRectCollection
	{
	public:
		std::vector<Rect> rects;

		/// 添加一个不重叠矩形到集合
		bool Add(const Rect &r)
		{
			// Degenerate rectangles are ignored.
			if (r.width == 0 || r.height == 0)
				return true;

			if (!Disjoint(r))
				return false;
			rects.push_back(r);
			return true;
		}

		void Clear()
		{
			rects.clear();
		}

		/// 判断矩形r是否与已放置矩形重叠：重叠-false，不重叠-true
		bool Disjoint(const Rect &r) const
		{
			// Degenerate rectangles are ignored.
			if (r.width == 0 || r.height == 0)
				return true;

			for (size_t i = 0; i < rects.size(); ++i)
				if (!Disjoint(rects[i], r))
					return false;
			return true;
		}

		/// 判断两个矩形重叠：重叠-false，不重叠-true
		static bool Disjoint(const Rect &a, const Rect &b)
		{
			if (a.x + a.width <= b.x ||
				b.x + b.width <= a.x ||
				a.y + a.height <= b.y ||
				b.y + b.height <= a.y)
				return true;
			return false;
		}
	};

}
