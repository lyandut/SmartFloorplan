/** @file Rect.h
	@author Jukka Jyl鋘ki

	This work is released to Public Domain, do whatever you want with it.
*/
#pragma once

#include <vector>
#include <string>

struct Rect {
	int id;
	int gid; // 分组id
	int x;
	int y;
	int width;
	int height;
	int area;
};

struct Block {
	std::string name;
	int width, height;
	int area;
	int x_coordinate, y_coordinate;
};

struct Terminal {
	std::string name;
	int x_coordinate, y_coordinate;
};

struct Net {
	int degree;
	std::vector<int> block_list;
	std::vector<int> terminal_list;
};

/// outline分组边界，划分时使用
struct Boundary {
	double x;
	double y;
	double width;
	double height;
};

/// Represents a single level (a horizontal line) of the skyline/horizon/envelope.
struct SkylineNode {
	/// The starting x-coordinate (leftmost).
	int x;
	/// The y-coordinate of the skyline level line.
	int y;
	/// The line width. The ending coordinate (inclusive) will be x+width-1.
	int width;
};

using Skyline = std::vector<SkylineNode>;

struct DisjointRects {
	std::vector<Rect> rects;

	void clear() { rects.clear(); }

	bool add(const Rect &r) {
		// Degenerate rectangles are ignored.
		if (r.width == 0 || r.height == 0) return true;

		if (!disjoint(r)) return false;

		rects.push_back(r);
		return true;
	}

	bool disjoint(const Rect &r) const {
		// Degenerate rectangles are ignored.
		if (r.width == 0 || r.height == 0) return true;

		for (const auto &rect : rects)
			if (!disjoint(rect, r))
				return false;

		return true;
	}

	static bool disjoint(const Rect &a, const Rect &b) {
		if (a.x + a.width <= b.x || b.x + b.width <= a.x ||
			a.y + a.height <= b.y || b.y + b.height <= a.y)
			return true;
		return false;
	}
};
