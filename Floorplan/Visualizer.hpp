//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "Instance.hpp"
#include "Config.hpp"

/// 基本数据结构定义
namespace utils_visualize {

	namespace bg = boost::geometry;

	using T = double;

	// n 维点坐标
	template<size_t dimension = 2>
	using point_base = bg::model::point<T, dimension, bg::cs::cartesian>;

	/*****************************
	* 以下定义全部针对二维点坐标 *
	******************************/

	// 二维坐标点
	using point_t = bg::model::d2::point_xy<T>;
	const point_t origin_point(0, 0);  // 坐标原点

	// 曲线
	using linestring_t = bg::model::linestring<point_t>;

	// 多边形（包括0个或者多个内环 inner rings，逆时针，起点=终点）
	using polygon_t = bg::model::polygon<point_t, false, true>;

	// 环（封闭曲线，逆时针，起点=终点）
	//using ring_t = bg::model::ring<point_t, false, true>;
	using ring_t = polygon_t::ring_type;

	// 点集合
	using multi_point_t = bg::model::multi_point<point_t>;

	// 曲线集合
	using multi_linestring_t = bg::model::multi_linestring<linestring_t>;

	// 多边形集合
	using multi_polygon_t = bg::model::multi_polygon<polygon_t>;

	// 矩形
	using box_t = bg::model::box<point_t>;

	// 线段（坐标点对）
	using segment_t = bg::model::segment<point_t>;

	// 多边形平移
	static void translatePolygon(const polygon_t &poly, polygon_t &translate_poly, T x, T y) {
		bg::strategy::transform::translate_transformer<double, 2, 2> translate_strategy(x, y);
		bg::transform(poly, translate_poly, translate_strategy);
	}

	// 多边形绕原点旋转一定角度
	static void rotatePolygon(const polygon_t &poly, polygon_t &rotate_poly, int angle) {
		bg::strategy::transform::rotate_transformer<bg::degree, double, 2, 2> rotate_strategy(angle);
		bg::transform(poly, rotate_poly, rotate_strategy);
	}

	// 求多边形的包络矩形
	using rectangle_t = std::pair<T, T>;
	static rectangle_t getEnvelope(const polygon_t &poly, box_t &envelope) {
		bg::envelope(poly, envelope);
		return {
			bg::get<bg::max_corner, 0>(envelope) - bg::get<bg::min_corner, 0>(envelope), // width
			bg::get<bg::max_corner, 1>(envelope) - bg::get<bg::min_corner, 1>(envelope)  // height
		};
	}
}

namespace utils_visualize_transform {

	using namespace utils_visualize;

	static vector<box_t> block_to_boxes(const Instance &ins) {
		vector<box_t> block_boxes;
		block_boxes.reserve(ins.get_block_num());
		for (auto &b : ins.get_blocks()) {
			block_boxes.emplace_back(point_t(b.x_coordinate, b.y_coordinate),
				point_t(b.x_coordinate + b.width, b.y_coordinate + b.height));
		}
		return block_boxes;
	}

	static vector<point_t> terminal_to_points(const Instance &ins) {
		vector<point_t> terminal_points;
		terminal_points.reserve(ins.get_terminal_num());
		for (auto &t : ins.get_terminals()) {
			terminal_points.emplace_back(t.x_coordinate, t.y_coordinate);
		}
		return terminal_points;
	}

	static vector<box_t> dst_to_boxes(const vector<rbp::Rect> &dst) {
		vector<box_t> dst_boxes;
		dst_boxes.reserve(dst.size());
		for (const auto &r : dst) {
			box_t new_box{ point_t(r.x, r.y), point_t(r.x + r.width, r.y + r.height) };
			for (const auto &old_box : dst_boxes) { assert(!bg::overlaps(new_box, old_box)); }
			dst_boxes.push_back(new_box);
		}
		return dst_boxes;
	}

	static vector<box_t> wire_to_boxes(const Instance &ins, vector<rbp::Rect> &rects, Config::LevelWireLength method) {
		vector<box_t> wire_boxes;
		wire_boxes.reserve(ins.get_net_num());
		sort(rects.begin(), rects.end(), [](auto &lhs, auto &rhs) { return lhs.id < rhs.id; }); // 确保输入顺序
		for (auto &net : ins.get_net_list()) {
			double max_x = 0, min_x = numeric_limits<double>::max();
			double max_y = 0, min_y = numeric_limits<double>::max();
			for (int b : net.block_list) {
				double pin_x = rects.at(b).x + rects.at(b).width / 2.0;
				double pin_y = rects.at(b).y + rects.at(b).height / 2.0;
				max_x = max(max_x, pin_x);
				min_x = min(min_x, pin_x);
				max_y = max(max_y, pin_y);
				min_y = min(min_y, pin_y);
			}
			if (method == Config::LevelWireLength::BlockAndTerminal) {
				for (int t : net.terminal_list) {
					double pad_x = ins.get_terminals().at(t).x_coordinate;
					double pad_y = ins.get_terminals().at(t).y_coordinate;
					max_x = max(max_x, pad_x);
					min_x = min(min_x, pad_x);
					max_y = max(max_y, pad_y);
					min_y = min(min_y, pad_y);
				}
			}
			wire_boxes.push_back({ {min_x,min_y},{max_x,max_y} });
		}
		return wire_boxes;
	}
}

namespace utils_visualize_drawer {

	class Random {
	public:
		using Generator = std::mt19937;


		Random(int seed) : rgen(seed) {}
		Random() : rgen(generateSeed()) {}


		static int generateSeed() {
			return static_cast<int>(std::time(nullptr) + std::clock());
		}

		Generator::result_type operator()() { return rgen(); }

		// pick with probability of (numerator / denominator).
		bool isPicked(unsigned numerator, unsigned denominator) {
			return ((rgen() % denominator) < numerator);
		}

		// pick from [min, max).
		int pick(int min, int max) {
			return ((rgen() % (max - min)) + min);
		}
		// pick from [0, max).
		int pick(int max) {
			return (rgen() % max);
		}

	protected:
		Generator rgen;
	};

	struct RandColor {
		static constexpr auto ColorCodeChar = "0123456789ABCDEF";
		static constexpr int ColorCodeBase = 16;
		static constexpr int ColorCodeLen = 6;

		void next() {
			for (int i = 0; i < ColorCodeLen; ++i) {
				int c = r.pick(ColorCodeBase);
				bcolor[i] = ColorCodeChar[c];
				fcolor[i] = ColorCodeChar[(c > (ColorCodeBase / 2)) ? 0 : (ColorCodeBase - 1)]; // (c + ColorCodeBase / 2) % ColorCodeBase
			}
		}

		char fcolor[ColorCodeLen + 1] = { 0 }; // front color.
		char bcolor[ColorCodeLen + 1] = { 0 }; // background color.
		Random r;
	};

	struct Drawer {
		static constexpr double W = 800;
		static constexpr double H = 800;

		Drawer(string path, double width, double height) : ofs(path), wx(W / width), hx(H / height) { begin(); }
		~Drawer() { end(); }

		void begin() {
			ofs << "<!DOCTYPE html>" << endl
				<< "<html>" << endl
				<< "  <head>" << endl
				<< "    <meta charset='utf-8'>" << endl
				<< "    <title>2D Packing/Cutting Visualization</title>" << endl
				<< "  </head>" << endl
				<< "  <body>" << endl // style='text-align:center;'
				<< "    <svg width='" << W << "' height='" << H << "' viewBox='0 0 " << W << " " << H << "'>" << endl;
		}
		void end() {
			ofs << "    </svg>" << endl
				<< "  </body>" << endl
				<< "</html>" << endl;
		}

		void rect(double x, double y, double w, double h, bool d, const string &label, const string &fcolor, const string &bcolor) {
			if (d) { swap(w, h); }
			x *= wx; y *= hx; w *= wx; h *= hx;
			ofs << "      <rect x='" << x << "' y='" << y << "' width='" << w << "' height='" << h << "' style='fill:#" << bcolor << "; stroke:black; stroke-width:2'/>" << endl
				<< "      <text x='" << x + w / 2 << "' y='" << y + h / 2 << "' text-anchor='middle' alignment-baseline='middle' style='fill:#" << fcolor << "'>" << label << "</text>" << endl << endl;
		}
		void rect(double x, double y, double w, double h, bool d, const string &label) {
			rc.next();
			rect(x, y, w, h, d, label, rc.fcolor, rc.bcolor);
		}

		void line(double x1, double y1, double x2, double y2, int layer) {
			static vector<int> cutWidth({ 8, 8, 6 });
			static vector<string> cutColor({ "red", "blue", "orange" });
			x1 *= wx; y1 *= hx; x2 *= wx; y2 *= hx;
			ofs << "      <line x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' stroke-dasharray='12, 4' stroke='" << cutColor[layer] << "' stroke-width='" << cutWidth[layer] << "'/>" << endl << endl;
		}

		void circle(double x, double y, double r) {
			x *= wx; y *= hx;
			ofs << "      <circle cx='" << x << "' cy='" << y << "' r='" << r << "' style='fill-opacity:0; stroke:#000000; stroke-width:2'/>" << std::endl;
		}

		double wx;
		double hx;
		ofstream ofs;
		RandColor rc;
	};
}