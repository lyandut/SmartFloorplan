//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <list>
#include <numeric>

#include "Config.hpp"
#include "Instance.hpp"

namespace fbp {

	using namespace std;

	class FloorplanPacker {
	public:
		FloorplanPacker() = delete;

		FloorplanPacker(const Instance &ins, const vector<Rect> &src, int bin_width, default_random_engine &gen) :
			_ins(ins), _src(src), _bin_width(bin_width), _gen(gen), _objective(numeric_limits<double>::max()) {}

		const vector<Rect>& get_dst() const { return _dst; }

		double get_objective() const { return _objective; }

		int get_area() const { return _obj_area; }

		double get_fill_ratio() const { return _obj_fillratio; }

		double get_wirelength() const { return _obj_wirelength; }

		virtual void run(int, double, double, Config::LevelWireLength, Config::LevelGroupSearch) = 0;

	protected:
		// 输入
		const Instance &_ins;
		const vector<Rect> &_src;
		const int _bin_width;
		default_random_engine &_gen;

		// 优化目标
		vector<Rect> _dst;
		double _objective;
		int _obj_area;
		double _obj_fillratio;
		double _obj_wirelength;
	};

}
