//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <string>
#include <fstream>
#include <random>

namespace visualizer {

	using namespace std;

	class Random {
	public:
		using Generator = mt19937;


		Random(int seed) : rgen(seed) {}
		Random() : rgen(generateSeed()) {}


		static int generateSeed() {
			return static_cast<int>(time(nullptr) + clock());
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
				<< "    <title>SmartFloorplan Visualization</title>" << endl
				<< "  </head>" << endl
				<< "  <body>" << endl // style='text-align:center;'
				<< "    <svg width='" << W << "' height='" << H << "' viewBox='-50 -50 " << W + 100 << " " << H + 100 << "'>" << endl;
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
		void rect(double x, double y, double w, double h, bool d = false, const string &label = "") {
			rc.next();
			rect(x, y, w, h, d, label, rc.fcolor, rc.bcolor);
		}

		void wire(double x1, double y1, double x2, double y2) {
			x1 *= wx; y1 *= hx; x2 *= wx; y2 *= hx;
			ofs << "      <line x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' style='stroke:#" << rc.bcolor << "; stroke-width:2'/>" << endl << endl;
		}

		void circle(double x, double y, double r = 2) {
			x *= wx; y *= hx;
			ofs << "      <circle cx='" << x << "' cy='" << y << "' r='" << r << "' style='fill-opacity:0; stroke:#000000; stroke-width:2'/>" << endl << endl;
		}

		double wx;
		double hx;
		ofstream ofs;
		RandColor rc;
	};
}