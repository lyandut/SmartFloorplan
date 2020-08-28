//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>

const int INF = 0x3f3f3f3f;

namespace utils {

	using namespace std;

	// 跳过文件中的n行 
	static void skip(FILE *file, int line_num) {
		char linebuf[1000];
		for (int i = 0; i < line_num; ++i) {
			fgets(linebuf, sizeof(linebuf), file); // skip
		}
	}

	class Date {
	public:
		// 返回表示日期格式的字符串，年-月-日
		static const string to_short_str() {
			ostringstream os;
			time_t now = time(0);
			os << put_time(localtime(&now), "%y-%m-%e");
			return os.str();
		}
		// 返回表示日期格式的字符串，连续数字
		static const string to_detail_str() {
			ostringstream os;
			time_t now = time(0);
			//os << put_time(localtime(&now), "%y-%m-%e-%H_%M_%S");
			os << put_time(localtime(&now), "%Y%m%d%H%M%S");
			return os.str();
		}
	};

	class Combination {
	public:
		Combination(const vector<int> &a, int k) : _a(a), _n(a.size()), _k(k), _index(a.size(), false), _first_comb(true) {}

		bool next_combination(vector<int> &comb, vector<int> &ncomb) {
			if (_first_comb) { // 第一个组合：选中前k个位置
				_first_comb = false;
				for (int i = 0; i < _k; ++i) { _index[i] = true; }
				record_combination(comb, ncomb);
				return true;
			}

			if (!has_done()) {
				for (int i = 0; i < _n - 1; ++i) {
					// 找到第一个“10”组合将其变成"01"组合
					if (_index[i] && !_index[i + 1]) {
						_index[i] = false;
						_index[i + 1] = true;
						// 将"01"组合左边的1移到最左边
						int count = 0;
						for (int j = 0; j < i; ++j) {
							if (_index[j]) {
								_index[j] = false;
								_index[count++] = true;
							}
						}
						record_combination(comb, ncomb);
						return true;
					}
				}
			}

			return false;
		}

	private:
		bool has_done() const { // 终止条件：最后k个位置全变成1
			for (int i = _n - 1; i >= _n - _k; --i)
				if (!_index[i])
					return false;
			return true;
		}

		void record_combination(vector<int> &comb, vector<int> &ncomb) const {
			comb.clear(); comb.reserve(_k);
			ncomb.clear(); ncomb.reserve(_n - _k);
			for (int i = 0; i < _n; ++i) {
				if (_index[i])
					comb.push_back(_a[i]);
				else
					ncomb.push_back(_a[i]);
			}
		}

	private:
		const vector<int> &_a;
		const int _n;
		const int _k;
		vector<bool> _index;
		bool _first_comb;
	};
}
