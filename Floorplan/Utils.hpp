//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <string>
#include <ctime>
#include <sstream>
#include <iomanip>

namespace utils {

	using namespace std;

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

	// 跳过文件中的n行 
	static void skip(FILE *file, int line_num) {
		char linebuf[1000];
		for (int i = 0; i < line_num; ++i) {
			fgets(linebuf, sizeof(linebuf), file); // skip
		}
	}

	// '01'法求组合问题 [todo] bug只能调用一次
	bool next_combination(const vector<int> &a, int n, int k, vector<int> &comb, vector<int> &ncomb) {
		static vector<bool> index(n, false);
		static bool first_comb = true;  // 第一个组合：选中前k个位置
		static auto has_done = [=]() {  // 终止条件：  最后k个位置全变成1
			for (int i = n - 1; i >= n - k; i--)
				if (!index[i])
					return false;
			return true;
		};

		if (first_comb) {
			for (int i = 0; i < k; ++i) { index[i] = true; }
			comb.clear(); comb.reserve(k);
			ncomb.clear(); ncomb.reserve(n - k);
			for (int i = 0; i < n; ++i) {
				if (index[i])
					comb.push_back(a[i]);
				else
					ncomb.push_back(a[i]);
			}
			first_comb = false;
			return true;
		}

		if (!has_done()) {
			for (int i = 0; i < n - 1; ++i) {
				// 找到第一个“10”组合将其变成"01"组合
				if (index[i] && !index[i + 1]) {
					index[i] = false;
					index[i + 1] = true;
					// 将"01"组合左边的1移到最左边
					int count = 0;
					for (int j = 0; j < i; ++j) {
						if (index[j]) {
							index[j] = false;
							index[count++] = true;
						}
					}
					// 记录当前组合
					comb.clear();
					ncomb.clear();
					for (int l = 0; l < n; ++l) {
						if (index[l])
							comb.push_back(a[l]);
						else
							ncomb.push_back(a[l]);
					}
						
					return true;
				}
			}
		}

		return false;
	}
}
