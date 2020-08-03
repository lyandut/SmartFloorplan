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
}
