//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <algorithm>
#include <cassert>

#include "Data.hpp"
#include "Utils.hpp"

using namespace std;

class Environment {
public:
	Environment(const string& bench, const string& type, const string& name) :
		_ins_bench(bench), _ins_type(type), _ins_name(name) {}

	string blocks_path() const { return instance_dir() + benchmark_dir() + type_dir() + _ins_name + ".blocks"; }
	string nets_path() const { return instance_dir() + benchmark_dir() + type_dir() + _ins_name + ".nets"; }
	string pl_path() const { return instance_dir() + benchmark_dir() + type_dir() + _ins_name + ".pl"; }
	string pl_html_path() const { return instance_dir() + benchmark_dir() + type_dir() + _ins_name + ".html"; }
	string fp_path() const { return solution_dir() + benchmark_dir() + _ins_name + ".fp"; }
	string fp_path_with_time() const { return solution_dir() + benchmark_dir() + _ins_name + "." + utils::Date::to_long_str() + ".fp"; }
	string fp_html_path() const { return solution_dir() + benchmark_dir() + _ins_name + ".html"; }
	string fp_html_path_with_time() const { return solution_dir() + benchmark_dir() + _ins_name + "." + utils::Date::to_long_str() + ".html"; }
	string log_path() const { return solution_dir() + _ins_bench + ".csv"; }

private:
	static string instance_dir() { return "Instance/"; }
	static string solution_dir() { return "Solution/"; }
	string benchmark_dir() const { return _ins_bench == "MCNC" ? "MCNC/" : "GSRC/"; }
	string type_dir() const { return _ins_type == "H" ? "HARD/" : "SOFT/"; }

public:
	const string _ins_bench;
	const string _ins_type;
	const string _ins_name;
};

class Instance {
public:
	Instance(const Environment& env) : _env(env), _fixed_width(0), _fixed_height(0) { read_instance(); }

	/// 待打包矩形列表，默认w≤h避免后续计算增加无用判断
	vector<Rect> get_rects(bool rotated = true) const {
		vector<Rect> rects;
		rects.reserve(_block_num);
		if (rotated) {
			for (int i = 0; i < _block_num; ++i) {
				rects.push_back({ i,
					_blocks[i].x_coordinate,
					_blocks[i].y_coordinate,
					min(_blocks[i].width, _blocks[i].height),
					max(_blocks[i].width, _blocks[i].height) });
			}
		}
		else {
			for (int i = 0; i < _block_num; ++i) {
				rects.push_back({ i,
					_blocks[i].x_coordinate,
					_blocks[i].y_coordinate,
					_blocks[i].width,
					_blocks[i].height });
			}
		}
		return rects;
	}

	const vector<Block>& get_blocks() const { return _blocks; }
	const vector<Terminal>& get_terminals() const { return _terminals; }
	const vector<Net>& get_netlist() const { return _nets; }

	int get_block_num() const { return _block_num; }
	int get_terminal_num() const { return _terminal_num; }
	int get_net_num() const { return _net_num; }
	int get_total_area() const { return _total_area; }

	int get_fixed_width() const { return _fixed_width; }
	int get_fixed_height() const { return _fixed_height; }

private:
	void read_instance() {
		read_blocks();
		read_nets();
		read_pl();
	}

	/// read .blocks file.
	void read_blocks() {
		auto file = fopen(_env.blocks_path().c_str(), "r");
		if (file == 0) {
			fprintf(stderr, "%s: no blocks file\n", _env.blocks_path().c_str());
			return;
		}

		utils::skip(file, 6);
		fscanf(file, "NumHardRectilinearBlocks : %d\n", &_block_num);
		fscanf(file, "NumTerminals : %d\n", &_terminal_num);

		_blocks.resize(_block_num);
		_total_area = 0;
		for (int i = 0; i < _block_num; ++i) {
			int a, b, c, d, e, f; // ignore parameter
			char block_name[10] = { 0 };
			fscanf(file, "%s hardrectilinear 4 (%d, %d) (%d, %d) (%d, %d) (%d, %d)\n",
				block_name, &a, &b, &c, &d, &_blocks[i].width, &_blocks[i].height, &e, &f);
			_blocks[i].name = block_name;
			_blocks[i].area = _blocks[i].width * _blocks[i].height;
			_total_area += _blocks[i].area;
		}

		_terminals.resize(_terminal_num);
		for (int i = 0; i < _terminal_num; ++i) {
			char terminal_name[10] = { 0 };
			fscanf(file, "%s terminal\n", terminal_name);
			_terminals[i].name = terminal_name;
		}

		fprintf(stdout, "%s: blocks file read success.\n", _env.blocks_path().c_str());
		fclose(file);
		return;
	}

	/// read .nets file.
	void read_nets() {
		auto file = fopen(_env.nets_path().c_str(), "r");
		if (file == 0) {
			fprintf(stderr, "%s: no nets file\n", _env.nets_path().c_str());
			return;
		}

		utils::skip(file, 5);
		fscanf(file, "NumNets : %d\n", &_net_num);
		fscanf(file, "NumPins : %d\n", &_pin_num);

		_nets.resize(_net_num);
		for (int i = 0; i < _net_num; ++i) {
			fscanf(file, "NetDegree : %d\n", &_nets[i].degree);
			for (int d = 0; d < _nets[i].degree; ++d) {
				char tmp_char[10];
				fscanf(file, "%s %*[^\n]%*c", tmp_char);
				if (tmp_char[0] == '#') { --d; }
				else {
					string tmp_name = tmp_char;
					auto block_iter = find_if(_blocks.begin(), _blocks.end(),
						[&tmp_name](Block& block) { return  block.name == tmp_name; });
					if (block_iter != _blocks.end()) { // it's a block
						block_iter->net_ids.push_back(i);
						_nets[i].block_list.push_back(distance(_blocks.begin(), block_iter));
					}
					else { // it's a terminal
						auto terminal_iter = find_if(_terminals.begin(), _terminals.end(),
							[&tmp_name](Terminal& terminal) { return terminal.name == tmp_name; });
						assert(terminal_iter != _terminals.end());
						terminal_iter->net_ids.push_back(i);
						_nets[i].terminal_list.push_back(distance(_terminals.begin(), terminal_iter));
					}
				}
			}
		}

		fprintf(stdout, "%s: nets file read success.\n", _env.nets_path().c_str());
		fclose(file);
		return;
	}

	/// read .pl file.
	void read_pl() {
		auto file = fopen(_env.pl_path().c_str(), "r");
		if (file == 0) {
			fprintf(stderr, "%s: no pl file\n", _env.pl_path().c_str());
			return;
		}

		utils::skip(file, 5);

		for (int i = 0; i < _block_num; ++i) {
			char block_name[10];
			int x, y;
			fscanf(file, "%s %d %d\n", block_name, &x, &y);
			auto iter = find_if(_blocks.begin(), _blocks.end(),
				[&block_name](Block& block) { return block.name == string(block_name); });
			assert(iter != _blocks.end());
			iter->x_coordinate = x;
			iter->y_coordinate = y;
			_fixed_width = max(_fixed_width, x);
			_fixed_height = max(_fixed_height, y);
		}

		for (int i = 0; i < _terminal_num; ++i) {
			char terminal_name[10];
			int x, y;
			fscanf(file, "%s %d %d\n", terminal_name, &x, &y);
			auto iter = find_if(_terminals.begin(), _terminals.end(),
				[&terminal_name](Terminal& terminal) { return terminal.name == string(terminal_name); });
			assert(iter != _terminals.end());
			iter->x_coordinate = x;
			iter->y_coordinate = y;
			_fixed_width = max(_fixed_width, x);
			_fixed_height = max(_fixed_height, y);
		}

		fprintf(stdout, "%s: pl file read success.\n", _env.pl_path().c_str());
		fclose(file);
		return;
	}

private:
	const Environment& _env;

	// fixed outline info
	int _fixed_width;
	int _fixed_height;

	// blocks info
	vector<Block> _blocks;
	int _block_num;
	int _total_area;

	// terminals info
	vector<Terminal> _terminals;
	int _terminal_num;

	// nets info
	vector<Net> _nets;
	int _net_num;
	int _pin_num;
};
