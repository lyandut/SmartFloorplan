//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#ifdef NDEBUG
#define debug_run(x)
#else
#define debug_run(x) x
#endif // !NDEBUG

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <random>
#include <unordered_map>

static constexpr int INF = 0x3f3f3f3f;

static std::vector<std::pair<std::string, std::string>> ins_list{
	{"MCNC", "apte"},{"MCNC", "xerox"},{"MCNC", "hp"},{"MCNC", "ami33"},{"MCNC", "ami49"},
	{"GSRC", "n10"},{"GSRC", "n30"},{"GSRC", "n50"},{"GSRC", "n100"},{"GSRC", "n200"},{"GSRC", "n300"}
};

static std::unordered_map<std::string, std::pair<int, int>> obj_map{
	{"apte", {47050000, 246000}}, {"xerox", {20340000, 379600}}, {"hp", {9460000, 149800}},
	{"ami33", {1220000, 59500}}, {"ami49", {37820000, 667000}},
	{"n10", {225242, 25788}}, {"n30", {216051, 79740}}, {"n50", {204295, 124326}},
	{"n100", {187880, 206269}}, {"n200", {183931, 389272}}, {"n300", {288702, 587739}}
};

/// 算法参数设置
struct Config {
	unsigned int random_seed = std::random_device{}(); // 随机种子，release版本使用`random_device{}()`

	double alpha = 0.5, beta = 0.5;        // 控制面积、线长权重
	double lb_scale = 0.8, ub_scale = 1.2; // 控制候选宽度数目、长宽比

	int ub_time = 3600; // ASA超时时间
	int ub_iter = 8192; // RLS最大迭代次数 or BS最大树宽度

	enum class LevelCandidateWidth {
		CombRotate, // [deprecated] 考虑组合及旋转的所有情况
		CombShort,  // [deprecated] 考虑短边的组合
		Interval,   // 以间隔等距划分宽度区间
		Sqrt        // 利用平方根控制长宽比
	} level_asa_cw = LevelCandidateWidth::Interval;

	enum class LevelFloorplanPacker {
		RandomLocalSearch,
		BeamSearch,
	} level_asa_fbp = LevelFloorplanPacker::BeamSearch;

	enum class LevelWireLength {
		Block,           // 计算block内部互连线长
		BlockAndTerminal // 计算block和terminal互连线长
	} level_fbp_wl = LevelWireLength::Block;

	// EDAthon2020P4目标函数：(Af+D)/Ac，https://blog.csdn.net/nobleman__/article/details/107880261
	enum class LevelObjDist {
		WireLengthDist,   // 线长
		SqrEuclideanDist, // 矩形对之间欧式平方距离的和
		SqrManhattanDist  // 矩形对之间曼哈顿平方距离的和
	} level_fbp_dist = LevelObjDist::SqrManhattanDist;
} cfg;

std::ostream& operator<<(std::ostream& os, const Config& cfg) {
	switch (cfg.level_asa_fbp) {
	case Config::LevelFloorplanPacker::RandomLocalSearch: os << "RandomLocalSearch,"; break;
	case Config::LevelFloorplanPacker::BeamSearch: os << "BeamSearch,"; break;
	default: break;
	}

	switch (cfg.level_fbp_wl) {
	case Config::LevelWireLength::Block: os << "Block,"; break;
	case Config::LevelWireLength::BlockAndTerminal: os << "BlockAndTerminal,"; break;
	default: break;
	}

	switch (cfg.level_fbp_dist) {
	case Config::LevelObjDist::WireLengthDist: os << "WireLengthDist"; break;
	case Config::LevelObjDist::SqrEuclideanDist: os << "SqrEuclideanDist"; break;
	case Config::LevelObjDist::SqrManhattanDist: os << "SqrManhattanDist"; break;
	default: break;
	}

	return os;
}
