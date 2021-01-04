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

static constexpr int INF = 0x3f3f3f3f;

static std::vector<std::pair<std::string, std::string>> ins_list{
	{"MCNC", "apte"},{"MCNC", "xerox"},{"MCNC", "hp"},{"MCNC", "ami33"},{"MCNC", "ami49"},
	{"GSRC", "n10"},{"GSRC", "n30"},{"GSRC", "n50"},{"GSRC", "n100"},{"GSRC", "n200"},{"GSRC", "n300"}
};

/// 算法参数设置
struct Config {
	unsigned int random_seed = std::random_device{}(); // 随机种子，release版本使用`random_device{}();`
	//unsigned int random_seed = 2231403657;

	double alpha = 0.8, beta = 0.2;        // 控制面积、线长优化权重
	double lb_scale = 0.8, ub_scale = 1.2; // 控制分支数目、长宽比

	int ub_time = 3600; // ASA超时时间
	int ub_iter = 8192; // RLS最大迭代次数/BS最大树宽度

	int dimension = 5;  // QAP最大分组维度. e.g.dimension=5则划分25个分组，不足25个block则求上界

	enum LevelCandidateWidth {
		CombRotate,            // [deprecated] 考虑组合及旋转的所有情况
		CombShort,             // [deprecated] 考虑短边的组合
		Interval,              // 以间隔等距划分宽度区间
		Sqrt                   // 利用平方根控制长宽比
	} level_asa_cw = Interval;

	enum LevelFloorplanPacker {
		RandomLocalSearch,
		BeamSearch
	} level_asa_fbp = BeamSearch;

	enum LevelWireLength {
		BlockOnly,         // 计算block内部互连线长
		BlockAndTerminal   // [deprecated] 计算block和terminal互连线长
	} level_fbp_wl = BlockOnly;

	enum LevelObjDist {   // EDAthon2020P4目标函数中dist计算方法
		SqrHpwlDist,      // 整个网表半周长的平方
		SqrEuclideanDist, // 矩形对之间欧式平方距离
		SqrManhattanDist  // 矩形对之间曼哈顿平方距离
	} level_fbp_dist = SqrManhattanDist;

	enum LevelQAPCluster {
		On,
		Off // 不使用QAP
	} level_rls_qapc = Off;

	enum LevelGroupSearch {
		NoGroup,		 // 不分组，LevelQAPCluster::Off <==> LevelGroupSearch::NoGroup
		NeighborNone,    // 当前分组
		NeighborAll,     // 当前分组和全部邻居分组
		NeighborPartial  // 当前分组和左下邻居分组，[todo]考虑按百分比选取一部分矩形?
	} level_rls_gs = NoGroup;

	enum LevelGraphConnect {
		Direct,  // 仅考虑直连的边，非直连的边权重定义为 0
		Indirect // 考虑非直连的边，非直连的边权重定义为 `1/最短路跳数`
	} level_qapc_gc = Direct;

	enum LevelFlow {
		Kway,
		Recursive
	} level_qapc_flow = Kway;

	enum LevelDist {
		EuclideanDist,   // 欧几里得距离
		ManhattanDist,   // 曼哈顿距离
		ChebyshevDist    // 切比雪夫距离
	} level_qapc_dist = ManhattanDist;
} cfg;

std::ostream& operator<<(std::ostream &os, const Config &cfg) {
	switch (cfg.level_asa_fbp) {
	case Config::LevelFloorplanPacker::RandomLocalSearch: os << "RandomLocalSearch,"; break;
	case Config::LevelFloorplanPacker::BeamSearch: os << "BeamSearch,"; break;
	default: break;
	}

	switch (cfg.level_fbp_wl) {
	case Config::LevelWireLength::BlockOnly: os << "BlockOnly,"; break;
	case Config::LevelWireLength::BlockAndTerminal: os << "BlockAndTerminal,"; break;
	default: break;
	}

	switch (cfg.level_fbp_dist) {
	case Config::LevelObjDist::SqrHpwlDist: os << "SqrHpwlDist,"; break;
	case Config::LevelObjDist::SqrEuclideanDist: os << "SqrEuclideanDist,"; break;
	case Config::LevelObjDist::SqrManhattanDist: os << "SqrManhattanDist,"; break;
	default: break;
	}

	switch (cfg.level_rls_qapc) {
	case Config::LevelQAPCluster::On: os << "On,"; break;
	case Config::LevelQAPCluster::Off: os << "Off,"; break;
	default: break;
	}

	switch (cfg.level_rls_gs) {
	case Config::LevelGroupSearch::NoGroup: os << "NoGroup,"; break;
	case Config::LevelGroupSearch::NeighborNone: os << "NeighborNone,"; break;
	case Config::LevelGroupSearch::NeighborAll: os << "NeighborAll,"; break;
	case Config::LevelGroupSearch::NeighborPartial: os << "NeighborPartial,"; break;
	default: break;
	}

	switch (cfg.level_qapc_gc) {
	case Config::LevelGraphConnect::Direct: os << "Direct,"; break;
	case Config::LevelGraphConnect::Indirect: os << "Indirect,"; break;
	default: break;
	}

	switch (cfg.level_qapc_flow) {
	case Config::LevelFlow::Kway: os << "Kway,"; break;
	case Config::LevelFlow::Recursive:os << "Recursive,"; break;
	default: break;
	}

	switch (cfg.level_qapc_dist) {
	case Config::LevelDist::EuclideanDist: os << "EuclideanDist"; break;
	case Config::LevelDist::ChebyshevDist: os << "ChebyshevDist"; break;
	case Config::LevelDist::ManhattanDist: os << "ManhattanDist"; break;
	default: break;
	}

	return os;
}
