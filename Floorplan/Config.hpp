//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

#include <iostream>

/// 算法参数设置
struct Config {
	unsigned int random_seed; // 随机种子，release版本使用`random_device{}();`
	double init_fill_ratio;   // 初始填充率，建议0.5
	double alpha;             // 面积权重
	double beta;			  // 线长权重
	int dimension;            // QAP最大分组维度. e.g.dimension=5则划分25个分组，不足25个block则求上界
	int ub_time;              // ASA超时时间
	int ub_iter;              // RLS最大迭代次数

	enum LevelCandidateWidth {
		CombRotate, // 考虑组合及旋转的所有情况
		CombShort,  // 考虑短边的组合
		Interval    // 以间隔等距划分宽度区间
	} level_asa_cw; // 建议Interval

	enum LevelGraphConnection {
		Direct,      // 仅考虑直连的边，非直连的边权重定义为 0
		Indirect     // 考虑非直连的边，非直连的边权重定义为 `1/最短路跳数`
	} level_qapc_gc; // [todo] 当前仅实现Direct

	enum LevelFlow {
		Kway,
		Recursive,
	} level_qapc_flow; // 建议Kway

	enum LevelDistance {
		EuclideanDis,   // 欧几里得距离
		ManhattanDis,   // 曼哈顿距离
		ChebyshevDis,   // 切比雪夫距离
		EuclideanSqrDis // 欧式平方距离
	} level_qapc_dis;   // 默认ManhattanDis

	enum LevelGroupSearch {
		NoGroup,		 // 不分组
		Selfishly,       // 仅当前分组
		NeighborAll,     // 当前分组和全部邻居分组
		NeighborPartial  // [todo] 当前分组和左/下邻居分组，或可考虑按百分比选取一部分矩形
	} level_fbp_gs;      // 默认NeighborAll

	enum LevelWireLength {
		BlockOnly,       // 仅计算block互连线长，用于测试[Chen.2017]论文结果
		BlockAndTerminal // 计算block和terminal互连线长
	} level_fbp_wl;      // 默认BlockAndTerminal

	enum LevelObjNorm {
		NoNorm,         // 不归一化，仅在纯优化面积时配合`LevelGroupSearch::NoGroup`使用
		Average,        // 使用平均值归一化
		Minimum			// 使用最小值归一化
	} level_fbp_norm;

	/// [deprecated] 仅对贪心算法生效
	enum LevelHeuristicSearch {
		MinHeightFit,    // 最小高度rectangle，不考虑靠skyline右侧放置
		MinWasteFit,     // 最小浪费rectangle，不考虑靠skyline右侧放置
		BottomLeftScore  // 最下最左skyline，考虑靠skyline右侧放置
	} level_fbp_hs;      // RLS仅支持BottomLeftScore
};

std::ostream& operator<<(std::ostream &os, const Config &cfg) {
	switch (cfg.level_asa_cw) {
	case Config::LevelCandidateWidth::CombRotate: os << "CombRotate,"; break;
	case Config::LevelCandidateWidth::CombShort: os << "CombShort,"; break;
	case Config::LevelCandidateWidth::Interval: os << "Interval,"; break;
	default: break;
	}

	switch (cfg.level_qapc_gc) {
	case Config::LevelGraphConnection::Direct: os << "Direct,"; break;
	case Config::LevelGraphConnection::Indirect: os << "Indirect,"; break;
	default: break;
	}

	switch (cfg.level_qapc_flow) {
	case Config::LevelFlow::Kway: os << "Kway,"; break;
	case Config::LevelFlow::Recursive:os << "Recursive,"; break;
	default: break;
	}

	switch (cfg.level_qapc_dis) {
	case Config::LevelDistance::EuclideanDis: os << "EuclideanDis,"; break;
	case Config::LevelDistance::ChebyshevDis: os << "ChebyshevDis,"; break;
	case Config::LevelDistance::ManhattanDis: os << "ManhattanDis,"; break;
	case Config::LevelDistance::EuclideanSqrDis: os << "EuclideanSqrDis,"; break;
	default: break;
	}

	switch (cfg.level_fbp_gs) {
	case Config::LevelGroupSearch::NoGroup: os << "NoGroup,"; break;
	case Config::LevelGroupSearch::Selfishly: os << "Selfishly,"; break;
	case Config::LevelGroupSearch::NeighborAll: os << "NeighborAll,"; break;
	case Config::LevelGroupSearch::NeighborPartial: os << "NeighborPartial,"; break;
	default: break;
	}

	switch (cfg.level_fbp_wl) {
	case Config::LevelWireLength::BlockOnly: os << "BlockOnly,"; break;
	case Config::LevelWireLength::BlockAndTerminal: os << "BlockAndTerminal,"; break;
	default: break;
	}

	switch (cfg.level_fbp_norm) {
	case Config::LevelObjNorm::NoNorm: os << "NoNorm"; break;
	case Config::LevelObjNorm::Average: os << "Average"; break;
	case Config::LevelObjNorm::Minimum: os << "Minimum"; break;
	default: break;
	}

	return os;
}