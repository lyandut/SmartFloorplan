//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#pragma once

/// 算法参数设置
struct Config {
	unsigned int random_seed; // 随机种子，release版本使用`random_device{}();`
	float init_fill_ratio;    // 初始填充率，建议0.5f
	int ub_time;              // ASA超时时间
	int dimension;            // QAP分组维度. e.g.dimension=5则划分25个分组
	int ub_iter;              // RLS最大迭代次数

	enum LevelCandidateWidth {
		CombRotate, // 考虑组合及旋转的所有情况
		CombShort,  // 考虑短边的组合
		Interval    // 以间隔等距划分宽度区间
	} level_asa_cw; // 建议Interval
	
	enum LevelGraphConnection {
		Direct,      // 仅考虑直连的边，非直连的边权重定义为 0
		Indirect     // 考虑非直连的边，非直连的边权重定义为 `1/最短路跳数`
	} level_qapc_gc; // 仅实现Direct

	enum LevelFlow {
		Kway,
		Recursive,
	} level_qapc_flow; // 建议Kway

	enum LevelDistance {
		EuclideanDis,   // 欧几里得距离
		ManhattanDis,   // 曼哈顿距离
		ChebyshevDis,   // 切比雪夫距离
		EuclideanSqrDis // 欧式平方距离
	} level_qapc_dis;

	/// 分组搜索策略
	enum LevelGroupSearch {
		None,		     // 不分组（纯优化面积）
		Selfishly,       // 仅当前分组
		NeighborAll,     // 当前分组和全部邻居分组
		NeighborPartial  // [todo] 当前分组和左/下邻居分组，或可考虑按百分比选取一部分矩形
	} level_fbp_gs;

	/// [deprecated] 放置策略
	enum LevelHeuristicSearch {
		MinHeightFit,    // 最小高度rectangle，不考虑靠skyline右侧放置
		MinWasteFit,     // 最小浪费rectangle，不考虑靠skyline右侧放置
		BottomLeftScore  // 最下最左skyline，考虑靠skyline右侧放置
	} level_fbp_hs;      // RLS仅支持BottomLeftScore
};