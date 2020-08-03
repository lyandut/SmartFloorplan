#include "Instance.hpp"
#include "FloorplanBinPack.hpp"
#include <cmath>

void test_floorplan_bin_pack(Instance &ins) {

	float dead_ratio = 0.14f;
	int bin_width = sqrt(ins.get_total_area() * (1 + dead_ratio));
	int bin_height = bin_width;

	printf("Initializing bin to size %dx%d.\n", bin_width, bin_height);
	vector<rbp::Rect> src = ins.get_rects();
	vector<rbp::Rect> dst;
	dst.reserve(src.size());
	fbp::FloorplanBinPack fbp_solver(src, 3000, 16000);

	printf("Perform the packing...\n");
	fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinHeightFit);
	//fbp_solver.insert_greedy_fit(dst, fbp::FloorplanBinPack::LevelSelfishly, fbp::FloorplanBinPack::LevelMinWasteFit);
	//fbp_solver.insert_bottom_left_score(dst, fbp::FloorplanBinPack::LevelGroupSearch::LevelSelfishly);

	// Test success or failure.
	if (src.empty()) { printf("Successful! Occupancy Ratio: %.2f%%\n", fbp_solver.Occupancy()*100.f); }
	else { printf("Failed!\n"); }
	printf("Done. All rectangles packed.\n");
}

int main(int argc, char **argv) {
	//Environment env("GSRC", "H", "n10");
	Environment env("MCNC", "H", "ami49");
	Instance ins(env);
	ins.read_instance();

	test_floorplan_bin_pack(ins);

	system("pause");
	return 0;
}
