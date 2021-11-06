//
// @author   liyan
// @contact  lyan_dut@outlook.com
//
#include "Tester.hpp"


void run_single_ins(const string& ins_bench, const string& ins_name) {
	Environment env(ins_bench, "H", ins_name);
	AdaptiveSelecter asa(env, cfg);
	asa.run();
	asa.draw_ins();
	asa.record_fp(env.fp_path());
	asa.record_fp(env.fp_path_with_time());
	asa.draw_fp(env.fp_html_path(), cfg.beta);
	asa.draw_fp(env.fp_html_path_with_time(), cfg.beta);
	asa.record_log();
}

void run_all_ins() { for_each(ins_list.begin(), ins_list.end(), [](auto& ins) { run_single_ins(ins.first, ins.second); }); }

int main(int argc, char** argv) {

	//qap::test();

	//metis::test();

	//test::record_gsrc_init_sol();

	//test::test_floorplan_packer();

	//run_single_ins("GSRC", "n10");

	run_all_ins();

	return 0;
}
