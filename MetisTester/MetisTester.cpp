#include <metis.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

vector<idx_t> func(vector<idx_t> &xadj, vector<idx_t> &adjncy, vector<idx_t> &adjwgt) {
	idx_t nVertices = xadj.size() - 1;   // 节点数
	idx_t nEdges = adjncy.size() / 2;    // 边数
	idx_t nWeights = 1;                  // 节点权重维数
	idx_t nParts = 2;                    // 子图个数
	idx_t objval;
	std::vector<idx_t> part(nVertices, 0);

	/* "Note
 	    This function should be used to partition a graph into a large number of partitions(greater than 8).
		If a small number of partitions is desired, the METIS_PartGraphRecursive should be used instead, 
		as it produces somewhat better partitions."
	*/
	int ret = METIS_PartGraphKway(&nVertices, &nWeights, xadj.data(), adjncy.data(),
		NULL, NULL, adjwgt.data(), &nParts, NULL,
		NULL, NULL, &objval, part.data());
	std::cout << ret << std::endl;
	std::cout << objval << std::endl;
	for (unsigned part_i = 0; part_i < part.size(); part_i++) {
		std::cout << part_i << " " << part[part_i] << std::endl;
	}

	ret = METIS_PartGraphRecursive(&nVertices, &nWeights, xadj.data(), adjncy.data(),
		NULL, NULL, adjwgt.data(), &nParts, NULL,
		NULL, NULL, &objval, part.data());
	std::cout << ret << std::endl;
	std::cout << objval << std::endl;
	for (unsigned part_i = 0; part_i < part.size(); part_i++) {
		std::cout << part_i << " " << part[part_i] << std::endl;
	}
	return part;
}


int main() {
	ifstream ingraph("Instance/metisgraphs/graph.txt");
	if (!ingraph) {
		cout << "打开文件失败！" << endl;
		exit(1);//失败退回操作系统    
	}
	int vexnum, edgenum;
	string line;
	getline(ingraph, line);
	istringstream tmp(line);
	tmp >> vexnum >> edgenum;
	vector<idx_t> xadj(0);
	vector<idx_t> adjncy(0); //点的id从0开始
	vector<idx_t> adjwgt(0);

	idx_t a, w;
	for (int i = 0; i < vexnum; i++) {
		xadj.push_back(adjncy.size());
		getline(ingraph, line);
		istringstream tmp(line);
		while (tmp >> a >> w) {
			adjncy.push_back(a);
			adjwgt.push_back(w);
		}
	}
	xadj.push_back(adjncy.size());

	ingraph.close();

	vector<idx_t> part = func(xadj, adjncy, adjwgt);
	ofstream outpartition("partition.txt");
	if (!outpartition) {
		cout << "打开文件失败！" << endl;
		exit(1);
	}
	for (int i = 0; i < part.size(); i++) {
		outpartition << i << " " << part[i] << endl;
	}
	outpartition.close();
}
