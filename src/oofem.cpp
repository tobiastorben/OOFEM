#include <Eigen/Core>
#include <iostream>
#include "Element.h"
#include "Truss.h"
#include "Model.h"


using namespace Eigen;

int main(int argc, char** argv) {

	std::vector<std::vector<int>> topology;
	topology = { { 0, 1 }, { 1, 2 }, { 2, 3 } };

	MatrixXd nodes(4, 3);
	nodes << 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0;

	Truss bar(210e9, 1);

	Element* eTypes;
	eTypes[0] = bar;

	VectorXi eTable;
	eTable << 0, 0, 0;

	MatrixXd loads(5, 1);
	loads << 0, 4, 0, 1000;

	MatrixXi supports(6, 2);
	supports << 1, 0, 1, 1, 1, 2, 1, 3, 1, 4, 1, 5;

	Model model(topology, nodes, eTypes, eTable, loads, supports);



}