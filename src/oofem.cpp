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

	Truss bar(210e9*0.01*0.01, 1);

	Element* eTypes = &bar;

	VectorXi eTable(3);
	eTable << 0, 0, 0;

	MatrixXd loads(1, 4);
	loads << 0, 3, 0, 1000;

	MatrixXi supports(3, 2);
	supports << 0, 0, 0, 1, 0, 2;

	Model model(topology, nodes, eTypes, eTable, loads, supports);

	model.assembleGlobalSystem();

	model.calcLoadVector();

	model.applyBCs();

	model.printModel();

	model.solveSystem();

	std::cout << model.getU();

	std::cin.ignore();

	return 0;
}