#include <Eigen/Core>
#include <iostream>
#include "Element.h"
#include "Truss.h"
#include "Model.h"
#include <math.h>

#define PI 3.14159265

using namespace Eigen;

int main(int argc, char** argv) {

	int dimensionality = 2;

	std::vector<std::vector<int>> topology;
	
	topology = {{0, 1} ,{ 1, 2 }, { 2, 3 }, { 4, 5 }, { 5, 6 }, { 0, 4 }, { 4, 1 }, { 1, 5 }, { 5, 2 }, { 2, 6 }, { 6, 3 }};

	MatrixXd nodes(7, 3);
	double h = sqrt(0.5);
	double w = 1;
	nodes << 0, 0, 0, w, 0, 0, 2 * w, 0, 0, 3 * w, 0, 0, 0.5*w, -h, 0, 1.5*w, -h, 0, 2.5*w, -h, 0;


	Truss bar(210e3, 0.1);

	Element* eTypes = &bar;

	VectorXi eTable = VectorXi::Zero(11);

	MatrixXd loads(2, 4);
	loads << 0, 1, 1, -500, 0,2,1,-500;

	MatrixXi supports(3, 2);
	supports << 0, 0, 0, 1, 3, 1;

	Model model(topology, nodes, eTypes, eTable, loads, supports, dimensionality);

	model.calcLoadVector();

	model.assembleGlobalSystem();

	model.applyDimensionality();

	model.applyBCs();

	model.printModel();

	model.solveSystem();

	std::cout << "\n\nNodal displacements:\n\n" << model.getU();

	auto R = model.calcReactions();

	std::cout << "\n\nSupport reactions: \n\n" << R;

	auto error = model.calcError();

	std::cout << "\n\nError: \n\n" << error;

	std::cin.ignore();

	return 0;
}