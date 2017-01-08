#include <Eigen/Core>
#include <iostream>
#include "Element.h"
#include "Truss.h"
#include "Frame.h"
#include "Model.h"
#include <math.h>

#define PI 3.14159265

using namespace Eigen;

int main(int argc, char** argv) {

	int dimensionality = 2;

	std::vector<std::vector<int>> topology;
	
	topology = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 }, { 5, 6 }, {6,7},{ 7, 8 }, { 8, 9 }, { 9, 10 } };

	MatrixXd nodes(11, 3);
	double h = sqrt(0.5);
	double w = 1000;
	nodes << 0, 0, 0, w, 0, 0, 2*w, 0, 0, 3*w, 0, 0, 4*w, 0, 0, 5*w, 0, 0, 6*w, 0, 0, 7*w, 0, 0, 8*w, 0, 0, 9*w, 0, 0, 10*w, 0, 0;


	Truss bar(210e3, 0.1);
	Frame beam(210e9, 0.1, 0.1, 0.1, 80e9, 5.0/6);

	Element* eTypes = &beam;

	VectorXi eTable = VectorXi::Zero(10);

	MatrixXd loads(1, 4);
	loads << 0, 5, 1, -5;

	MatrixXi supports(6, 2);
	supports << 0, 0, 0, 1,0,2,0,3,0,4,10, 1;

	Model model(topology, nodes, eTypes, eTable, loads, supports, dimensionality);

	model.calcLoadVector();

	model.assembleGlobalSystem();

//	model.applyDimensionality();

	model.applyBCs();

	//model.printModel();

	model.solveSystem();

	std::cout << "\n\nNodal displacements:\n\n" << model.getU();

	auto R = model.calcReactions();

	std::cout << "\n\nSupport reactions: \n\n" << R;

	auto error = model.calcError();

	std::cout << "\n\nError: \n\n" << error;

	std::cin.ignore();

	return 0;
}