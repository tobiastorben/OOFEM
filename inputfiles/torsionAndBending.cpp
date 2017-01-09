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
	
	topology = { { 0, 1 }, { 1, 2 } };

	MatrixXd nodes(3, 3);
	double L = 1;
	nodes << 0,0,0,
			 L, 0, 0,
			 L, 0, L;
			
	Frame beam(210e9, 0.1*0.03, 2.5e-6, 2.5e-6, 80e9, 5.0/6);

	Element* eTypes = &beam;

	VectorXi eTable = VectorXi::Zero(2);

	MatrixXd loads(1, 4);
	loads << 0, 2, 1, -5e3;

	MatrixXi supports(6, 2);
	supports << 0, 0,
		0, 1,
		0, 2,
		0, 3,
		0, 4,
		0, 5;
				


	Model model(topology, nodes, eTypes, eTable, loads, supports, dimensionality);

	model.calcLoadVector();

	model.assembleGlobalSystem();

	model.applyBCs();

	//model.printModel();

	model.solveSystem();

	printVecInCols("Nodal displacements:", model.getU());

	auto R = model.calcReactions();

	printVecInCols("Support Reactions:", R);

	auto error = model.calcError();

	printVecInCols("Error in nodal displacements", error);

	std::cin.ignore();

	return 0;
}