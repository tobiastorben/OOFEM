#ifndef MODEL_H
#define MODEL_H

#include <Eigen\Core>
#include "Element.h"
#include <vector>
#include <iostream>
#include "utilities.h"
using namespace Eigen;
class Model {

private:
	MatrixXd nodes;//x,y,z for each node in each row
	Element* elementTypes;//Array of each element type
	VectorXi elementTable;//Each entry contains the element type of the corresponding element
	std::vector<std::vector<int>> topology;//Each row contains the nodes for the corresponding element
	MatrixXi supports;//Each row contain the node number, and the direction to be fixed [0 = x, 1 = y, 2 = z, 3 = rotX, 4 = rotY, 5 = rotZ]
	MatrixXd loads;//Each row contains load type: [0 = point load, 1 = uniform loads], node/element, direction [0-5] and intensity [N,Nm,N/m ....]
	int nElem;
	int nNodes;
	int nDof;
	int dimensionality;

	MatrixXd Kglob;
	MatrixXd constrainedKglob;
	VectorXd b;
	VectorXd constrainedB;
	VectorXd u;


public:
	Model(std::vector<std::vector<int>> someTopology, MatrixXd someNodes, Element* someEtypes, VectorXi someEtable, MatrixXd loads, MatrixXi someSupports, int dim);
	void calcLoadVector();
	void solveSystem();
	MatrixXd getElementCoords(int eNum);//Return matrix where each row contains the coordanate of a node
	void assembleGlobalSystem();
	void applyBCs();
	VectorXd getU();
	void printModel();
	VectorXd calcReactions();
	VectorXd calcError();
	void reduceSystem();
	void applyDimensionality();

	


};

#endif