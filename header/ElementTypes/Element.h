#ifndef ELEMENT_H
#define ELEMENT_H

#include <Eigen/Core>
#include <vector>
using namespace Eigen;

class Element {

public:
	static const bool hasRotDof;
	//virtual VectorXd calcStress(VectorXd u) = 0;
	//virtual VectorXd calcStrain(VectorXd u) = 0;
	virtual void assembleToGlobal(std::vector<int> dofs, MatrixXd coords, MatrixXd& Kglob) = 0;
	virtual MatrixXd calcNodalLoads(MatrixXd coords, int dof, double intensity) = 0;//Return a matrix where each row contains the forcevector [0-5] for each node
};

#endif