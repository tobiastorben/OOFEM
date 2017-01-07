#ifndef TRUSS_H
#define TRUSS_H

#include "Element.h"
#include <Eigen/Core>
#include <vector>

class Truss : public Element {
	private :
		double E;
		double A;
		
	public:
		Truss(double someE, double someA);
		VectorXd calcStress(VectorXd u);
		VectorXd calcStrain(VectorXd u);
		void assembleToGlobal(std::vector<int> eNodes, MatrixXd coords, MatrixXd& Kglob);
		MatrixXd calcNodalLoads(MatrixXd coords, int dof, double intensity);
		
};

#endif