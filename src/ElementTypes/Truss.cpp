#include "Truss.h"
#include "utilities.h"
#include <Eigen/Core>

using namespace Eigen;

const bool Truss::hasRotDof = false;

Truss::Truss(double someE, double someA) {
	E = someE;
	A = someA;
}


MatrixXd Truss::calcNodalLoads(MatrixXd coords, int dof, double intensity) {
	double L = calcL(coords);
	double force = 0.5*L*intensity;
	VectorXd forceVec = force*calcUnitVec(coords);
	MatrixXd nodalLoads(2, 3);
	nodalLoads << forceVec, forceVec;	
	return nodalLoads;
}

void Truss::assembleToGlobal(std::vector<int> eNodes, MatrixXd coords, MatrixXd& Kglob) {
	MatrixXd T(2, 6);
	double L = calcL(coords);
	double cx = (coords(1, 0) - coords(0, 0)) / L;
	double cy = (coords(1, 1) - coords(0, 1)) / L;
	double cz = (coords(1, 2) - coords(0, 2)) / L;

	T << cx, cy, cz, 0, 0, 0, 0, 0, 0, cx, cy, cz;

	MatrixXd Kaxial(2, 2);
	Kaxial << 1, -1, -1, 1;
	MatrixXd K = (E*A/L)*T.transpose()*Kaxial*T;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 3; l++)
					Kglob(6 * eNodes[i] + j, 6 * eNodes[k] + l) = Kglob(6 * eNodes[i] + j, 6 * eNodes[k] + l) + K(3 * i + j, 3 * k + l);
			}
		}
	}
}


