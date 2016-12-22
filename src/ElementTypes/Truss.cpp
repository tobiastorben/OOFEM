#include "Truss.h"
#include "utilities.h"
#include <Eigen/Core>

using namespace Eigen;

Truss::Truss(double someEA, double someL) {
	EA = someEA;
	L = someL;
	rotDof = true;
}

//void Truss::calcStiffness() {
//	MatrixXd K(2, 2);
//	K << 1, -1, 1, -1;
//	K = K*(EA / L);
//}

MatrixXd Truss::calcNodalLoads(MatrixXd coords, int dof, double intensity) {
	double force = 0.5*L*intensity;
	VectorXd forceVec = force*calcUnitVec(coords);
	MatrixXd nodalLoads(2, 3);
	nodalLoads << forceVec, forceVec;	
	return nodalLoads;
}

void Truss::assembleToGlobal(std::vector<int> Dofs, MatrixXd coords, MatrixXd& Kglob) {
	MatrixXd T(2, 6);

	double cx = (coords[1, 0] - coords[0, 0]) / L;
	double cy = (coords[1, 1] - coords[0, 1]) / L;
	double cz = (coords[1, 2] - coords[0, 2]) / L;

	T << cx, cy, cz, 0, 0, 0, 0, 0, 0, cx, cy, cz;

	MatrixXd Kaxial(2, 2);
	Kaxial << 1, -1, -1, 1;
	MatrixXd K = (EA/L)*T.transpose()*Kaxial*T;
	int globDof;

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 3; l++)
					Kglob[Dofs[i] + j, Dofs[k] + l] = K[3 * i + j, 3 * k + l];
			}
		}
	}
}