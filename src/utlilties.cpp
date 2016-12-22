#include "utilities.h"


Vector3d calcUnitVec(MatrixXd coords) {
	Vector3d vec = (coords.row(1) - coords.row(0));
	return vec.normalized();
}

VectorXi stdToEigenVec(std::vector<int> stdVec) {
	int n = stdVec.size();
	VectorXi eigenVec(n);
	for (int i = 0; i < n; i++) {
		eigenVec << stdVec[i];
	}
	return eigenVec;
}

VectorXd stdToEigenVec(std::vector<double> stdVec) {
	int n = stdVec.size();
	VectorXd eigenVec(n);
	for (int i = 0; i < n; i++) {
		eigenVec << stdVec[i];
	}
	return eigenVec;
}