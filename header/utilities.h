#ifndef UTILITIES_H
#define UTILITIES_H

#include <Eigen/Core>
#include <vector>

using namespace Eigen;

Vector3d calcUnitVec(MatrixXd coords);
VectorXi stdToEigenVec(std::vector<int> stdVec);
VectorXd stdToEigenVec(std::vector<double> stdVec);

#endif