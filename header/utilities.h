#ifndef UTILITIES_H
#define UTILITIES_H

#include <Eigen/Core>
#include <vector>
#include <math.h>

using namespace Eigen;

Vector3d calcUnitVec(MatrixXd coords);
VectorXi stdToEigenVec(std::vector<int> stdVec);
VectorXd stdToEigenVec(std::vector<double> stdVec);
double calcL(MatrixXd coords);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
void removeRow(Eigen::VectorXd& vector, unsigned int rowToRemove);

#endif