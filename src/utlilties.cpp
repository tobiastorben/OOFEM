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

double calcL(MatrixXd coords) {
	VectorXd diff = coords.row(1) - coords.row(0);
	return diff.norm();
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows)
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

	matrix.conservativeResize(numRows, numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;

	if (colToRemove < numCols)
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

	matrix.conservativeResize(numRows, numCols);
}

void removeRow(Eigen::VectorXd& vector, unsigned int rowToRemove)
{
	int numRows = vector.rows();
	VectorXd tmp(numRows - 1);
	int j = 0;
	for (int i = 0; i < numRows; i++) {
		if (i != rowToRemove) {
			tmp(j) = vector(i);
			j++;
		}
	}

	vector = tmp;
}