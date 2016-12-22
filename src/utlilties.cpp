#include "utilities.h"

MatrixXd getCoords(int eNum, MatrixXd topology, MatrixXd nodes) {
	
}

Vector3d calcUnitVec(MatrixXd coords) {
	return (coords.row(1) - coords.row(0)).normalize;
}