#include <Eigen/Core>

using namespace Eigen;

MatrixXd getCoords(int eNum, MatrixXd topology, MatrixXd nodes);
Vector3d calcUnitVec(MatrixXd coords);