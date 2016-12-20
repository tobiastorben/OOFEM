#include <Eigen\Core>
#include <vector>
using namespace Eigen;
class Model {

private:
	MatrixXd nodes;
	std::vector<std::vector<char[]>> elementTypes;
	Element elements[];
	MatrixXd topology;
	MatrixXd supports;
	MatrixXd loads;

	MatrixXd Kglob;
	VectorXd b;


public:
	VectorXd solve(MatrixXd K, VectorXd b);
	


};