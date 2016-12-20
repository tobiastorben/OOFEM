#include <Eigen/Core>

using namespace Eigen;

class Element {

private:
	MatrixXd K;
	
public:
	void calcStiffness();
	VectorXd calcStress(VectorXd u);
	VectorXd calcStrain(VectorXd u);
	void assembleToGlobal(VectorXd nodes);
	void calcNodalLoads(VectorXd nodes, VectorXd load);



};