#include <Eigen/Core>

using namespace Eigen;

class Element {

protected:
	//MatrixXd K;
	bool rotDof;

	
public:
	//virtual void calcStiffness() = 0;
	virtual VectorXd calcStress(VectorXd u) = 0;
	virtual VectorXd calcStrain(VectorXd u) = 0;
	virtual void assembleToGlobal(std::vector<int> dofs, MatrixXd coords, MatrixXd& Kglob) = 0;
	virtual MatrixXd calcNodalLoads(MatrixXd coords, int dof, double intensity) = 0;//Return a matrix where each row contains the forcevector [0-5] for each node
	bool hasRotDof();
};