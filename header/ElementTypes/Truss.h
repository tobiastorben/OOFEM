#include "Element.h"
#include <Eigen/Core>
#include <vector>

class Truss : public Element {
	private :
		double EA;
		double L;
		
	public:
		Truss(double someEA, double someL);
		//void calcStiffness();
		VectorXd calcStress(VectorXd u);
		VectorXd calcStrain(VectorXd u);
		void assembleToGlobal(std::vector<int> dofs, MatrixXd coords, , MatrixXd& Kglob);
		MatrixXd calcNodalLoads(MatrixXd coords, int dof, double intensity);
};