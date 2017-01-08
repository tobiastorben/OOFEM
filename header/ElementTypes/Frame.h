#include "Element.h"
#include <Eigen/Core>
//2-noded, 12DOF Timoschenko beam element, free from shear locking, based on this model: http://www.m-hikari.com/atam/atam2008/atam1-4-2008/luoATAM1-4-2008-1
//Local coordinate system is through area center or shear center
//Local y-axis is defined to always align with the global XY-plane
class Frame : public Element {

public:
	Frame(double someE, double someA, double someIy, double someIz,double someG, double someK);
	//VectorXd calcStress(VectorXd u);
	//VectorXd calcStrain(VectorXd u);
	void assembleToGlobal(std::vector<int> eNodes, MatrixXd coords, MatrixXd& Kglob);
	MatrixXd calcNodalLoads(MatrixXd coords, int dof, double intensity);

private:
	double E;
	double A;
	double Iy;
	double Iz;
	double G;
	double k;


};