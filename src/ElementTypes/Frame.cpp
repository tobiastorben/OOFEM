#include "Frame.h"
#include "utilities.h"

Frame::Frame(double someE, double someA, double someIy, double someIz, double someG, double someK) {
	E = someE;
	A = someA;
	Iy = someIy;
	Iz = someIz;
	G = someG;
	k = someK;
}

//VectorXd Frame::calcStress(VectorXd u) {}

//VectorXd Frame::calcStrain(VectorXd u) {}

void Frame::assembleToGlobal(std::vector<int> eNodes, MatrixXd coords, MatrixXd& Kglob) {
	
	int globIndI, globIndJ;
	double L = calcL(coords);

	//Coordinate transformation matrix
	double lx = (coords(1, 0) - coords(0, 0))/L;
	double mx = (coords(1, 1) - coords(0, 1))/L;
	double nx = (coords(1, 2) - coords(0, 2))/L;
	double Lxy = sqrt(lx*lx + mx*mx);

	MatrixXd R(3, 3);

	if (lx == 0 && mx == 0 && nx == 1) {
		R << 0, 0, -1,
			0, 1, 0,
			1, 0, 0;
	}

	else if (lx == 0 && mx == 0 && nx == -1) {
		R << 0, 0, 1,
			 0, 1, 0,
			-1, 0, 0;
	}

	else {
		R << lx, mx, nx,
			-mx / Lxy, lx / Lxy, 0,
			-lx*nx / Lxy, -mx*nx / Lxy, Lxy;
	}

	MatrixXd T = MatrixXd::Zero(12, 12);
	T.block(0, 0, 3, 3) = R;
	T.block(3, 3, 3, 3) = R;
	T.block(6, 6, 3, 3) = R;
	T.block(9, 9, 3, 3) = R;

	//Element stiffness matrix in local coordinates
	MatrixXd K = MatrixXd::Zero(12, 12);
	K(0, 0) = (A*E) / L;
	K(0, 6) = -(A*E) / L;
	K(1, 1) = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K(1, 5) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(1, 7) = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K(1, 11) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(2, 2) = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K(2, 4) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(2, 8) = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K(2, 10) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(3, 3) = (G*k*(Iy + Iz)) / L;
	K(3, 9) = -(G*k*(Iy + Iz)) / L;
	K(4, 2) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(4, 4) = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*4.0) / L;
	K(4, 8) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(4, 10) = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*2.0) / L;
	K(5, 1) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(5, 5) = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*4.0) / L;
	K(5, 7) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(5, 11) = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*2.0) / L;
	K(6, 0) = -(A*E) / L;
	K(6, 6) = (A*E) / L;
	K(7, 1) = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K(7, 5) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(7, 7) = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K(7, 11) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(8, 2) = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K(8, 4) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(8, 8) = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K(8, 10) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(9, 3) = -(G*k*(Iy + Iz)) / L;
	K(9, 9) = (G*k*(Iy + Iz)) / L;
	K(10, 2) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(10, 4) = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*2.0) / L;
	K(10, 8) = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(10, 10) = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*4.0) / L;
	K(11, 1) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K(11, 5) = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*2.0) / L;
	K(11, 7) = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K(11, 11) = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*4.0) / L;

	//Element stiffness matrix in global coordinates
	MatrixXd Kg = T.transpose()*K*T;

	//Assemlby to global system
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			if (i < 6) {
				globIndI = 6 * eNodes[0] + i;
			}
			else {
				globIndI = 6 * eNodes[1] + i - 6;
			}

			if (j < 6) {
				globIndJ = 6 * eNodes[0] + j;
			}
			else {
				globIndJ = 6 * eNodes[1] + j - 6;
			}
			Kglob(globIndI, globIndJ) = Kglob(globIndI, globIndJ) + Kg(i,j);
		}
	}
}

MatrixXd Frame::calcNodalLoads(MatrixXd coords, int dof, double intensity){ return MatrixXd::Zero(1, 3); }
