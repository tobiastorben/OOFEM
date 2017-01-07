#include "Frame.h"
#include "utilities.h"

Frame::Frame(double someE, double someA, double someIy, double someIz, double someJ, double someG, double someK) {
	E = someE;
	A = someA;
	Iy = someIy;
	Iz = someIz;
	J = someJ;
	G = someG;
	k = someK;
}

//VectorXd Frame::calcStress(VectorXd u) {}

//VectorXd Frame::calcStrain(VectorXd u) {}

void Frame::assembleToGlobal(std::vector<int> eNodes, MatrixXd coords, MatrixXd& Kglob) {
	//double entry;
	//int globIndI, globIndJ,i2,j2;
	//double L = calcL(coords);
	//double az = (12 * E*Iz - k*G*A*L*L)*(12 * E*Iz - k*G*A*L*L);
	//double ay = (12 * E*Iy - k*G*A*L*L)*(12 * E*Iy - k*G*A*L*L);
	//for (int i = 0; i < 12; i++) {
	//	for (int j = 0; j < 12; j++) {
	//		if ((i == 0 && j == 0) || (i == 6 && j == 6)) {
	//			entry = E*A / L;
	//		}
	//		else if (i == 0 && j == 6) {
	//			entry = -E*A / L;
	//		}
	//		else if ((i == 1 && j == 1) || (i == 7 && j == 7)) {
	//			entry = (12 * k*G*A*E*Iy*(12 * E*Iy + k*G*A*L*L)) / (L*ay);
	//		}
	//		else if (i == 1 && j == 7) {
	//			entry = -(12 * k*G*A*E*Iy*(12 * E*Iy + k*G*A*L*L)) / (L*ay);
	//		}
	//		else if (i == 1 && (j == 5 || j == 11)) {
	//			entry = (6 * k*G*A*E*Iy*(12 * E*Iy + k*G*A*L*L)) / (ay);
	//		}
	//		else if ((i == 2 && (j == 2 || j == 8)) || (i == 8 && j == 8)) {
	//			entry = (12 * k*G*A*E*Iz*(12 * E*Iz + k*G*A*L*L)) / (L*az);
	//		}
	//		else if (i == 2 && (j == 4 || j == 10)) {
	//			entry = -(6 * k*G*A*E*Iz*(12 * E*Iz + k*G*A*L*L)) / (az);
	//		}
	//		else if ((i == 3 && j == 3) || (i == 9 && j == 9)) {
	//			entry = G*(Iy + Iz) / L;
	//		}
	//		else if (i == 3 && j == 9) {
	//			entry = -G*(Iy + Iz) / L;
	//		}
	//		else if (i == 4 && j == 4 || (i == 10 && j == 10)) {
	//			entry = (4 * E*Iz*((k*G*A)*(k*G*A)*L*L*L*L + 3 * k*G*A*L*L*E*Iz + 36 * (E*E*Iz*Iz))) / (L*az);
	//		}
	//		else if ((i == 4 && j == 8) || (i == 8 && j == 10)) {
	//			entry = (6 * k*G*A*E*Iz*(12 * E*Iz + k*G*A*L*L)) / (az);
	//		}
	//		else if (i == 4 && j == 10) {
	//			entry = -(2 * E*Iz*(72 * E*E*Iz*Iz - k*k*G*G*A*A*L*L*L*L - 30 * k*G*A*L*L*E*Iz)) / (L*az);
	//		}
	//		else if ((i == 5 && j == 5) || (i == 11 && j == 11)) {
	//			entry = (4 * E*Iy*((k*G*A)*(k*G*A)*L*L*L*L + 3 * k*G*A*L*L*E*Iy + 36 * (E*E*Iy*Iy))) / (L*ay);
	//		}
	//		else if ((i == 5 && j == 7) || (i == 7 && j == 11)) {
	//			entry = -(6 * k*G*A*E*Iy*(12 * E*Iy + k*G*A*L*L)) / (ay);
	//		}
	//		else if (i == 5 && j == 11) {
	//			entry = -(2 * E*Iy*(72 * E*E*Iy*Iy - k*k*G*G*A*A*L*L*L*L - 30 * k*G*A*L*L*E*Iy)) / (L*ay);
	//		}

	//		else entry = 0;

	//		if (i < 6) {
	//			globIndI = 6*eNodes[0] + i;
	//		}
	//		else {
	//			globIndI = 6*eNodes[1] + i - 6;
	//		}

	//		if (j < 6) {
	//			globIndJ = 6 * eNodes[0] + j;
	//		}
	//		else {
	//			globIndJ = 6 * eNodes[1] + j - 6;
	//		}
	//		Kglob(globIndI, globIndJ) = Kglob(globIndI, globIndJ) + entry;
	//		i2 = j;
	//		j2 = i;
	//		if (i != j){
	//			if (i < 6) {
	//				globIndI = 6 * eNodes[0] + i;
	//			}
	//			else {
	//				globIndI = 6 * eNodes[1] + i - 6;
	//			}

	//			if (j < 6) {
	//				globIndJ = 6 * eNodes[0] + j;
	//			}
	//			else {
	//				globIndJ = 6 * eNodes[1] + j - 6;
	//			}
	//			Kglob(globIndI, globIndJ) = Kglob(globIndI, globIndJ) + entry;
	//		}

	//	}
	//}


	int globIndI, globIndJ;
	double L = calcL(coords);

	double K[12][12];

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			K[i][j] = 0;
		}
	}

	K[0][0] = (A*E) / L;
	K[0][6] = -(A*E) / L;
	K[1][1] = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K[1][5] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[1][7] = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K[1][11] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[2][2] = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K[2][4] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[2][8] = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K[2][10] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[3][3] = (G*k*(Iy + Iz)) / L;
	K[3][9] = -(G*k*(Iy + Iz)) / L;
	K[4][2] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[4][4] = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*4.0) / L;
	K[4][8] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[4][10] = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*2.0) / L;
	K[5][1] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[5][5] = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*4.0) / L;
	K[5][7] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[5][11] = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*2.0) / L;
	K[6][0] = -(A*E) / L;
	K[6][6] = (A*E) / L;
	K[7][1] = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K[7][5] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[7][7] = (A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K[7][11] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[8][2] = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-1.2E1) / L;
	K[8][4] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[8][8] = (A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*1.2E1) / L;
	K[8][10] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[9][3] = -(G*k*(Iy + Iz)) / L;
	K[9][9] = (G*k*(Iy + Iz)) / L;
	K[10][2] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[10][4] = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*2.0) / L;
	K[10][8] = A*E*G*Iz*k*(E*Iz*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[10][10] = -(1.0 / pow(E*Iz*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iz*Iz*Iz)*4.32E2 - A*(E*E)*G*(Iz*Iz)*(L*L)*k*1.08E2)) / L + (E*Iz*4.0) / L;
	K[11][1] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*6.0;
	K[11][5] = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*2.0) / L;
	K[11][7] = A*E*G*Iy*k*(E*Iy*1.2E1 + A*G*(L*L)*k)*1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*-6.0;
	K[11][11] = -(1.0 / pow(E*Iy*1.2E1 - A*G*(L*L)*k, 2.0)*((E*E*E)*(Iy*Iy*Iy)*4.32E2 - A*(E*E)*G*(Iy*Iy)*(L*L)*k*1.08E2)) / L + (E*Iy*4.0) / L;

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
			Kglob(globIndI, globIndJ) = Kglob(globIndI, globIndJ) + K[i][j];
		}
	}
}

MatrixXd Frame::calcNodalLoads(MatrixXd coords, int dof, double intensity){ return MatrixXd::Zero(1, 3); }
