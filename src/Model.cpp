#include "Model.h"
#include <vector>
#include <Eigen/Dense>

void Model::calcLoadVector() {

	int dof;
	int eType;
	int eNum;
	VectorXi nodes;
	this->b = VectorXd::Zero(nDof);
	MatrixXd coords;
	MatrixXd nodalLoads;//Each row contains the load [x,y,z] for the corresponding node

	for (int i = 0; i < loads.rows(); i++) {
		if (loads(i,0) == 0) {
			dof = 6*loads(i, 1) + loads(i, 2);
			b(dof) = loads(i, 3);
		}
		else {
			eNum = loads(i, 1);
			coords = this->getElementCoords(eNum);
			eType = elementTable(eNum);
			nodes = stdToEigenVec(topology[eNum]);
			nodalLoads = elementTypes[eType].calcNodalLoads(coords, loads(i,2), loads(i,3));
			
			for (int j = 0; j < nodalLoads.rows(); j++) {
				dof = 6*nodes[j];
				for (int k = 0; k < nodalLoads.cols(); k++) {
					b(dof + k) = nodalLoads(j, k);
				}
			}
		}
	}

}

Model::Model(std::vector<std::vector<int>> someTopology, MatrixXd someNodes, Element* someEtypes, VectorXi someEtable, MatrixXd someLoads, MatrixXi someSupports, int dim) {
	this->topology = someTopology;
	this->nodes = someNodes;
	this->elementTypes = someEtypes;
	this->elementTable = someEtable;
	this->loads = someLoads;
	this->supports = someSupports;
	this->dimensionality = dim;

	this->nElem = topology.size();
	this->nNodes = nodes.rows();

	this->nDof = 6 * nNodes;
}

void Model::assembleGlobalSystem() {
	this->Kglob = MatrixXd::Zero(nDof, nDof);
	std::vector<int> eNodes;
	MatrixXd coords;
	int n;
	for (int i = 0; i < nElem; i++) {
		eNodes = topology[i];
		coords = getElementCoords(i);
		elementTypes[elementTable(i)].assembleToGlobal(eNodes, coords, Kglob);
	}
}

MatrixXd Model::getElementCoords(int eNum) {
	std::vector<int> elementNodes = topology[eNum];
	int n = elementNodes.size();
	MatrixXd coords(n, 3);

	for (int i = 0; i < n; i++) {
		coords.row(i) = nodes.row(elementNodes[i]);
	}

	return coords;
}

void Model::applyBCs() {
	this->constrainedKglob = Kglob;
	this->constrainedB = b;
	int dof;
	for (int i = 0; i < supports.rows(); i++) {
		dof = 6*supports(i, 0) + supports(i, 1);
		for (int j = 0; j < nDof; j++) {
			constrainedKglob(dof, j) = 0;
			constrainedKglob(j, dof) = 0;
		}
		constrainedKglob(dof, dof) = 1;
		constrainedB(dof) = 0;
	}
}

void Model::solveSystem() {
	auto QRdecomposition = constrainedKglob.colPivHouseholderQr();
	if (QRdecomposition.isInvertible()) {
		this-> u = QRdecomposition.solve(constrainedB);
	}

	else {
		std::cerr << "ERROR: SYSTEM IS NOT INVERTIBLE! ABORTING" << std::endl;
		std::cin.ignore();
		abort();
	}
}

VectorXd Model::getU() { return u; }

void Model::printModel() {

	std::cout << "\n\nLoad vector\n\n" << b << std::endl;
	std::cout << "\n\nStiffness matrix\n\n" << Kglob << std::endl;
	std::cout << "\n\nConstrained stiffness matrix\n\n" << constrainedKglob << std::endl;
	std::cin.ignore();
}

VectorXd Model::calcReactions() {
	return Kglob*u - b;
}

VectorXd Model::calcError() {
	return constrainedKglob*u - b;
}

void Model::reduceSystem() {

	for (int i = 0; i < nDof; i++) {
		if (constrainedKglob.row(i).isZero() && constrainedKglob.col(i).isZero()) {
			constrainedKglob(i, i) = 1;
		}
	}

}

void Model::applyDimensionality() {//TODO : Take rotdof into account
	for (int i = 0; i < nDof;) {
		if (dimensionality == 2) {
			removeRow(Kglob, i+2);
			removeColumn(Kglob, i+2);
			removeRow(b, i + 2);
			nDof -= 1;
			i += 2;
		}
		if (dimensionality == 1) {
			removeRow(Kglob, i + 1);
			removeColumn(Kglob, i + 1);
			removeRow(b, i + 1);
			removeRow(Kglob, i + 2);
			removeColumn(Kglob, i + 2);
			removeRow(b, i + 2);
			nDof -= 2;
			i += 1;
		}
	}
}