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
			dof = nodeIndices(loads(i, 1)) + loads(i, 2);
			b(dof) = loads(i, 3);
		}
		else {
			eNum = loads(i, 1);
			coords = this->getElementCoords(eNum);
			eType = elementTable(eNum);
			nodes = stdToEigenVec(topology[eNum]);
			nodalLoads = elementTypes[eType].calcNodalLoads(coords, loads(i,2), loads(i,3));
			
			for (int j = 0; j < nodalLoads.rows(); j++) {
				dof = nodeIndices(nodes[j]);
				for (int k = 0; k < nodalLoads.cols(); k++) {
					b(dof + k) = nodalLoads(j, k);
				}
			}
		}
	}

}

Model::Model(std::vector<std::vector<int>> someTopology, MatrixXd someNodes, Element* someEtypes, VectorXi someEtable, MatrixXd someLoads, MatrixXi someSupports) {
	this->topology = someTopology;
	this->nodes = someNodes;
	this->elementTypes = someEtypes;
	this->elementTable = someEtable;
	this->loads = someLoads;
	this->supports = someSupports;

	this->nElem = topology.size();
	this->nNodes = nodes.rows();

	this->rotNodes = new bool[nNodes];
	for (int i = 0; i < nNodes; i++) rotNodes[i] = false;

	for (int i = 0; i < nElem; i++) {
		for (int j = 0; j < topology[i].size(); j++){
			if (elementTypes[elementTable(i)].hasRotDof()) rotNodes[i] = true;
		}
	}

	this->nodeIndices = VectorXi(nNodes);

	int index = 0;
	for (int i = 0; i < nNodes; i++) {
		nodeIndices(i) = index;
		if (rotNodes[i]) index += 6;
		else index += 3;
	}

	this->nDof  = index;
}

void Model::assembleGlobalSystem() {
	this->Kglob = MatrixXd::Zero(nDof, nDof);
	std::vector<int> dofs, eNodes;
	MatrixXd coords;
	int n;
	for (int i = 0; i < nElem; i++) {
		eNodes = topology[i];
		n = eNodes.size();
		for (int j = 0; j < n; j++) dofs.push_back(nodeIndices(eNodes[j]));
		coords = getElementCoords(i);
		elementTypes[elementTable(i)].assembleToGlobal(dofs, coords, Kglob);
		dofs.clear();
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
	int dof;
	for (int i = 0; i < supports.rows(); i++) {
		dof = nodeIndices(supports(i, 0)) + supports(i, 1);
		for (int j = 0; j < nDof; j++) {
			constrainedKglob(dof, j) = 0;
			constrainedKglob(j, dof) = 0;
		}
		constrainedKglob(dof, dof) = 1;
	}
}

void Model::solveSystem() {
	this->u = constrainedKglob.colPivHouseholderQr().solve(b);
}

VectorXd Model::getU() { return u; }

void Model::printModel() {

	std::cout << "\n\nLoad vector\n\n" << b << std::endl;
	std::cout << "\n\nStiffness matrix\n\n" << Kglob << std::endl;
	std::cout << "\n\nConstrained stiffness matrix\n\n" << constrainedKglob << std::endl;
	std::cin.ignore();
}