#include "Model.h"
#include <vector>

void Model::calcLoadVector() {

	int dof;
	int eType;
	int eNum;
	VectorXi nodes;
	VectorXd b(nDof);
	MatrixXd coords;
	MatrixXd nodalLoads;//Each row contains the load [x,y,z] for the corresponding node

	for (int i = 0; i < loads.size(); i++) {
		if (loads[i,0] == 0) {
			dof = nodeIndices[loads[i, 1]] + loads[i, 2];
			b[dof] = loads[i, 3];
		}
		else {
			eNum = loads[i, 1];
			coords = this->getElementCoords(eNum);
			eType = elementTable[eNum];
			nodes = VectorXi(topology[eNum]);
			nodalLoads = elementTypes[eType].calcNodalLoads(coords, loads[i,2], loads[i,3]);
			
			for (int j = 0; j < nodalLoads.size[0]; j++) {
				dof = nodeIndices[nodes[j]];
				for (int k = 0; k < nodalLoads.size[1]; k++) {
					b[dof + k] = nodalLoads[j, k];
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

	this->nElem = topology.size[0];
	this->nNodes = nodes.size[0];

	this->rotNodes = new bool[nNodes];
	for (int i = 0; i < nNodes; i++) rotNodes[i] = false;

	for (int i = 0; i < nElem; i++) {
		for (int j = 0; j < topology[i].size(); j++){
			if (elementTypes[elementTable[i]].hasRotDof()) rotNodes[i] = true;
		}
	}

	this->nodeIndices = VectorXi();

	int index = 0;
	for (int i = 0; i < nNodes; i++) {
		nodeIndices[i] = index;
		if (rotNodes[i]) index += 6;
		else index += 3;
	}

	this->nDof = nodeIndices[nNodes-1] + 1;
}

void Model::assembleGlobalSystem() {
	this->Kglob = MatrixXd::Zero(nDof, nDof);
	std::vector<int> dofs, nodes;
	MatrixXd coords;
	int n;
	for (int i = 0; i < nElem; i++) {
		nodes = topology[i];
		n = nodes.size();
		for (int j = 0; j < n; j++) dofs.push_back(nodeIndices(nodes[j]));
		coords = getElementCoords(i);
		elementTypes[elementTable[i]].assembleToGlobal(dofs, coords, Kglob);
	}
}

MatrixXd Model::getElementCoords(int eNum) {
	std::vector<int> elementNodes = topology[eNum];
	int n = elementNodes.size();
	MatrixXd coords(n, 3);

	for (int i = 0; i < n; i++) {
		coords << nodes.row(elementNodes[i]);
	}
}
