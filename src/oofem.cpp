#include <Eigen/Core>
#include <iostream>

int main(int argc, char** argv) {

	Eigen::MatrixXd m(3, 3);
	m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	std::cout << m;
	std::cin.ignore();
}