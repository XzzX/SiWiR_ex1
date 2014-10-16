#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>

#include	"Timer.h"

static const unsigned int	SIZE = 2048;
static const double			EPSILON = 1e-5;

/**
Converts a string to an arbitrary type. >> operator must be defined for the target type.
@param string string which should be converted
@return converted string
**/
template<typename T>
T StringTo(const std::string& string){
    T valor;

    std::stringstream stream(string);
    stream >> valor;
    return valor;
}

double a[SIZE * SIZE];
double b[SIZE * SIZE];

int main(int argc, char **argv) {

	if (argc != 3) {
		std::cout << "Invalid number of arguments!" << std::endl;
		std::cout << "./compare A.out B.out" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::ifstream	fIn(argv[1]);
	if (!fIn) {
		std::cout << "Error opening file: " << argv[1] << std::endl;
		exit(EXIT_FAILURE);
	}

	unsigned int	rowCount = 0;
	unsigned int	colCount = 0;

	fIn >> rowCount >> colCount;
	if (( rowCount > SIZE ) || ( colCount > SIZE )) {
		std::cout << "Matrix too big!" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = 0; i < rowCount; i++) {
		for (unsigned int j = 0; j < colCount; j++) {
			fIn >> a[i*rowCount + j];
		}
	}

	fIn.close();

	fIn.open(argv[2]);
	if (!fIn) {
		std::cout << "Error opening file: " << argv[2] << std::endl;
		exit(EXIT_FAILURE);
	}

	unsigned int	tempRowCount = 0;
	unsigned int	tempColCount = 0;

	fIn >> tempRowCount >> tempColCount;
	if (( tempRowCount != rowCount ) || ( tempColCount != colCount )) {
		std::cout << "Matrix size is not equal!" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = 0; i < rowCount; i++) {
		for (unsigned int j = 0; j < colCount; j++) {
			fIn >> b[i*rowCount + j];
		}
	}

	fIn.close();

	for (unsigned int i = 0; i < rowCount; i++) {
		for (unsigned int j = 0; j < colCount; j++) {
			if ((fabs(a[i*rowCount + j] - b[i*rowCount + j])) > EPSILON) {
				std::cout << "Matrices NOT equal!" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	std::cout << "Matrices ARE equal!" << std::endl;
	exit(EXIT_SUCCESS);
};