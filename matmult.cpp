#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>

extern "C" {
#include	<cblas.h>
}

#ifdef USE_LIKWID
extern "C" {
#include	<likwid.h>
}
#endif

#include	"Timer.h"

static const unsigned int	SIZE = 2048;

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

double a[SIZE * SIZE]; /// KxL - input
double b[SIZE * SIZE]; /// LxM - input
double c[SIZE * SIZE]; /// KxM - output

int main(int argc, char **argv) {

	///******************************************************
	///********************** INPUT *************************
	///******************************************************

	if (argc != 4) {
		std::cout << "Invalid number of arguments!" << std::endl;
		std::cout << "./compare A.out B.out" << std::endl;
		exit(EXIT_FAILURE);
	}

	//dimension of matrices
	// KxL   LxM
	unsigned int	dimK = 0;
	unsigned int	dimL = 0;
	unsigned int	dimM = 0;

	std::ifstream	fIn(argv[1]);
	if (!fIn) {
		std::cout << "Error opening file: " << argv[1] << std::endl;
		exit(EXIT_FAILURE);
	}

	fIn >> dimK >> dimL;
	if (( dimK > SIZE ) || ( dimL > SIZE )) {
		std::cout << "Matrix too big!" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = 0; i < dimK; i++) {
		for (unsigned int j = 0; j < dimL; j++) {
			fIn >> a[i*dimK + j];
		}
	}

	fIn.close();

	fIn.open(argv[2]);
	if (!fIn) {
		std::cout << "Error opening file: " << argv[2] << std::endl;
		exit(EXIT_FAILURE);
	}

	unsigned int	tempDim = 0;

	fIn >> tempDim >> dimM;
	if ( tempDim != dimL ) {
		std::cout << "Matrix dimensions not correct!" << std::endl;
		exit(EXIT_FAILURE);
	}
	if ( dimM > SIZE ) {
		std::cout << "Matrix too big!" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned int i = 0; i < dimL; i++) {
		for (unsigned int j = 0; j < dimM; j++) {
			fIn >> b[i*dimL + j];
		}
	}

	fIn.close();

	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************

#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("blas");
#endif

	siwir::Timer	timer;

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dimK, dimM, dimL, 1.0, &a[0], dimL, &b[0], dimM, 0.0, &c[0], dimM);

	std::cout << "calculation took " << timer.elapsed() << "s" << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("blas");
	likwid_markerClose();
#endif

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************

	std::ofstream	fOut(argv[3]);
	if (!fOut) {
		std::cout << "Error opening file: " << argv[3] << std::endl;
		exit(EXIT_FAILURE);
	}

	fOut << dimK << " " << dimM << std::endl;

	for (unsigned int i = 0; i < dimK; i++) {
		for (unsigned int j = 0; j < dimM; j++) {
			fOut << c[i*dimL + j] << std::endl;
		}
	}

	fOut.close();
};