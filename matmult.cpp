#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

#include	"Timer.h"
#include	"matrix.hpp"

static const int	SIZE = 2048;

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

double bT[SIZE * SIZE]; /// KxM - output

inline
void	transpose(const Matrix& M, Matrix& MT){
	//transpose b
	for (int m = 0; m < M.getDimM(); ++m){				///rows of b
		for (int n = 0; n < M.getDimN(); ++n){			///cols of b	
			MT(n, m) = M(m, n);
		}
	}
}

///calculates AxB=C
inline
void	naive(const Matrix& A, const Matrix& B, Matrix& C){
	for (int m = 0; m < C.getDimM(); m++){				///rows of c
		for (int n = 0; n < C.getDimN(); n++){			///cols of c	
			for (int l = 0; l < A.getDimN(); l++){
				C(m, n) += A(m, l) * B.T(l, n);
			}
		}
	}
}

void	zeroMemory(double* data){
	for (int i = 0; i < SIZE; ++i)
		for (int j = 0; j < SIZE; ++j)
			data[i*SIZE + j] = 0;
}

Matrix	loadMatrix(std::string filename, double* data){
	int dimM = 0;
	int dimN = 0;
	
	std::ifstream	fIn(filename);
	if (!fIn) {
		std::cout << "Error opening file: " << filename << std::endl;
		exit(EXIT_FAILURE);
	}

	if(!(fIn >> dimM >> dimN))
	{
		std::cout << "Error in reading matrix entries!" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (( dimM > SIZE ) || ( dimN > SIZE )) {
		std::cout << "Matrix too big!" << std::endl;
		exit(EXIT_FAILURE);
	}

	Matrix	temp(data, nullptr, dimM, dimN, 0, 0);

	zeroMemory(data);

	for (int m = 0; m < dimM; m++) {
		for (int n = 0; n < dimN; n++) {
			if(!(fIn >> temp(m, n)))
			{
				std::cout << "Error in reading matrix entries!" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	fIn.close();

	return temp;
}

void	saveMatrix(std::string filename, const Matrix& mat){
	std::ofstream	fOut(filename);
	if (!fOut) {
		std::cout << "Error opening file: " << filename << std::endl;
		exit(EXIT_FAILURE);
	}

	if(!(fOut << mat.getDimM() << " " << mat.getDimN() << std::endl))
	{
		std::cout << "Error in writing matrix entries!" << std::endl;
		exit(EXIT_FAILURE);
	}

	for (int m = 0; m < mat.getDimM(); m++) {
		for (int n = 0; n < mat.getDimN(); n++) {
			if(!(fOut << mat(m, n) << std::endl))
			{
				std::cout << "Error in writing matrix entries!" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	fOut.close();
}

int main(int argc, char **argv) {

	///******************************************************
	///********************** INPUT *************************
	///******************************************************

	if (argc != 4) {
		std::cout << "Invalid number of arguments!" << std::endl;
		std::cout << "./compare A.out B.out" << std::endl;
		exit(EXIT_FAILURE);
	}

	Matrix	A = loadMatrix(argv[1], &a[0]);
	Matrix	B = loadMatrix(argv[2], &b[0]);
	B.dataT = &bT[0];

	Matrix	C(&c[0], nullptr, A.getDimM(), B.getDimN(), 0, 0);

	Matrix	BT(&bT[0], nullptr, B.getDimN(), B.getDimM(), 0, 0);
	
	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************

#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

	transpose(B, BT);

	naive(A, B, C);

	std::cout << A.getDimM() << "\t" << A.getDimN() << "\t" << C.getDimN() << "\t" << timer.elapsed() << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************

	saveMatrix(argv[3], C);
};
