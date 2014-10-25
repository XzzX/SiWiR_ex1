#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>
#include	<immintrin.h>

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

#include	"Timer.h"
#include	"matrix.hpp"

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

double a[LD * LD] __attribute__((aligned(0x40)));/// KxL - input
double b[LD * LD] __attribute__((aligned(0x40)));/// LxM - input
double bT[LD * LD] __attribute__((aligned(0x40)));
double c[LD * LD] __attribute__((aligned(0x40)));/// KxM - output


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
	for (int m = 0; m < C.getDimM(); m+=2){				///rows of c
		for (int n = 0; n < C.getDimN(); n+=2){			///cols of c	
			//std::cout << m << "\t" << n << std::endl;
			__m256d*	pA = A.get(m, 0);
			__m256d*	pB = A.get(m+1, 0);
			__m256d*	pC = B.getT(0, n);
			__m256d*	pD = B.getT(0, n+1);
			//std::cout << pA << "\t" << pB << "\t" << pC << "\t" << pD << std::endl;
			__m256d		W = _mm256_setzero_pd();
			__m256d		X = _mm256_setzero_pd();
			__m256d		Y = _mm256_setzero_pd();
			__m256d		Z = _mm256_setzero_pd();
			for (int l = 0; l < A.getDimN(); l+=4){
				//std::cout <<"mul" << std::endl;
				W = W + (*pA) * (*pC);
				X = X + (*pA) * (*pD);
				Y = Y + (*pB) * (*pC);
				Z = Z + (*pB) * (*pD);
				//std::cout << "inc" <<std::endl;
				pA++;
				pB++;
				pC++;
				pD++;
				//pE++;
			}
			//__m256d s = {W[0]+W[1]+W[2]+W[3], X[0]+X[1]+X[2]+X[3], Y[0]+Y[1]+Y[2]+Y[3], Z[0]+Z[1]+Z[2]+Z[3]};

			C(m, n) = W[0] + W[1] + W[2] + W[3];
			C(m, n+1) = X[0] + X[1] + X[2] + X[3];
			C(m+1, n) = Y[0] + Y[1] + Y[2] + Y[3];
			C(m+1, n+1) = Z[0] + Z[1] + Z[2] + Z[3];
		}
	}
}

void	zeroMemory(double* data){
	for (int i = 0; i < LD; ++i)
		for (int j = 0; j < LD; ++j)
			data[i*LD + j] = 0;
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
	if (( dimM > LD ) || ( dimN > LD )) {
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
	double time = 0;
	
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

	transpose(B, BT);

	naive(A, B, C);

	time = timer.elapsed();
	std::cout << A.getDimM() << "\t" << A.getDimN() << "\t" << C.getDimN() << "\t" << time << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************

	saveMatrix(argv[3], C);
};
