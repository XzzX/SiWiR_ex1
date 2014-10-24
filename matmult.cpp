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

double a[SIZE * SIZE] __attribute__((aligned(32))); /// KxL - input
double b[SIZE * SIZE] __attribute__((aligned(32))); /// LxM - input
double c[SIZE * SIZE] __attribute__((aligned(32))); /// KxM - output

double bT[SIZE * SIZE] __attribute__((aligned(32))); /// KxM - output

void	zeroMemory(double* data){
	for (unsigned int i = 0; i < SIZE; ++i)
		for (unsigned int j = 0; j < SIZE; ++j)
			data[i*SIZE + j] = 0;
}

int main(int argc, char **argv) {
	zeroMemory(a);
	zeroMemory(b);
	zeroMemory(c);
	zeroMemory(bT);

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

	for (unsigned int k = 0; k < dimK; k++) {
		for (unsigned int l = 0; l < dimL; l++) {
			fIn >> a[k*dimL + l];
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

	for (unsigned int l = 0; l < dimL; l++) {
		for (unsigned int m = 0; m < dimM; m++) {
			fIn >> b[l*dimM + m];
		}
	}

	fIn.close();
	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************

	

#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;
	
	//transpose b
	for (unsigned int l = 0; l < dimL; l++){				///rows of b
		for (unsigned int m = 0; m < dimM; m++){			///cols of b	
			bT[m * dimL + l] = b[l * dimM + m];
		}
	}

	for (size_t k = 0; k < dimK; k++){				///rows of c
		for (size_t m = 0; m < dimM; m+=4){			///cols of c		
			__m256d*	pA = (__m256d*) &a[k*dimL];
			__m256d*	pB = (__m256d*) &bT[m * dimL];
			__m256d*	pC = (__m256d*) &bT[m * dimL + 1];
			__m256d*	pD = (__m256d*) &bT[m * dimL + 2];
			__m256d*	pE = (__m256d*) &bT[m * dimL + 3];
			//__m256d		pC;
			__m256d		W = _mm256_setzero_pd();
			__m256d		X = _mm256_setzero_pd();
			__m256d		Y = _mm256_setzero_pd();
			__m256d		Z = _mm256_setzero_pd();
			for (size_t l = 0; l < dimL; l += 4){
				//pC = _mm256_mul_pd(*pA, *pB);
				//pD = _mm256_add_pd(pC, pD);
				W = W + (*pA) * (*pB);
				X = X + (*pA) * (*pC);
				Y = Y + (*pA) * (*pD);
				Z = Z + (*pA) * (*pE);
				pA++;
				pB++;
				pC++;
				pD++;
				pE++;
			}
			__m256d s = {W[0]+W[1]+W[2]+W[3], X[0]+X[1]+X[2]+X[3], Y[0]+Y[1]+Y[2]+Y[3], Z[0]+Z[1]+Z[2]+Z[3]};

			//std::cout << k << "\t" << m << std::endl;
			// 0 972
			_mm256_store_pd(&c[k * dimM + m], s);
		}
	}
		std::cout << "calculation" << std::endl;
	/*
	double temp[4] __attribute__((aligned(32)));
	
	//calculate product
	for (size_t k = 0; k < dimK; k++){				///rows of c
		for (size_t m = 0; m < dimM; m++){			///cols of c	
			__m256d*	pA = (__m256d*) &a[k*dimL];
			__m256d*	pB = (__m256d*) &bT[m * dimL];
			__m256d		pC;
			__m256d		pD = _mm256_setzero_pd();
			for (size_t l = 0; l < dimL; l += 4){*/
				/*for (size_t t = 0; t < 2; t++){
					temp[t] += a[k*dimL + l + t] * bT[m * dimL + l + t];
				}*/
				/*pC = _mm256_mul_pd(*pA, *pB);
				pD = _mm256_add_pd(pC, pD);
				pA++;
				pB++;
			}
			_mm256_store_pd(&temp[0], pD);
			c[k * dimM + m] = temp[0] + temp[1] + temp[2] + temp[3];
		}
	}*/

	std::cout << dimK << "\t" << dimL << "\t" << dimM << "\t" << timer.elapsed() << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
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

	for (unsigned int k = 0; k < dimK; k++) {
		for (unsigned int m = 0; m < dimM; m++) {
			fOut << c[k*dimM + m] << std::endl;
		}
	}

	fOut.close();
};
