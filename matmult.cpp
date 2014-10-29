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

static const int	BLOCK_SIZE = 128; ///multiple of 4!!!!!!!!!!!!!!!!!!!!!!!

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
	int dimM = M.getDimM();
	int dimN = M.getDimN();
	//transpose b

	for (int m = 0; m < dimM; m+=BLOCK_SIZE){				///rows of b
		for (int n = 0; n < dimN; n+=4){			///cols of b	
			for (int i = m; i<m+BLOCK_SIZE; ++i){
				__m256d*	pM = M.get(i, n);
				MT(n, i) = (*pM)[0];
				MT(n+1, i) = (*pM)[1];
				MT(n+2, i) = (*pM)[2];
				MT(n+3, i) = (*pM)[3];
				pM++;
				
			}
		}
	}
}

///calculates AxB=C
inline
void	naive(const Matrix& A, const Matrix& B, Matrix& C){
	int dimM = C.getDimM();
	int dimN = C.getDimN();
	int dimL = A.getDimN();
	
	for (int m = 0; m < dimM; m+=4){				///rows of c
		for (int n = 0; n < dimN; n+=4){			///cols of c	
			//std::cout << m << "\t" << n << std::endl;
			__m256d*	pA = A.get(m, 0);
			__m256d*	pB = A.get(m+1, 0);
			__m256d*	pC = A.get(m+2, 0);
			__m256d*	pD = A.get(m+3, 0);
			__m256d*	pK = B.getT(0, n);
			__m256d*	pL = B.getT(0, n+1);
			__m256d*	pM = B.getT(0, n+2);
			__m256d*	pN = B.getT(0, n+3);
			//std::cout << pA << "\t" << pB << "\t" << pC << "\t" << pD << std::endl;
			__m256d		K = _mm256_setzero_pd();
			__m256d		L = _mm256_setzero_pd();
			__m256d		M = _mm256_setzero_pd();
			__m256d		N = _mm256_setzero_pd();
			__m256d		O = _mm256_setzero_pd();
			__m256d		P = _mm256_setzero_pd();
			__m256d		Q = _mm256_setzero_pd();
			__m256d		R = _mm256_setzero_pd();
			__m256d		S = _mm256_setzero_pd();
			__m256d		T = _mm256_setzero_pd();
			__m256d		U = _mm256_setzero_pd();
			__m256d		V = _mm256_setzero_pd();
			__m256d		W = _mm256_setzero_pd();
			__m256d		X = _mm256_setzero_pd();
			__m256d		Y = _mm256_setzero_pd();
			__m256d		Z = _mm256_setzero_pd();
			for (int l = 0; l < dimL; l+=4){
				//std::cout <<"mul" << std::endl;
				K = K + (*pA) * (*pK);
				L = L + (*pA) * (*pL);
				M = M + (*pA) * (*pM);
				N = N + (*pA) * (*pN);
				O = O + (*pB) * (*pK);
				P = P + (*pB) * (*pL);
				Q = Q + (*pB) * (*pM);
				R = R + (*pB) * (*pN);
				S = S + (*pC) * (*pK);
				T = T + (*pC) * (*pL);
				U = U + (*pC) * (*pM);
				V = V + (*pC) * (*pN);
				W = W + (*pD) * (*pK);
				X = X + (*pD) * (*pL);
				Y = Y + (*pD) * (*pM);
				Z = Z + (*pD) * (*pN);
				//std::cout << "inc" <<std::endl;
				pA++;
				pB++;
				pC++;
				pD++;
				pK++;
				pL++;
				pM++;
				pN++;
			}
			C(m  , n)     = K[0] + K[1] + K[2] + K[3];
			C(m  , n+1)   = L[0] + L[1] + L[2] + L[3];
			C(m  , n+2)   = M[0] + M[1] + M[2] + M[3];
			C(m  , n+3)   = N[0] + N[1] + N[2] + N[3];
			C(m+1, n  )   = O[0] + O[1] + O[2] + O[3];
			C(m+1, n+1)   = P[0] + P[1] + P[2] + P[3];
			C(m+1, n+2)   = Q[0] + Q[1] + Q[2] + Q[3];
			C(m+1, n+3)   = R[0] + R[1] + R[2] + R[3];
			C(m+2, n  )   = S[0] + S[1] + S[2] + S[3];
			C(m+2, n+1)   = T[0] + T[1] + T[2] + T[3];
			C(m+2, n+2)   = U[0] + U[1] + U[2] + U[3];
			C(m+2, n+3)   = V[0] + V[1] + V[2] + V[3];
			C(m+3, n  )   = W[0] + W[1] + W[2] + W[3];
			C(m+3, n+1)   = X[0] + X[1] + X[2] + X[3];
			C(m+3, n+2)   = Y[0] + Y[1] + Y[2] + Y[3];
			C(m+3, n+3)   = Z[0] + Z[1] + Z[2] + Z[3];
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
