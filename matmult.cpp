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
double p[LD * LD] __attribute__((aligned(0x40)));/// KxM - output
double ps[LD * LD] __attribute__((aligned(0x40)));/// KxM - output
double s[LD * LD] __attribute__((aligned(0x40)));/// KxM - output
double t[LD * LD] __attribute__((aligned(0x40)));/// KxM - output


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

void	strassen(Matrix& A, Matrix& B, Matrix& C, Matrix& P, Matrix& Ps, Matrix& S, Matrix& T){
	int	dim = A.getDimM(); //equal for all matrix dimensions
	int	dim2 = dim * 0.5;
	//std::cout << dim << "\t" << dim2 << std::endl;	
	
	Matrix	A1 = A.getSubMatrix(0, 0, dim2, dim2);
	Matrix	A2 = A.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	A3 = A.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	A4 = A.getSubMatrix(dim2, dim2, dim2, dim2);
	
	Matrix	B1 = B.getSubMatrix(0, 0, dim2, dim2);
	Matrix	B2 = B.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	B3 = B.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	B4 = B.getSubMatrix(dim2, dim2, dim2, dim2);
	
	Matrix	C1 = C.getSubMatrix(0, 0, dim2, dim2);
	Matrix	C2 = C.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	C3 = C.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	C4 = C.getSubMatrix(dim2, dim2, dim2, dim2);
	
	Matrix	S1 = S.getSubMatrix(0, 0, dim2, dim2);
	Matrix	S2 = S.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	S3 = S.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	S4 = S.getSubMatrix(dim2, dim2, dim2, dim2);
	
	
	Matrix	T1 = T.getSubMatrix(0, 0, dim2, dim2);
	Matrix	T2 = T.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	T3 = T.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	T4 = T.getSubMatrix(dim2, dim2, dim2, dim2);
		
	Matrix	P1 = P.getSubMatrix(0, 0, dim2, dim2);
	Matrix	P2 = P.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	P3 = P.getSubMatrix(dim2, 0, dim2, dim2);
	Matrix	P4 = P.getSubMatrix(dim2, dim2, dim2, dim2);

	Matrix	P5 = Ps.getSubMatrix(0, 0, dim2, dim2);
	Matrix	P6 = Ps.getSubMatrix(0, dim2, dim2, dim2);
	Matrix	P7 = Ps.getSubMatrix(dim2, 0, dim2, dim2);
	//Matrix	P8 = Ps.getSubMatrix(dim2, dim2, dim2, dim2);
	//std::cout << "submatrices" << std::endl;
	
	//compute temporary S and T matrices
	for (int i=0; i<dim2; ++i){
		for (int j=0; j<dim2; j+=4){
			//std::cout << "S" << std::endl;
			__m256d*	pA1 = A1.get(i, j);
			__m256d*	pA2 = A2.get(i, j);
			__m256d*	pA3 = A3.get(i, j);
			__m256d*	pA4 = A4.get(i, j);
			__m256d*	pS2 = S2.get(i, j);
			S1.set(i, j, (*pA3)+(*pA4));
			S2.set(i, j, (*pA3)+(*pA4)-(*pA1));
			S3.set(i, j, (*pA1)-(*pA3));
			S4.set(i, j, (*pA2)-(*pS2));
			pA1++;
			pA2++;
			pA3++;
			pA4++;
			pS2++;
			//S1(i, j) = A3(i, j) + A4(i, j);
			//S2(i, j) = S1(i, j) - A1(i, j);
			//S3(i, j) = A1(i ,j) - A3(i, j);
			//S4(i, j) = A2(i, j) - S2(i, j);
			//std::cout << "T" << std::endl;
			__m256d*	pB1 = B1.getT(j, i);
			__m256d*	pB2 = B2.getT(j, i);
			__m256d*	pB3 = B3.getT(j, i);
			__m256d*	pB4 = B4.getT(j, i);
			__m256d*	pT2 = T2.getT(j, i);
			T1.setT(j, i, (*pB2)-(*pB1));
			T2.setT(j, i, (*pB4)-((*pB2)-(*pB1)));
			T3.setT(j, i, (*pB4)-(*pB2));
			T4.setT(j, i, (*pB3)-(*pT2));
			pB1++;
			pB2++;
			pB3++;
			pB4++;
			pT2++;
			//T1.T(j, i) = B2.T(j, i) - B1.T(j, i);
			//T2.T(j, i) = B4.T(j, i) - T1.T(j, i);
			//T3.T(j, i) = B4.T(j ,i) - B2.T(j, i);
			//T4.T(j, i) = B3.T(j, i) - T2.T(j, i);

		}
	
	}

	//calculate products
	if (dim2<600) {
		naive(A1, B1, P1);
		naive(A2, B3, P2);
		naive(S1, T1, P3);
		naive(S2, T2, P4);
		naive(S3, T3, P5);
		naive(S4, B4, P6);
		naive(A4, T4, P7);
	} else {
		strassen(A1, B1, P1, C1, C2, C3, B2);
		strassen(A2, B3, P2, C1, C2, C3, B2);
		strassen(S1, T1, P3, C1, C2, C3, B2);
		strassen(S2, T2, P4, C1, C2, C3, B2);
		strassen(S3, T3, P5, C1, C2, C3, B2);
		strassen(S4, B4, P6, C1, C2, C3, B2);
		strassen(A4, T4, P7, C1, C2, C3, B2);
	}
	
	//assemble final matrix
	for (int i=0; i<dim2; ++i){
		for (int j=0; j<dim2; j+=4){
			__m256d*	pP1 = P1.get(i, j);
			__m256d*	pP2 = P2.get(i, j);
			__m256d*	pP3 = P3.get(i, j);
			__m256d*	pP4 = P4.get(i, j);
			__m256d*	pP5 = P5.get(i, j);
			__m256d*	pP6 = P6.get(i, j);
			__m256d*	pP7 = P7.get(i, j);
			C1.set(i, j, (*pP1) + (*pP2));
			C3.set(i, j, (*pP1) + (*pP4) + (*pP5) + (*pP7));
			C4.set(i, j, (*pP1) + (*pP4) + (*pP5) + (*pP3));
			C2.set(i, j, (*pP1) + (*pP4) + (*pP3) + (*pP6));
			//C1(i, j) = P1(i, j) + P2(i, j);
			//C3(i, j) = P1(i, j) + P4(i, j) + P5(i, j) + P7(i, j);
			//C4(i, j) = P1(i, j) + P4(i, j) + P5(i, j) + P3(i, j);
			//C2(i, j) = P1(i, j) + P3(i, j) + P4(i, j) + P6(i, j);
			pP1++;
			pP2++;
			pP3++;
			pP4++;
			pP5++;
			pP6++;
			pP7++;
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

void	printMatrix(const Matrix& mat){
	for (int m = 0; m < mat.getDimM(); m++) {
		for (int n = 0; n < mat.getDimN(); n++) {
			std::cout << mat(m, n) << " ";
		}
		std::cout << std::endl;
	}
}

void	printMatrixT(const Matrix& mat){
	for (int m = 0; m < mat.getDimM(); m++) {
		for (int n = 0; n < mat.getDimN(); n++) {
			std::cout << mat.T(m, n) << " ";
		}
		std::cout << std::endl;
	}
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
	
	Matrix	P(&p[0], nullptr, A.getDimM(), A.getDimM(), 0, 0);
	Matrix	PS(&ps[0], nullptr, A.getDimM(), A.getDimM(), 0, 0);
	Matrix	S(&s[0], nullptr, A.getDimM(), A.getDimM(), 0, 0);
	Matrix	T(nullptr, &t[0], A.getDimM(), A.getDimM(), 0, 0);
	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;
	
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

	A.dimRows = 2048;

	transpose(B, BT);
	
	//naive(A, B, C);

	strassen(A, B, C, P, PS, S, T);

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
