#ifndef	MATRIX_HPP_INCLUDED
#define	MATRIX_HPP_INCLUDED

#include	<immintrin.h>

static const int	LD = 2144;

class	Matrix{
public:
	double*		data;
	double*		dataT;
	
	//matrix dimensions
	int	dimM;
	int	dimN;

	//access offset
	int	oM;
	int	oN;

	Matrix(double*	dat, double* datT, const int m, const int n, const int offsetM, const int offsetN) : data(dat), dataT(datT), dimM(m), dimN(n), oM(offsetM), oN(offsetN){
	}

	
	inline
	double& operator()(const int row, const int col) {
		return data[row * LD + col];
	}
	
	inline
	const double& operator()(const int row, const int col) const {
		return data[row * LD + col];
	}
	
	inline
	void	set(const int row, const int col, const double& v){
		data[row * LD + col] = v;
	}

	inline
	void 	setT(const int row, const int col, const double& v){
		dataT[col * LD + row] = v;
	}

	inline
	__m256d* get(const int row, const int col) const {
		return (__m256d*)(&data[row * LD + col]);
	}

	inline
	__m256d* getT(const int row, const int col) const {
		return (__m256d*)(&dataT[col * LD + row]);
	}

	inline
	int getDimM() const {
		return dimM;
	}
	inline
	int getDimN() const {
		return dimN;
	}
};

#endif
