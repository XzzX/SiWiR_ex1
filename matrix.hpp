#ifndef	MATRIX_HPP_INCLUDED
#define	MATRIX_HPP_INCLUDED

#include	<immintrin.h>

int	LD = 2144;

class	Matrix{
public:
	double*		data;
	double*		dataT;
	
	//matrix dimensions
	int	dimRows;
	int	dimCols;

	//access offset
	int	oRows;
	int	oCols;

	Matrix(double*	dat, double* datT, const int m, const int n, const int offsetM, const int offsetN) : data(dat), dataT(datT), dimRows(m), dimCols(n), oRows(offsetM), oCols(offsetN){
	}

	
	inline
	double& operator()(const int row, const int col) {
		return data[(row+oRows) * LD + col + oCols];
	}
	
	inline
	const double& operator()(const int row, const int col) const {
		return data[(row + oRows) * LD + col + oCols];
	}
	
	inline
	double& T(const int row, const int col) {
		return dataT[(col+oCols) * LD + row + oRows];
	}
	
	inline
	const double& T(const int row, const int col) const {
		return dataT[(col + oCols) * LD + row + oRows];
	}
	
	inline
	void	set(const int row, const int col, const __m256d& v){
		_mm256_store_pd(&data[(row + oRows) * LD + col + oCols], v);
	}

	inline
	void 	setT(const int row, const int col, const __m256d& v){
		_mm256_store_pd(&dataT[(col + oCols) * LD + row + oRows], v);
	}

	inline
	__m256d* get(const int row, const int col) const {
		return (__m256d*)(&data[(row + oRows) * LD + col + oCols]);
	}

	inline
	__m256d* getT(const int row, const int col) const {
		return (__m256d*)(&dataT[(col + oCols) * LD + row + oRows]);
	}

	inline
	Matrix getSubMatrix(const int row, const int col, const int dimR, const int dimC){
		return Matrix(data, dataT, dimR, dimC, row+oRows, col + oCols);
	}

	inline
	int getDimM() const {
		return dimRows;
	}
	inline
	int getDimN() const {
		return dimCols;
	}
};

#endif
