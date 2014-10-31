#ifndef	MATRIX_HPP_INCLUDED
#define	MATRIX_HPP_INCLUDED

#include	<immintrin.h>

///global leading dimension of all matrices
int	LD = 2144;

/*
Matrix access object. Not responsible for memory management. Just handles access. Possibility for easy block matrix creation.
*/
class	Matrix{
public:
	///raw data for ROW MAJOR matrix
	double*		data;
	///raw data for COLUMN MAJOR matrix
	double*		dataT;
	
	///number of rows
	int	dimRows;
	///number of cols
	int	dimCols;

	///access offset for rows
	int	oRows;
	///access offset for cols
	int	oCols;

	///basic constructor
	Matrix(double*	dat, double* datT, const int m, const int n, const int offsetM, const int offsetN) : data(dat), dataT(datT), dimRows(m), dimCols(n), oRows(offsetM), oCols(offsetN){
	}

	/*
	Normal access operator for double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	double& operator()(const int row, const int col) {
		return data[(row+oRows) * LD + col + oCols];
	}
	
	/*
	Normal access operator for double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	const double& operator()(const int row, const int col) const {
		return data[(row + oRows) * LD + col + oCols];
	}
	
	/*
	Access operation for column major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	double& T(const int row, const int col) {
		return dataT[(col+oCols) * LD + row + oRows];
	}
	
	/*
	Access operation for column major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	const double& T(const int row, const int col) const {
		return dataT[(col + oCols) * LD + row + oRows];
	}
	
	/*
	AVX write operation for row major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	void	set(const int row, const int col, const __m256d& v){
		_mm256_store_pd(&data[(row + oRows) * LD + col + oCols], v);
	}

	/*
	AVX write operation for column major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	void 	setT(const int row, const int col, const __m256d& v){
		_mm256_store_pd(&dataT[(col + oCols) * LD + row + oRows], v);
	}

	/*
	AVX read operation for row major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	__m256d* get(const int row, const int col) const {
		return (__m256d*)(&data[(row + oRows) * LD + col + oCols]);
	}

	/*
	AVX read operation for column major double access

	@param row Row
	@param col Col
	@return return value
	*/
	inline
	__m256d* getT(const int row, const int col) const {
		return (__m256d*)(&dataT[(col + oCols) * LD + row + oRows]);
	}

	/*
	Creates access object for sub block. No data is copied!

	@param row starting row
	@param col starting col
	@param dimR row dimension
	@param dimC col dimension
	@return block matrix object
	*/
	inline
	Matrix getSubMatrix(const int row, const int col, const int dimR, const int dimC){
		return Matrix(data, dataT, dimR, dimC, row+oRows, col + oCols);
	}
	///returns row dimension
	inline
	int getDimM() const {
		return dimRows;
	}
	///returns col dimension
	inline
	int getDimN() const {
		return dimCols;
	}
};

#endif
