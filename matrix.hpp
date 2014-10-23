#ifndef	MATRIX_HPP_INCLUDED
#define	MATRIX_HPP_INCLUDED

static const int	LDM = 2144;
static const int	LDN = 2144;

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
		return data[row * LDN + col];
	}
	
	inline
	const double& operator()(const int row, const int col) const {
		return data[row * LDN + col];
	}
	
	inline
	void	set(const int row, const int col, const double& v){
		data[row * LDN + col] = v;
	}

	inline
	void 	setT(const int row, const int col, const double& v){
		dataT[col * LDM + row] = v;
	}

	inline
	const double& get(const int row, const int col) const {
		return data[row * LDN + col];
	}

	inline
	const double& getT(const int row, const int col) const {
		return dataT[col * LDM + row];
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
