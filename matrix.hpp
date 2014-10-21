class	Matrix{
public:
	double*		data;
	double*		dataT;
	
	//matrix dimensions
	int	dimM;
	int	dimN;

	//dimensions of data lines
	int	ldm;
	int	ldn;

	//access offset
	int	oM;
	int	oN;

	Matrix(double*	dat, double* datT, const int m, const int n, const int offsetM, const int offsetN) : data(dat), dataT(datT), dimM(m), dimN(n), ldm(m), ldn(n), oM(offsetM), oN(offsetN){
	}

	inline
	double& operator()(int row, int col){
		return data[row * ldn + col];
	}

	inline
	const double& operator()(int row, int col) const {
		return data[row * ldn + col];
	}

	inline
	double& T(int row, int col){
		return dataT[col * ldm + row];
	}

	inline
	const double& T(int row, int col) const {
		return dataT[col * ldm + row];
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
