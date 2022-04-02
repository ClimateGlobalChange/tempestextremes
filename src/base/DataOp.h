///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataOp.h
///	\author  Paul Ullrich
///	\version July 22, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _DATAOP_H_
#define _DATAOP_H_

#include "DataArray1D.h"
#include "SparseMatrix.h"

#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class VariableRegistry;
class SimpleGrid;
class VariableIndexVector;

///////////////////////////////////////////////////////////////////////////////

class DataOp;

///////////////////////////////////////////////////////////////////////////////

class DataOpManager : protected std::map<std::string, DataOp*> {

protected:
	typedef std::map<std::string, DataOp *> DataOpMap;
	typedef DataOpMap::iterator DataOpMapIterator;
	typedef DataOpMap::value_type DataOpMapPair;

public:
	///	<summary>
	///		Destructor.
	///	</summary>
	~DataOpManager();

public:
	///	<summary>
	///		Add a new DataOp.
	///	</summary>
	DataOp * Add(DataOp * pdo);

	///	<summary>
	///		Attempt to dynamically generate an operator from the given name.
	///	</summary>
	DataOp * Add(const std::string & strName);

	///	<summary>
	///		Find a DataOp.
	///	</summary>
	DataOp * Find(const std::string & strName);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp() :
		m_strName("")
	{ }

	///	<summary>
	///		Constructor with name.
	///	</summary>
	DataOp(
		const std::string & strName
	) :
		m_strName(strName)
	{ }

public:
	///	<summary>
	///		Get the name of this operator.
	///	</summary>
	const std::string & GetName() {
		return m_strName;
	}

	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Name of this DataOp.
	///	</summary>
	std::string m_strName;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_VECMAG : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_VECMAG() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_ABS : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_ABS() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_SIGN : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_SIGN() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_ALLPOS : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_ALLPOS() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_SUM : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_SUM() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_AVG : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_AVG() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_DIFF : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_DIFF() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_PROD : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_PROD() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_DIV : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_DIV() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_MIN : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_MIN() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_MAX : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_MAX() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_COND : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_COND() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_SQRT : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_SQRT() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_POW : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_POW() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_LAT : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_LAT() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_F : public DataOp {

public:
	///	<summary>
	///		Operartor name.
	///	</summary>
	static const char * name;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_F() :
		DataOp(name)
	{ }

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_LAPLACIAN : public DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_LAPLACIAN(
		const std::string & strName,
		int nLaplacianPoints,
		double dLaplacianDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Number of points in this Laplacian.
	///	</summary>
	int m_nLaplacianPoints;

	///	<summary>
	///		Evaluation distance for the Laplacian.
	///	</summary>
	double m_dLaplacianDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator.
	///	</summary>
	SparseMatrix<float> m_opLaplacian;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_CURL : public DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_CURL(
		const std::string & strName,
		int nCurlPoints,
		double dCurlDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Number of points in this curl operator.
	///	</summary>
	int m_nCurlPoints;

	///	<summary>
	///		Evaluation distance for the curl operator.
	///	</summary>
	double m_dCurlDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator (eastward component).
	///	</summary>
	SparseMatrix<float> m_opCurlE;

	///	<summary>
	///		Sparse matrix operator (northward component).
	///	</summary>
	SparseMatrix<float> m_opCurlN;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_DIVERGENCE : public DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_DIVERGENCE(
		const std::string & strName,
		int nDivPoints,
		double dDivDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Number of points in this curl operator.
	///	</summary>
	int m_nDivPoints;

	///	<summary>
	///		Evaluation distance for the curl operator.
	///	</summary>
	double m_dDivDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator (eastward component).
	///	</summary>
	SparseMatrix<float> m_opDivE;

	///	<summary>
	///		Sparse matrix operator (northward component).
	///	</summary>
	SparseMatrix<float> m_opDivN;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_GRADMAG : public DataOp {

public:
	///	<summary>
	///		A pair that supports in-place addition and subtraction.
	///	</summary>
	template<typename T>
	class pair_with_plus_minus : public std::pair<T,T> {
		public:
			pair_with_plus_minus(
				const T & a_first,
				const T & a_second
			) :
				std::pair<T,T>(a_first, a_second)
			{ }

			pair_with_plus_minus(
				const T & a_val
			) :
				std::pair<T,T>(a_val, a_val)
			{ }

			pair_with_plus_minus<T> & operator+=(const pair_with_plus_minus<T> & pr) {
				std::pair<T,T>::first += pr.first;
				std::pair<T,T>::second += pr.second;
				return (*this);
			}

			pair_with_plus_minus<T> & operator-=(const pair_with_plus_minus<T> & pr) {
				std::pair<T,T>::first -= pr.first;
				std::pair<T,T>::second -= pr.second;
				return (*this);
			}

			pair_with_plus_minus<T> & operator=(const T & val) {
				std::pair<T,T>::first = val;
				std::pair<T,T>::second = val;
				return (*this);
			}
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_GRADMAG(
		const std::string & strName,
		int nGradPoints,
		double dGradDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Number of points in this curl operator.
	///	</summary>
	int m_nGradPoints;

	///	<summary>
	///		Evaluation distance for the curl operator.
	///	</summary>
	double m_dGradDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator (eastward component).
	///	</summary>
	SparseMatrix< pair_with_plus_minus<float> > m_opGrad;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_VECDOTGRAD : public DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_VECDOTGRAD(
		const std::string & strName,
		int nGradPoints,
		double dGradDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Number of points in this curl operator.
	///	</summary>
	int m_nGradPoints;

	///	<summary>
	///		Evaluation distance for the curl operator.
	///	</summary>
	double m_dGradDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator (eastward component).
	///	</summary>
	SparseMatrix< DataOp_GRADMAG::pair_with_plus_minus<float> > m_opGrad;
};

///////////////////////////////////////////////////////////////////////////////

class DataOp_MEAN : public DataOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DataOp_MEAN(
		const std::string & strName,
		double dMeanDist
	);

public:
	///	<summary>
	///		Apply the operator.
	///	</summary>
	virtual bool Apply(
		const SimpleGrid & grid,
		const std::vector<std::string> & strArg,
		const std::vector<DataArray1D<float> const *> & vecArgData,
		DataArray1D<float> & dataout
	);

protected:
	///	<summary>
	///		Evaluation distance for the mean operator.
	///	</summary>
	double m_dMeanDist;

	///	<summary>
	///		Flag indicating the sparse matrix operator is initialized.
	///	</summary>
	bool m_fInitialized;

	///	<summary>
	///		Sparse matrix operator.
	///	</summary>
	SparseMatrix<float> m_opMean;
};

///////////////////////////////////////////////////////////////////////////////

#endif

