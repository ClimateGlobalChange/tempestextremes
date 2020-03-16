///////////////////////////////////////////////////////////////////////////////
///
///	\file    Variable.h
///	\author  Paul Ullrich
///	\version November 18, 2015
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _VARIABLE_H_
#define _VARIABLE_H_

#include "netcdfcpp.h"

#include "DataArray1D.h"
#include "SimpleGrid.h"
#include "DataOp.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

class NcFileVector : public std::vector<NcFile *> {

private:
	///	<summary>
	///		Copy constructor.
	///	</summary>
	NcFileVector(const NcFileVector & vecFiles) {
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	NcFileVector & operator= (const NcFileVector & vecFiles) {
		return (*this);
	}

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	NcFileVector() :
		std::vector<NcFile *>()
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~NcFileVector() {
		NcFileVector::clear();
	}

	///	<summary>
	///		Clear the contents of this NcFileVector.
	///	</summary>
	void clear() {
		for (size_t i = 0; i < size(); i++) {
			(*this)[i]->close();
			delete (*this)[i];
		}
	}

	///	<summary>
	///		Parse from a semi-colon delineated string of file names.
	///	</summary>
	void ParseFromString(
		const std::string & strFiles
	) {
		int iLast = 0;
		for (int i = 0; i <= strFiles.length(); i++) {
			if ((i == strFiles.length()) ||
			    (strFiles[i] == ';')
			) {
				std::string strFile =
					strFiles.substr(iLast, i - iLast);

				NcFile * pNewFile = new NcFile(strFile.c_str());

				if (!pNewFile->is_valid()) {
					_EXCEPTION1("Cannot open input file \"%s\"",
						strFile.c_str());
				}

				push_back(pNewFile);
				iLast = i+1;
			}
		}

		if (size() == 0) {
			_EXCEPTION1("No input files found in \"%s\"",
				strFiles.c_str());
		}
	}

	///	<summary>
	///		Get an iterator to the file containing the specified variable.
	///	</summary>
	const_iterator FindContainingVariable(
		const std::string & strVariable
	) const {
		for (const_iterator iter = begin(); iter != end(); iter++) {
			NcVar * var = (*iter)->get_var(strVariable.c_str());
			if (var == NULL) {
				return iter;
			}
		}
		return end();
	}
};

///////////////////////////////////////////////////////////////////////////////

class Variable;

///	<summary>
///		A vector of pointers to Variable.
///	</summary>
typedef std::vector<Variable *> VariableVector;

///	<summary>
///		A unique index assigned to each Variable.
///	</summary>
typedef int VariableIndex;

///	<summary>
///		The invalid variable index.
///	</summary>
static const VariableIndex InvalidVariableIndex = (-1);

///	<summary>
///		A vector for storing VariableIndex.
///	</summary>
class VariableIndexVector : public std::vector<VariableIndex> {};

///////////////////////////////////////////////////////////////////////////////

class VariableRegistry {

public:
	///	<summary>
	///		Constructor (build an empty registry).
	///	</summary>
	VariableRegistry();

	///	<summary>
	///		Destructor.
	///	</summary>
	~VariableRegistry();

public:
/*
	///	<summary>
	///		Register a Variable.  Or return an index if the Variable already
	///		exists in the registry.
	///	</summary>
	VariableIndex FindOrRegister(
		const Variable & var
	);
*/
	///	<summary>
	///		Parse the given recursively string and register all relevant
	///		Variables.  At the end of the Variable definition return
	///		the current position in the string.
	///	</summary>
	VariableIndex FindOrRegisterSubStr(
		const std::string & strIn,
		int * piFinalStringPos
	);

	///	<summary>
	///		Parse the given string recursively and register all
	///		relevant Variables.
	///	</summary>
	VariableIndex FindOrRegister(
		const std::string & strIn
	);

	///	<summary>
	///		Get the Variable with the specified index.
	///	</summary>
	Variable & Get(VariableIndex varix);

	///	<summary>
	///		Get the descriptor for the Variable with the specified index.
	///	</summary>
	std::string GetVariableString(VariableIndex varix);

	///	<summary>
	///		Unload all data.
	///	</summary>
	void UnloadAllGridData();

public:
	///	<summary>
	///		Get the DataOp with the specified name.
	///	</summary>
	DataOp * GetDataOp(const std::string & strName);

private:
	///	<summary>
	///		Array of variables.
	///	</summary>
	VariableVector m_vecVariables;

	///	<summary>
	///		Map of data operators.
	///	</summary>
	DataOpManager m_domDataOp;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An object holding dimension indices for a given Variable.
///	</summary>
typedef std::vector<long> VariableDimIndex;

///	<summary>
///		An object holding auxiliary indices for a given Variable.
///	</summary>
typedef VariableDimIndex VariableAuxIndex;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class describing an iterator over the auxiliary indices of a
///		Variable.
///	</summary>
class VariableAuxIndexIterator {

friend class Variable;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableAuxIndexIterator() :
		m_fEnd(false)
	{ }

protected:
	///	<summary>
	///		Initializer, only accessible from Variable.
	///	</summary>
	void Initialize(
		const VariableAuxIndex & vecSize,
		bool fEnd
	) {
		m_fEnd = fEnd;
		m_vecSize = vecSize;
		m_vecValue.resize(m_vecSize.size());
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			if (m_vecSize[d] <= 0) {
				_EXCEPTIONT("Invalid auxiliary index size entry");
			}
		}
		if ((fEnd) && (vecSize.size() > 0)) {
			m_vecValue[0] = m_vecSize[0];
		}
	}

public:
	///	<summary>
	///		Remove the dimension with specified index.
	///	</summary>
	void RemoveDim(
		size_t dim
	) {
		if (dim >= m_vecSize.size()) {
			_EXCEPTIONT("Logic errror");
		}
		m_vecSize.erase(m_vecSize.begin() + dim);
		m_vecValue.erase(m_vecValue.begin() + dim);
	}

public:
	///	<summary>
	///		Prefix incrementor.
	///	</summary>
	VariableAuxIndexIterator & operator++() {
		if (m_vecValue.size() == 0) {
			if (m_fEnd) {
				_EXCEPTIONT("Iterator exceeded bounds");
			} else {
				m_fEnd = true;
			}
			return (*this);
		}
		for (size_t d = 0; d < m_vecValue.size(); d++) {
			size_t dx = m_vecValue.size()-d-1;
			if (m_vecValue[dx] >= m_vecSize[dx]) {
				_EXCEPTIONT("Iterator exceeded bounds");
			} else if (m_vecValue[dx] == m_vecSize[dx]-1) {
				m_vecValue[dx] = 0;
			} else {
				m_vecValue[dx]++;
				return (*this);
			}
		}
		m_fEnd = true;
		return (*this);
	}

	///	<summary>
	///		Postfix incrementor.
	///	</summary>
	VariableAuxIndexIterator operator++(int) {
		this->operator++();
		return (*this);
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator==(const VariableAuxIndexIterator & iter) const {
		if (m_vecSize.size() != iter.m_vecSize.size()) {
			_EXCEPTIONT("Invalid comparison");
		}
		if ((m_fEnd) && (iter.m_fEnd)) {
			return true;
		}
		if (m_fEnd != iter.m_fEnd) {
			return false;
		}
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			if (m_vecValue[d] != iter.m_vecValue[d]) {
				return false;
			}
		}
		return true;
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator!=(const VariableAuxIndexIterator & iter) const {
		return !((*this)==(iter));
	}

	///	<summary>
	///		Cast to VariableAuxIndex.
	///	</summary>
	operator VariableAuxIndex() const {
		return m_vecValue;
	}

	///	<summary>
	///		Get the value of the iterator.
	///	</summary>
	const VariableAuxIndex & Value() const {
		return m_vecValue;
	}

	///	<summary>
	///		Convert to a std::string.
	///	</summary>
	std::string ToString() const {
		std::string strOut;
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			strOut +=
				std::string("(")
				+ std::to_string(m_vecValue[d])
				+ std::string("/")
				+ std::to_string(m_vecSize[d])
				+ std::string(")");
		}
		if (m_fEnd) {
			strOut += std::string("(end)");
		} else if (m_vecSize.size() == 0) {
			strOut += std::string("(begin)");
		}

		return strOut;
	}

	///	<summary>
	///		Get the length of this auxiliary index.
	///	</summary>
	size_t size() const {
		return m_vecSize.size();
	}

protected:
	///	<summary>
	///		The sizes of the auxiliary indices.
	///	</summary>
	VariableAuxIndex m_vecSize;

	///	<summary>
	///		The values of the auxiliary indices.
	///	</summary>
	VariableAuxIndex m_vecValue;

	///	<summary>
	///		At the end.
	///	</summary>
	bool m_fEnd;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A map between auxiliary indices and the data stored therein.
///	</summary>
typedef std::map<VariableAuxIndex, DataArray1D<float> *> DataMap;

///	<summary>
///		A class storing a parsed variable name.
///	</summary>
class Variable {

friend class VariableRegistry;

public:
	///	<summary>
	///		Maximum number of arguments in variable.
	///	</summary>
	static const int MaxArguments = 4;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Variable() :
		m_strName(),
		m_fOp(false),
		m_nSpecifiedDim(0),
		m_fNoTimeInNcFile(false),
		m_iTime(-2)
	{
		memset(m_iDim, 0, MaxArguments * sizeof(int));
	}

public:
	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Variable & var);

public:
	///	<summary>
	///		Get a string representation of this variable.
	///	</summary>
	std::string ToString(
		VariableRegistry & varreg
	) const;

	///	<summary>
	///		Get this variable in the given NcFile.
	///	</summary>
	NcVar * GetFromNetCDF(
		NcFileVector & vecFiles,
		int iTime = (-1)
	);

	///	<summary>
	///		Load a data block from the NcFileVector.
	///	</summary>
	void LoadGridData(
		VariableRegistry & varreg,
		NcFileVector & vecFiles,
		const SimpleGrid & grid,
		int iTime = (-1)
	);

	///	<summary>
	///		Unload the current data block.
	///	</summary>
	void UnloadGridData();

	///	<summary>
	///		Get the data associated with this variable.
	///	</summary>
	const DataArray1D<float> & GetData() const {
		return m_data;
	}

	///	<summary>
	///		Get the data associated with this variable.
	///	</summary>
	DataArray1D<float> & GetData() {
		return m_data;
	}

protected:
	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strName;
/*
	///	<summary>
	///		Variable units.
	///	</summary>
	std::string m_strUnits;

	///	<summary>
	///		Auxiliary dimensions associated with this variable.
	///	</summary>
	std::vector<std::string> m_vecAuxDimNames;

	///	<summary>
	///		Index of the record dimension.
	///	</summary>
	int m_iTimeDimIx;

	///	<summary>
	///		Index of the vertical dimension.
	///	</summary>
	int m_iVerticalDimIx;

	///	<summary>
	///		(+1) if the vertical coordinate is bottom-up, (-1) if top-down.
	///	</summary>
	int m_nVerticalDimOrder;
*/
protected:
	///	<summary>
	///		Flag indicating this is an operator.
	///	</summary>
	bool m_fOp;

	///	<summary>
	///		Number of dimensions specified.
	///	</summary>
	int m_nSpecifiedDim;

	///	<summary>
	///		Specified dimension values.
	///	</summary>
	int m_iDim[MaxArguments];

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	std::vector<std::string> m_strArg;

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	VariableIndexVector m_varArg;

public:
	///	<summary>
	///		Flag indicating this Variable has no time index in NetCDF file.
	///	</summary>
	bool m_fNoTimeInNcFile;

	///	<summary>
	///		Time index associated with data loaded in this Variable.
	///	</summary>
	int m_iTime;

	///	<summary>
	///		Data associated with this Variable.
	///	</summary>
	DataArray1D<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

