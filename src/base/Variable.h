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

class NcFileVector : protected std::vector<NcFile *> {

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
	///		Size of this NcFileVector.
	///	</summary>
	size_type size() const {
		return std::vector<NcFile *>::size();
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
	///		Populate from a vector of file names.
	///	</summary>
	void InsertFile(
		const std::string & strFile
	) {
		NcFile * pNewFile = new NcFile(strFile.c_str());
		if (pNewFile == NULL) {
			_EXCEPTIONT("Unable to allocate new NcFile");
		}
		if (!pNewFile->is_valid()) {
			_EXCEPTION1("Cannot open input file \"%s\"", strFile.c_str());
		}
		push_back(pNewFile);
		m_vecFilenames.push_back(strFile);
	}

	///	<summary>
	///		Parse from a semi-colon delineated string of file names.
	///	</summary>
	void ParseFromString(
		const std::string & strFiles,
		bool fAppend = true
	) {
		if (!fAppend) {
			clear();
		}

		int iLast = 0;
		for (int i = 0; i <= strFiles.length(); i++) {
			if ((i == strFiles.length()) ||
			    (strFiles[i] == ';')
			) {
				std::string strFile =
					strFiles.substr(iLast, i - iLast);

				InsertFile(strFile);

				iLast = i+1;
			}
		}

		if (size() == 0) {
			_EXCEPTION1("No input files found in \"%s\"",
				strFiles.c_str());
		}
		_ASSERT(size() == m_vecFilenames.size());
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

	///	<summary>
	///		Get the filename at the specified position.
	///	</summary>
	const std::string & GetFilename(size_t pos) const {
		_ASSERT(pos < m_vecFilenames.size());
		return m_vecFilenames[pos];
	}

public:
	///	<summary>
	///		Accessor.
	///	</summary>
	NcFile * operator[](size_type pos) const {
		return at(pos);
	}

protected:
	///	<summary>
	///		Vector of file names.
	///	</summary>
	std::vector<std::string> m_vecFilenames;
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

///	<summary>
///		A structure containing both a dimension name and size.
///	</summary>
class DimInfo {
public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	DimInfo() : name(), size(0) { }

	///	<summary>
	///		Constructor.
	///	</summary>
	DimInfo(
		const std::string & _name,
		long _size
	) :
		name(_name),
		size(_size)
	{ }

	///	<summary>
	///		Comparator (needed to have sets of DimInfo).
	///	</summary>
	bool operator< (const DimInfo & di) const {
		if (name < di.name) {
			return true;
		}
		if (size < di.size) {
			return true;
		}
		return false;
	}

public:
	///	<summary>
	///		Dimension name.
	///	</summary>
	std::string name;

	///	<summary>
	///		Dimension size.
	///	</summary>
	long size;
};

///	<summary>
///		A vector of DimInfo.
///	</summary>
typedef std::vector<DimInfo> DimInfoVector;

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
				+ std::to_string((long long)m_vecValue[d])
				+ std::string("/")
				+ std::to_string((long long)m_vecSize[d])
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
	///		Maximum number of arguments in variable.
	///	</summary>
	static const int MaxVariableArguments = 10;

public:
	///	<summary>
	///		Constructor (build an empty registry).
	///	</summary>
	VariableRegistry();

	///	<summary>
	///		Destructor.
	///	</summary>
	~VariableRegistry();

protected:
	///	<summary>
	///		Find this Variable in the registry.  If it does exist then
	///		return the corresponding VariableIndex and delete pvar.  If it
	///		does not exist then insert it into the registry.
	///	</summary>
	void InsertUniqueOrDelete(
		Variable * pvar,
		VariableIndex * pvarix
	);

public:
	///	<summary>
	///		Parse the given recursively string and register all relevant
	///		Variables.  At the end of the Variable definition return
	///		the current position in the string.
	///	</summary>
	int FindOrRegisterSubStr(
		const std::string & strIn,
		VariableIndex * pvarix
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
	std::string GetVariableString(VariableIndex varix) const;

	///	<summary>
	///		Unload all data.
	///	</summary>
	void UnloadAllGridData();

protected:
	///	<summary>
	///		Get the list of base variable indices.
	///	</summary>
	void GetDependentVariableIndicesRecurse(
		VariableIndex varix,
		std::vector<VariableIndex> & vecDependentIxs
	) const;

public:
	///	<summary>
	///		Get the list of base variable indices.
	///	</summary>
	void GetDependentVariableIndices(
		VariableIndex varix,
		std::vector<VariableIndex> & vecDependentIxs
	) const;

	///	<summary>
	///		Get the list of base variable names.
	///	</summary>
	void GetDependentVariableNames(
		VariableIndex varix,
		std::vector<std::string> & vecDependentVarNames
	) const;

public:
	///	<summary>
	///		Obtain auxiliary dimension sizes for the specified variables.
	///	</summary>
	static void GetAuxiliaryDimInfo(
		const NcFileVector & vecncDataFiles,
		const SimpleGrid & grid,
		const std::string & strVarName,
		DimInfoVector & vecAuxDimInfo
	);

	///	<summary>
	///		Obtain auxiliary dimension sizes for the specified variables,
	///		and verify consistency among all variables.
	///	</summary>
	static void GetAuxiliaryDimInfoAndVerifyConsistency(
		const NcFileVector & vecncDataFiles,
		const SimpleGrid & grid,
		const std::vector<std::string> & vecVariables,
		DimInfoVector & vecAuxDimInfo
	);

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
///		A map between auxiliary indices and the data stored therein.
///	</summary>
typedef std::map<VariableAuxIndex, DataArray1D<float> *> DataMap;

///	<summary>
///		A class for storing a 2D slice of data, and associated metadata.
///	</summary>
class Variable {

friend class VariableRegistry;

public:
	///	<summary>
	///		An invalid time index.
	///	</summary>
	static const long InvalidTimeIndex;

	///	<summary>
	///		No time index.
	///	</summary>
	static const long NoTimeIndex;

	///	<summary>
	///		Invalid argument.
	///	</summary>
	static const long InvalidArgument;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Variable() :
		m_strName(),
		m_fOp(false),
		m_fNoTimeInNcFile(false),
		m_lTime(InvalidTimeIndex)
	{ }

public:
	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Variable & var);

	///	<summary>
	///		Get the name of this Variable.
	///	</summary>
	const std::string & GetName() const {
		return m_strName;
	}

	///	<summary>
	///		Check if this Variable is an operator.
	///	</summary>
	bool IsOp() const {
		return m_fOp;
	}

	///	<summary>
	///		Get the number of specified arguments.
	///	</summary>
	size_t GetArgumentCount() const {
		_ASSERT(m_strArg.size() == m_varArg.size());
		_ASSERT(m_strArg.size() == m_lArg.size());
		return m_strArg.size();
	}

	///	<summary>
	///		Get the array of specified operator argument strings.
	///	</summary>
	const std::vector<std::string> & GetArgumentStrings() const {
		return m_strArg;
	}

	///	<summary>
	///		Get the array of specified operator arguments (as long).
	///	</summary>
	const std::vector<long> & GetArgumentLongs() const {
		return m_lArg;
	}

	///	<summary>
	///		Get the array of specified operator Variables.
	///	</summary>
	const VariableIndexVector & GetArgumentVarIxs() const {
		return m_varArg;
	}

public:
	///	<summary>
	///		Get a string representation of this variable.
	///	</summary>
	std::string ToString(
		const VariableRegistry & varreg
	) const;

protected:
	///	<summary>
	///		Get the first instance of this variable in the given NcFileVector.
	///	</summary>
	NcVar * GetFromNetCDF(
		VariableRegistry & varreg,
		const NcFileVector & vecFiles,
		const SimpleGrid & grid,
		long lTime = (-1)
	);

public:
	///	<summary>
	///		Load a data block from the NcFileVector.
	///	</summary>
	void LoadGridData(
		VariableRegistry & varreg,
		const NcFileVector & vecFiles,
		const SimpleGrid & grid,
		long lTime = (-1)
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

	///	<summary>
	///		Flag indicating this is an operator.
	///	</summary>
	bool m_fOp;

	///	<summary>
	///		Specified operator arguments (as std::string).
	///	</summary>
	std::vector<std::string> m_strArg;

	///	<summary>
	///		Specified operator arguments (as long, for NetCDF file indexing).
	///	</summary>
	std::vector<long> m_lArg;

	///	<summary>
	///		Specified operator arguments (as VariableIndex).
	///	</summary>
	VariableIndexVector m_varArg;

protected:
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

public:
	///	<summary>
	///		Flag indicating this Variable has no time index in NetCDF file.
	///	</summary>
	bool m_fNoTimeInNcFile;

	///	<summary>
	///		Time index associated with data loaded in this Variable.
	///	</summary>
	long m_lTime;

	///	<summary>
	///		Data associated with this Variable.
	///	</summary>
	DataArray1D<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

