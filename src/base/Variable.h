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
#include "NcFileVector.h"

#include <vector>
#include <string>
#include <limits>

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

	///	<summary>
	///		Equality comparator.
	///	</summary>
	bool operator== (const DimInfo & di) const {
		if ((name == di.name) && (size == di.size)) {
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
class DimInfoVector : public std::vector<DimInfo> {

public:
	///	<summary>
	///		Get the total size of this DimInfoVector.
	///	</summary>
	size_t GetTotalSize() const {
		size_t sTotalSize = 1;
		for (size_t i = 0; i < size(); i++) {
			sTotalSize *= static_cast<size_t>((*this)[i].size);
		}
		return sTotalSize;
	}

	///	<summary>
	///		Convert this to a string.
	///	</summary>
	std::string ToString() const {
		std::string str;
		for (size_t d = 0; d < size(); d++) {
			str += "[" + (*this)[d].name + "," + std::to_string((*this)[d].size) + "]";
		}
		if (str.length() == 0) {
			str = "[]";
		}
		return str;
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

///	<summary>
///		A class describing an iterator over the auxiliary indices of a
///		Variable.
///	</summary>
class VariableAuxIndexIterator {

friend class VariableRegistry;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	VariableAuxIndexIterator() :
		m_varix(InvalidVariableIndex),
		m_fEnd(false)
	{ }

protected:
	///	<summary>
	///		Initializer, only accessible from VariableRegistry.
	///	</summary>
	void Initialize(
		VariableIndex varix,
		const VariableAuxIndex & vecSize,
		bool fEnd
	) {
		m_varix = varix;
		m_fEnd = fEnd;
		m_vecSize = vecSize;
		m_vecValue.resize(m_vecSize.size(), 0);
		for (size_t d = 0; d < m_vecSize.size(); d++) {
			if (m_vecSize[d] <= 0) {
				_EXCEPTIONT("Invalid auxiliary index size entry");
			}
		}
		if ((fEnd) && (vecSize.size() > 0)) {
			m_vecValue[0] = m_vecSize[0];
		}
	}

	///	<summary>
	///		Reset the value to zero.
	///	</summary>
	void Reset() {
		for (size_t d = 0; d < m_vecValue.size(); d++) {
			m_vecValue[d] = 0;
		}
		m_fEnd = false;
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
	///		Get the flattened index of the iterator.
	///	</summary>
	size_t Offset() const {
		size_t sValue = 0;
		size_t sAccumulatedSize = 1;
		for (long d = static_cast<long>(m_vecValue.size())-1; d >= 0; d--) {
			sValue += m_vecValue[d] * sAccumulatedSize;
			sAccumulatedSize *= m_vecSize[d];
		}
		return sValue;
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

	///	<summary>
	///		Returns true if this iterator is at the end.
	///	</summary>
	bool at_end() const {
		return m_fEnd;
	}

protected:
	///	<summary>
	///		VariableIndex associated with this iterator.
	///	</summary>
	VariableIndex m_varix;

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

private:
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
	///		Obtain auxiliary dimension sizes for the specified variables.
	///	</summary>
	void GetAuxiliaryDimInfo(
		const NcFileVector & vecncDataFiles,
		const SimpleGrid & grid,
		VariableIndex varix,
		DimInfoVector & vecAuxDimInfo
	) const;

private:
	///	<summary>
	///		Recursively assign auxiliary indices to Variables.
	///	</summary>
	void AssignAuxiliaryIndicesRecursive(
		Variable & var,
		const std::vector<std::string> & vecArg
	);

public:
	///	<summary>
	///		Clear the processing queue.
	///	</summary>
	void ClearProcessingQueue();

	///	<summary>
	///		Add a new variable and its auxiliary indices to the processing
	///		queue for this VariableRegistry.
	///	</summary>
	void AppendVariableToProcessingQueue(
		const NcFileVector & vecncDataFiles,
		const SimpleGrid & grid,
		VariableIndex varix,
		DimInfoVector * pvecAuxDimInfo = NULL
	);

	///	<summary>
	///		Get the current position in the processing queue.
	///	</summary>
	size_t GetProcessingQueueVarPos() const;

	///	<summary>
	///		Get the current processing queue variable index.
	///	</summary>
	VariableIndex GetProcessingQueueVarIx() const;

	///	<summary>
	///		Get the variable currently in the processing queue.
	///	</summary>
	Variable & GetProcessingQueueVariable();

	///	<summary>
	///		Get the current processing queue auxiliary index.
	///	</summary>
	const VariableAuxIndex & GetProcessingQueueAuxIx() const;

	///	<summary>
	///		Get the current processing queue auxiliary dimension size.
	///	</summary>
	const VariableAuxIndex & GetProcessingQueueAuxSize() const;

	///	<summary>
	///		Get the current processing queue auxiliary index offset.
	///	</summary>
	size_t GetProcessingQueueOffset() const;

	///	<summary>
	///		Advance the processing queue.
	///	</summary>
	bool AdvanceProcessingQueue();

	///	<summary>
	///		Reset the processing queue.
	///	</summary>
	void ResetProcessingQueue();

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

private:
	///	<summary>
	///		Current variable index in the processing queue.
	///	</summary>
	size_t m_sProcessingQueueVarPos;

	///	<summary>
	///		Map of auxiliary index iterators.
	///	</summary>
	std::vector<VariableAuxIndexIterator> m_vecProcessingQueue;
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
	///		Default constructor.
	///	</summary>
	Variable() :
		m_strName(),
		m_fOp(false),
		m_dFillValueFloat(-std::numeric_limits<float>::max()),
		m_fNoTimeInNcFile(false),
		m_timeStored(Time::CalendarUnknown)
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
		return m_strArg.size();
	}

	///	<summary>
	///		Get the array of specified operator argument strings.
	///	</summary>
	const std::vector<std::string> & GetArgumentStrings() const {
		return m_strArg;
	}

	///	<summary>
	///		Get the array of specified operator Variables.
	///	</summary>
	const VariableIndexVector & GetArgumentVarIxs() const {
		return m_varArg;
	}

	///	<summary>
	///		Get the _FillValue for this variable.
	///	</summary>
	float GetFillValueFloat() const {
		return m_dFillValueFloat;
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
	NcVar * GetNcVarFromNcFileVector(
		const NcFileVector & ncfilevec,
		const SimpleGrid & grid
	);

public:
	///	<summary>
	///		Load a data block from the NcFileVector.
	///	</summary>
	void LoadGridData(
		VariableRegistry & varreg,
		const NcFileVector & vecFiles,
		const SimpleGrid & grid
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
	///		Free operator arguments.
	///	</summary>
	std::vector<bool> m_fFreeArg;

	///	<summary>
	///		Specified operator arguments or auxiliary indices (as std::string).
	///	</summary>
	std::vector<std::string> m_strArg;

	///	<summary>
	///		Specified operator arguments (as VariableIndex).
	///	</summary>
	VariableIndexVector m_varArg;

	///	<summary>
	///		_FillValue for this Variable.
	///	</summary>
	float m_dFillValueFloat;

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
	///		Time currently stored in this Variable.
	///	</summary>
	Time m_timeStored;

	///	<summary>
	///		Data associated with this Variable.
	///	</summary>
	DataArray1D<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

