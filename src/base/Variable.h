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

#include "DataVector.h"
#include "SimpleGrid.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<NcFile *> NcFileVector;

///////////////////////////////////////////////////////////////////////////////

class Variable;

typedef std::vector<Variable> VariableVector;

typedef int VariableIndex;

typedef std::vector<VariableIndex> VariableIndexVector;

///////////////////////////////////////////////////////////////////////////////

class VariableRegistry {

public:
	///	<summary>
	///		Register a variable.  Or return an index if the Variable already
	///		exists in the registry.
	///	</summary>
	int FindOrRegister(const Variable & var);

	///	<summary>
	///		Get the variable with the specified index.
	///	</summary>
	Variable & Get(VariableIndex varix);

	///	<summary>
	///		Unload all data.
	///	</summary>
	void UnloadAllGridData();

private:
	///	<summary>
	///		Array of variables.
	///	</summary>
	VariableVector m_vecVariables;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a parsed variable name.
///	</summary>
class Variable {

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
		m_fOp(false),
		m_strName(),
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
	///		Parse the variable information from a string.  Return the index
	///		of the first character after the variable information.
	///	</summary>
	int ParseFromString(
		VariableRegistry & varreg,
		const std::string & strIn
	);

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
	const DataVector<float> & GetData() const {
		return m_data;
	}

public:
	///	<summary>
	///		Flag indicating this is an operator.
	///	</summary>
	bool m_fOp;

	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strName;

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
	DataVector<float> m_data;
};

///////////////////////////////////////////////////////////////////////////////

#endif

