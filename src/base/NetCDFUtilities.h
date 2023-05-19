///////////////////////////////////////////////////////////////////////////////
///
///	\file    NetCDFUtilities.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#ifndef _NETCDFUTILITIES_H_
#define _NETCDFUTILITIES_H_

#include "netcdfcpp.h"
#include "TimeObj.h"

#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class NcTimeDimension : public std::vector<Time> {

public:
	///	<summary>
	///		Type of time dimension.
	///	</summary>
	enum TimeDimType {
		TimeDimType_Standard,
		TimeDimType_DailyMean,
		TimeDimType_MonthlyMean,
		TimeDimType_SeasonalMean,
		TimeDimType_AnnualMean
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	NcTimeDimension() :
		m_nctype(ncDouble),
		m_units(""),
		m_dimtype(TimeDimType_Standard)
	{ }

public:
	///	<summary>
	///		Get the type.
	///	</summary>
	NcType nctype() const {
		return m_nctype;
	}

	///	<summary>
	///		Get the units.
	///	</summary>
	const std::string & units() const {
		return m_units;
	}

	///	<summary>
	///		Get the dimension type.
	///	</summary>
	TimeDimType dimtype() const {
		return m_dimtype;
	}

public:
	///	<summary>
	///		Associated NcType.
	///	</summary>
	NcType m_nctype;

	///	<summary>
	///		Associated units.
	///	</summary>
	std::string m_units;

	///	<summary>
	///		Associated type of time.
	///	</summary>
	TimeDimType m_dimtype;
};

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load one of several possible names for a given variable.
///	</summary>
///	<returns>
///		The index in the vecVarNames array of the variable found, or
///		vecVarNames.size() if not found.
///	</returns>
size_t NcGetVarFromList(
	NcFile & ncFile,
	const std::vector<std::string> & vecVarNames,
	NcVar ** pvar = NULL,
	NcDim ** pdim = NULL
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the name of the latitude and longitude dimension from a file.
///	</summary>
void NcGetLatitudeLongitudeName(
	NcFile & ncFile,
	std::string & strLatitudeName,
	std::string & strLongitudeName
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if a given NcDim is a time dimension.
///	</summary>
bool NcIsTimeDimension(
	NcDim * dim
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the time dimension from the NetCDF file.
///	</summary>
NcDim * NcGetTimeDimension(
	NcFile & ncFile
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the time variable from the NetCDF file.
///	</summary>
NcVar * NcGetTimeVariable(
	NcFile & ncFile
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy NetCDF attribute metadata from one file to another.
///	</summary>
void CopyNcFileAttributes(
	NcFile * fileIn,
	NcFile * fileOut
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy NetCDF attribute metadata from one variable to another.
///	</summary>
void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Insert a new dimension into the NcFile, or use existing.
///	</summary>
NcDim * AddNcDimOrUseExisting(
	NcFile & ncFile,
	const std::string & strDimName,
	long lDimSize
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy a NetCDF variable from one file to another.
///	</summary>
void CopyNcVar(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes = true,
	bool fCopyData = true
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy a NetCDF varaible from one file to another if it exists.
///	</summary>
void CopyNcVarIfExists(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes = true,
	bool fCopyData = true
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy a NetCDF variable from one file to another.
///	</summary>
void CopyNcVarTimeSubset(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	const std::vector<Time> & vecOutputTimes
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Read the time data from a NetCDF file.
///	</summary>
void ReadCFTimeDataFromNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	NcTimeDimension & vecTimes,
	bool fWarnOnMissingCalendar
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Write the time data to a NetCDF file.
///	</summary>
void WriteCFTimeDataToNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	NcTimeDimension & vecTimes,
	bool fRecordDim = true
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a value-based index into an integer index.
///	</summary>
long GetIntegerIndexFromValueBasedIndex(
	NcFile * ncfile,
	const std::string & strFilename,
	const std::string & strVariableName,
	long lDim,
	const std::string & strValueIndex
);

////////////////////////////////////////////////////////////////////////////////

#endif

