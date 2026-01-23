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
///		An object holding sizes for a given Variable.
///	</summary>
typedef VariableDimIndex VariableAuxSize;

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

	///	<summary>
	///		Time bounds.
	///	</summary>
	std::vector< std::pair<Time,Time> > m_vecTimeBounds;
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
///		Copy NetCDF attribute metadata from one variable to another.
///	</summary>
void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut,
	const std::vector<std::string> & vecDoNotCopyNames
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
///		Copy a NetCDF variable from one file to another if it exists.
///	</summary>
bool CopyNcVarIfExists(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes = true,
	bool fCopyData = true
);

////////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Copy latitude and longitude variables over and check for consistency.
///	</summary>
void CopyNcLatitudeLongitude(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName,
	const std::vector<size_t> & nGridDim,
	NcDim ** pdimGrid0 = NULL,
	NcDim ** pdimGrid1 = NULL
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
///		Read Time data from a WRF output NetCDF file.
///	</summary>
void ReadWRFTimeDataFromNcFile(
	NcFile * ncfile,
	const std::string & strFilename,
	NcTimeDimension & vecTimes
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
	bool fRecordDim = true,
	bool fAppend = false
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

