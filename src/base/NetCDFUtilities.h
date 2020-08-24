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
	///		Get the type.
	///	</summary>
	NcType type() const {
		return m_nctype;
	}

	///	<summary>
	///		Get the units.
	///	</summary>
	const std::string & units() const {
		return m_units;
	}

public:
	///	<summary>
	///		Associated type.
	///	</summary>
	NcType m_nctype;

	///	<summary>
	///		Associated units.
	///	</summary>
	std::string m_units;
};

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

