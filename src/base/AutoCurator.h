///////////////////////////////////////////////////////////////////////////////
///
///	\file    AutoCurator.h
///	\author  Paul Ullrich
///	\version August 15, 2018
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

#ifndef _AUTOCURATOR_H_
#define _AUTOCURATOR_H_

#include <vector>
#include <map>
#include <string>
#include "netcdfcpp.h"

#include "TimeObj.h"

///////////////////////////////////////////////////////////////////////////////

class NcFileVector;

///////////////////////////////////////////////////////////////////////////////

class AutoCuratorDataset {

public:
	///	<summary>
	///		A local (file index, local time index) pair.
	///	</summary>
	typedef std::pair<int, int> FileTimeIx;

	///	<summary>
	///		A vector of (file index, local time index) pairs.
	///	</summary>
	typedef std::vector<FileTimeIx> FileTimeIxVector;

	///	<summary>
	///		A map from a Time to a vector of (file index, local time index) pairs.
	///	</summary>
	typedef std::map<Time, FileTimeIxVector> TimeToFileTimeIxMap;

	///	<summary>
	///		A (filename, local time index) pair.
	///	</summary>
	typedef std::pair<std::string, int> FilenameTimePair;

	///	<summary>
	///		A vector of (filename, local time index) pairs.
	///	</summary>
	typedef std::vector<FilenameTimePair> FilenameTimePairVector;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	AutoCuratorDataset() :
		m_eCalendarType(Time::CalendarUnknown),
		m_eNcTimeType(ncNoType),
		m_strNcTimeUnits("")
	{ }

public:
	///	<summary>
	///		Create a new entry in m_mapTimeToTimeFileIx for the given Time.
	///	</summary>
	void InsertTimeToFileTimeIx(
		const Time & time,
		const int & iFileIx,
		const int & iTimeIx
	);

	///	<summary>
	///		Find the file/time pair associated with a given Time.
	///	</summary>
	FilenameTimePairVector Find(
		const Time & time
	) const;

public:
	///	<summary>
	///		Write the AutoCuratorDataset to a file.
	///	</summary>
	void ToYAMLFile(
		const std::string & strFile
	);

	///	<summary>
	///		Read from a file.
	///	</summary>
	void FromYAMLFile(
		const std::string & strFile
	);

public:
	///	<summary>
	///		CalendarType for all Time objects.
	///	</summary>
	Time::CalendarType m_eCalendarType;

	///	<summary>
	///		NcType used for time in the NetCDF files.
	///	</summary>
	NcType m_eNcTimeType;

	///	<summary>
	///		CF-compliant time units.
	///	</summary>
	std::string m_strNcTimeUnits;

	///	<summary>
	///		List of filenames.
	///	</summary>
	std::vector<std::string> m_vecFiles;

	///	<summary>
	///		A map between Time and (file,time) index.
	///	</summary>
	TimeToFileTimeIxMap m_mapTimeToTimeFileIx;
};

///////////////////////////////////////////////////////////////////////////////

class AutoCurator : public AutoCuratorDataset {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	AutoCurator()
	{ }

public:
	///	<summary>
	///		Change the calendar type.
	///	</summary>
	void SetCalendar(
		const Time::CalendarType eCalendarType
	) {
		m_eCalendarType = eCalendarType;
	}

public:
	///	<summary>
	///		Index the contents of a single file or list of files delimited
	///		by semicolons.
	///	</summary>
	void IndexFiles(
		const std::string & strFile
	);

	///	<summary>
	///		Get the CalendarType.
	///	</summary>
	Time::CalendarType GetCalendarType() const {
		return m_eCalendarType;
	}

	///	<summary>
	///		Get the vector of filenames.
	///	</summary>
	const std::vector<std::string> & GetFilenames() const {
		return m_vecFiles;
	}

	///	<summary>
	///		Get the total count of time slices.
	///	</summary>
	size_t GetTimeCount() const {
		return m_mapTimeToTimeFileIx.size();
	}

	///	<summary>
	///		Get the type used to store time.
	///	</summary>
	NcType GetNcTimeType() const {
		return m_eNcTimeType;
	}

	///	<summary>
	///		Get the units used to store time.
	///	</summary>
	const std::string & GetNcTimeUnits() const {
		return m_strNcTimeUnits;
	}

public:
	///	<summary>
	///		Generate a NcFileVector for the given Time.
	///	</summary>
	bool FindFilesAtTime(
		const Time & time,
		NcFileVector & vecncDataFiles
	) const;

	///	<summary>
	///		Generate a NcFileVector for a Time near the given Time.
	///	</summary>
	bool FindFilesNearTime(
		const Time & time,
		NcFileVector & vecncDataFiles,
		const Time & timeMaxDelta,
		bool fVerbose = true
	) const;

protected:
	///	<summary>
	///		Daily mean dataset.
	///	</summary>
	AutoCuratorDataset m_acdDailyMean;
};

///////////////////////////////////////////////////////////////////////////////


#endif

