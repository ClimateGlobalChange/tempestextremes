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

#include "TimeObj.h"

///////////////////////////////////////////////////////////////////////////////

class NcFileVector;

///////////////////////////////////////////////////////////////////////////////

class AutoCurator {

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
	AutoCurator() :
		m_eCalendarType(Time::CalendarUnknown)
	{ }

protected:
	///	<summary>
	///		Create a new entry in m_mapTimeToTimeFileIx for the given Time.
	///	</summary>
	void InsertTimeToFileTimeIx(
		const Time & time,
		const int & iFileIx,
		const int & iTimeIx
	);

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
	///		Find the file/time pair associated with a given Time.
	///	</summary>
	FilenameTimePairVector Find(
		const Time & time
	) const;

	///	<summary>
	///		Get the total count of time slices.
	///	</summary>
	size_t GetTimeCount() const {
		return m_mapTimeToTimeFileIx.size();
	}

public:
	///	<summary>
	///		Generate a NcFileVector and local time index for the given Time.
	///	</summary>
	bool FindFilesAtTime(
		const Time & time,
		NcFileVector & vecncDataFiles,
		int & iTime
	) const;

protected:
	///	<summary>
	///		CalendarType for all Time objects.
	///	</summary>
	Time::CalendarType m_eCalendarType;

	///	<summary>
	///		List of filenames.
	///	</summary>
	std::vector<std::string> m_vecFiles;

	///	<summary>
	///		A map between Time, and (file,time) index.
	///	</summary>
	TimeToFileTimeIxMap m_mapTimeToTimeFileIx;
};

///////////////////////////////////////////////////////////////////////////////


#endif

