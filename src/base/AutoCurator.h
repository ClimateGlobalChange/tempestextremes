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

class AutoCurator {

public:
	///	<summary>
	///		A local (file index, local time index) pair.
	///	</summary>
	typedef std::pair<int, int> LocalFileTimePair;

	///	<summary>
	///		A (filename, local time index) pair.
	///	</summary>
	typedef std::pair<std::string, int> FilenameTimePair;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	AutoCurator() :
		m_eCalendarType(Time::CalendarUnknown)
	{ }

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
	///		Find the file/time pair associated with a given Time.
	///	</summary>
	FilenameTimePair Find(
		const Time & time
	);

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
	std::map<Time, LocalFileTimePair> m_mapTimeToTimeFileIx;
};

///////////////////////////////////////////////////////////////////////////////


#endif

