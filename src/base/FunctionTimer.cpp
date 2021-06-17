///////////////////////////////////////////////////////////////////////////////
///
///	\file    FunctionTimer.cpp
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "FunctionTimer.h"
#include "Exception.h"

#include <iostream>
#include <sys/time.h>

///////////////////////////////////////////////////////////////////////////////

FunctionTimer::GroupDataMap FunctionTimer::m_mapGroupData;

///////////////////////////////////////////////////////////////////////////////

FunctionTimer::FunctionTimer(const char *szGroup) {

	// Assign group name
	if (szGroup == NULL) {
		m_strGroup = "";
	} else {
		m_strGroup = szGroup;
	}

	// Assign start time
	gettimeofday(&m_tvStartTime, NULL);
}

///////////////////////////////////////////////////////////////////////////////

void FunctionTimer::Reset() {
	gettimeofday(&m_tvStartTime, NULL);
}

///////////////////////////////////////////////////////////////////////////////

unsigned long FunctionTimer::Time(bool fDone) {
	timeval tv;
	gettimeofday(&tv, NULL);

	unsigned long iTime =
	    MICROSECONDS_PER_SECOND * (tv.tv_sec - m_tvStartTime.tv_sec)
	    + (tv.tv_usec - m_tvStartTime.tv_usec);

	// If no name associated with this timer, ignore fDone.
	if (m_strGroup == "") {
		return iTime;
	}

	// Add the time to the group record
	if (fDone) {
		GroupDataMap::iterator iter;

		iter = m_mapGroupData.find(m_strGroup);

		// Add to existing group record
		if (iter != m_mapGroupData.end()) {
			iter->second.iTotalTime += iTime;
			iter->second.nEntries++;

		// Create new group record
		} else {
			GroupDataPair gdp;
			gdp.first = m_strGroup;
			gdp.second.iTotalTime = iTime;
			gdp.second.nEntries = 1;

			m_mapGroupData.insert(gdp);
		}
	}

	return iTime;
}

///////////////////////////////////////////////////////////////////////////////

unsigned long FunctionTimer::StopTime() {
	return Time(true);
}

///////////////////////////////////////////////////////////////////////////////

const FunctionTimer::TimerGroupData & FunctionTimer::GetGroupTimeRecord(
	const char *szName
) {
	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		return iter->second;

	// Group record does not exist
	} else {
		_EXCEPTION1("Group time record %s does not exist.", szName);
	}
}

///////////////////////////////////////////////////////////////////////////////

unsigned long FunctionTimer::GetAverageGroupTime(const char *szName) {

	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		const TimerGroupData & tgd = iter->second;

		return (tgd.iTotalTime / tgd.nEntries);

	// Group record does not exist
	} else {
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////////

unsigned int FunctionTimer::GetNumberOfEntries(const char *szName) {

	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		const TimerGroupData & tgd = iter->second;

		return (tgd.nEntries);

	// Group record does not exist
	} else {
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void FunctionTimer::ResetGroupTimeRecord(const char *szName) {
	GroupDataMap::iterator iter;

	iter = m_mapGroupData.find(szName);

	// Retrieve existing group record
	if (iter != m_mapGroupData.end()) {
		m_mapGroupData.erase(iter);

	// Group record does not exist
	} else {
		_EXCEPTION1("Group time record %s does not exist.", szName);
	}
}

///////////////////////////////////////////////////////////////////////////////

