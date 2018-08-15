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

#include "AutoCurator.h"
#include "Exception.h"
#include "DataVector.h"

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

void AutoCurator::IndexFiles(
	const std::string & strFile
) {
	// Search for semi-colons in strFile and break up accordingly
	int iLast = 0;
	for (int i = 0; i < strFile.length(); i++) {
		if (strFile[i] == ';') {
			iLast = i+1;
			IndexFiles(strFile.substr(iLast, i-iLast));
		}
	}
	if (iLast != 0) {
		IndexFiles(strFile.substr(iLast));
		return;
	}

	// Make sure file doesn't exist already
	for (int f = 0; f < m_vecFiles.size(); f++) {
		if (m_vecFiles[f] == strFile) {
			_EXCEPTION1("Duplicated filename \"%s\"", strFile.c_str());
		}
	}

	// Add file to list of files indexed
	int iFileIx = static_cast<int>(m_vecFiles.size());
	m_vecFiles.push_back(strFile);

	// Open file
	NcFile ncFile(strFile.c_str());
	if (!ncFile.is_valid()) {
		_EXCEPTION1("Unable to open input file \"%s\"", strFile.c_str());
	}

	// Check for time dimension and variable
	NcDim * dimTime = ncFile.get_dim("time");
	NcVar * varTime = ncFile.get_var("time");

	if (dimTime == NULL) {
		if (varTime != NULL) {
			_EXCEPTION1("File \"%s\" has \"time\" variable but no \"time\" dimension",
				strFile.c_str());
		}
		return;
	}

	// Check for time variable; if not present use CalendarNone and
	// fill in FileTimePairs.
	if (varTime == NULL) {
		if (m_eCalendarType == Time::CalendarUnknown) {
			m_eCalendarType = Time::CalendarNone;

		} else if (m_eCalendarType != Time::CalendarNone) {
			_EXCEPTION1("\"time\" variable in \"%s\" does not exist "
				"although \"time\" dimension is present", strFile.c_str());
		}

		for (int t = 0; t < dimTime->size(); t++) {
			Time time(m_eCalendarType);
			time.SetYear(iFileIx);
			time.SetMonth(t);

			m_mapTimeToTimeFileIx.insert(
				std::map<Time, LocalFileTimePair>::value_type(
					time, LocalFileTimePair(iFileIx, t)));
		}

	// If time variable is present make sure it is CF-compliant and
	// fill in FileTimePairs.
	} else {
		if (varTime->num_dims() != 1) {
			_EXCEPTION1("\"time\" variable requires one dimension "
				"with name \"time\" in file \"%s\"",
				strFile.c_str());
		}
		if (strcmp(varTime->get_dim(0)->name(), "time") != 0) {
			_EXCEPTION1("\"time\" variable requires one dimension "
				"with name \"time\" in file \"%s\"",
				strFile.c_str());
		}

		NcAtt * attCalendar = varTime->get_att("calendar");
		NcAtt * attUnits = varTime->get_att("units");
		if (attCalendar == NULL) {
			_EXCEPTION1("\"time\" variable in \"%s\" missing "
				"\"calendar\" attribute", strFile.c_str());
		}
		if (attUnits == NULL) {
			_EXCEPTION1("\"time\" variable in \"%s\" missing "
				"\"units\" attribute", strFile.c_str());
		}

		Time::CalendarType eCalendarType =
			Time::CalendarTypeFromString(attCalendar->as_string(0));

		if (m_eCalendarType == Time::CalendarUnknown) {
			m_eCalendarType = eCalendarType;
		} else if (m_eCalendarType != eCalendarType) {
			_EXCEPTION2("CalendarType mismatch in \"%s\", found \"%s\"",
				strFile.c_str(), attCalendar->as_string(0));
		}

		DataVector<int> vecTimeInt;
		DataVector<float> vecTimeFloat;
		DataVector<double> vecTimeDouble;
		if (varTime->type() == ncInt) {
			vecTimeInt.Initialize(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeInt[0]), dimTime->size());

		} else if (varTime->type() == ncFloat) {
			vecTimeFloat.Initialize(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeFloat[0]), dimTime->size());

		} else if (varTime->type() == ncDouble) {
			vecTimeDouble.Initialize(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeDouble[0]), dimTime->size());

		} else {
			_EXCEPTION1("Variable \"time\" has invalid type in file \"%s\"",
				strFile.c_str());
		}

		for (int t = 0; t < dimTime->size(); t++) {
			Time time(m_eCalendarType);
			if (varTime->type() == ncInt) {
				time.FromCFCompliantUnitsOffsetInt(
					attUnits->as_string(0),
					vecTimeInt[t]);

			} else if (varTime->type() == ncFloat) {
				time.FromCFCompliantUnitsOffsetDouble(
					attUnits->as_string(0),
					static_cast<double>(vecTimeFloat[t]));

			} else if (varTime->type() == ncDouble) {
				time.FromCFCompliantUnitsOffsetDouble(
					attUnits->as_string(0),
					vecTimeDouble[t]);
			}

			m_mapTimeToTimeFileIx.insert(
				std::map<Time, LocalFileTimePair>::value_type(
					time, LocalFileTimePair(iFileIx, t)));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

AutoCurator::FilenameTimePair AutoCurator::Find(
	const Time & time
) {
	std::map<Time, LocalFileTimePair>::const_iterator iter =
		m_mapTimeToTimeFileIx.find(time);

	if (iter == m_mapTimeToTimeFileIx.end()) {
		return FilenameTimePair("",-1);	
	} else {
		return
			FilenameTimePair(
				m_vecFiles[iter->second.first],
				iter->second.second);
	}
}

///////////////////////////////////////////////////////////////////////////////

