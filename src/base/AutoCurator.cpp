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
#include "DataArray1D.h"
#include "Variable.h"

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

void AutoCurator::InsertTimeToFileTimeIx(
	const Time & time,
	const int & iFileIx,
	const int & iTimeIx
) {
	TimeToFileTimeIxMap::iterator iter =
		m_mapTimeToTimeFileIx.find(time);

	if (iter == m_mapTimeToTimeFileIx.end()) {
		iter = m_mapTimeToTimeFileIx.insert(
			TimeToFileTimeIxMap::value_type(
				time, FileTimeIxVector())).first;
	}

	//std::cout << time.ToString() << " " << iFileIx << " " << iTimeIx << std::endl;

	iter->second.push_back(FileTimeIx(iFileIx, iTimeIx));
}

///////////////////////////////////////////////////////////////////////////////

void AutoCurator::IndexFiles(
	const std::string & strFile
) {
	// Search for semi-colons in strFile and break up accordingly
	int iLast = 0;
	for (int i = 0; i < strFile.length(); i++) {
		if (strFile[i] == ';') {
			IndexFiles(strFile.substr(iLast, i-iLast));
			iLast = i+1;
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

			InsertTimeToFileTimeIx(time, iFileIx, t);
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

		DataArray1D<int> vecTimeInt;
		DataArray1D<float> vecTimeFloat;
		DataArray1D<double> vecTimeDouble;
		if (varTime->type() == ncInt) {
			vecTimeInt.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeInt[0]), dimTime->size());

		} else if (varTime->type() == ncFloat) {
			vecTimeFloat.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeFloat[0]), dimTime->size());

		} else if (varTime->type() == ncDouble) {
			vecTimeDouble.Allocate(dimTime->size());
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

			InsertTimeToFileTimeIx(time, iFileIx, t);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

AutoCurator::FilenameTimePairVector AutoCurator::Find(
	const Time & time
) const {
	FilenameTimePairVector vec;

	TimeToFileTimeIxMap::const_iterator iter =
		m_mapTimeToTimeFileIx.find(time);

	if (iter == m_mapTimeToTimeFileIx.end()) {
		return vec;
	}

	const FileTimeIxVector & vecFileTimeIx = iter->second;

	for (int f = 0; f < vecFileTimeIx.size(); f++) {
		vec.push_back(
			FilenameTimePair(
				m_vecFiles[vecFileTimeIx[f].first],
				vecFileTimeIx[f].second));
	}

	return vec;
}

///////////////////////////////////////////////////////////////////////////////

bool AutoCurator::FindFilesAtTime(
	const Time & time,
	NcFileVector & vecncDataFiles,
	int & iTime
) const {
	FilenameTimePairVector vec = Find(time);

	vecncDataFiles.clear();
	iTime = (-1);

	for (int i = 0; i < vec.size(); i++) {
		NcFile * pncfile = new NcFile(vec[i].first.c_str());
		if (pncfile == NULL) {
			_EXCEPTIONT("Unable to allocate new NcFile");
		}
		if (!pncfile->is_valid()) {
			_EXCEPTION1("Unable to open data file \"%s\"", vec[i].first.c_str());
		}

		vecncDataFiles.push_back(pncfile);

		if (iTime == (-1)) {
			iTime = vec[i].second;
		} else if (iTime != vec[i].second) {
			_EXCEPTIONT("Data files have different local time indices (unsupported)");
		}
	}

	if (vecncDataFiles.size() == 0) {
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

