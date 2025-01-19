///////////////////////////////////////////////////////////////////////////////
///
///	\file    AutoCurator.cpp
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

#include "Defines.h"
#include "AutoCurator.h"
#include "Exception.h"
#include "DataArray1D.h"
#include "Variable.h"
#include "Announce.h"

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////
// AutoCuratorDataset
///////////////////////////////////////////////////////////////////////////////

void AutoCuratorDataset::InsertTimeToFileTimeIx(
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

	iter->second.push_back(FileTimeIx(iFileIx, iTimeIx));
}

///////////////////////////////////////////////////////////////////////////////

AutoCuratorDataset::FilenameTimePairVector AutoCuratorDataset::Find(
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

void AutoCuratorDataset::ToYAMLFile(
	const std::string & strYAMLFile
) {
	std::ofstream ofAC(strYAMLFile);
	if (!ofAC.is_open()) {
		_EXCEPTION1("Unable to open file \"%s\" for reading", strYAMLFile.c_str());
	}

	ofAC << "---" << std::endl;
	ofAC << "calendar: " << Time::StringFromCalendarType(m_eCalendarType) << std::endl;
	ofAC << "time_nctype: " << (int)(m_eNcTimeType) << std::endl;
	ofAC << "time_units: \"" << m_strNcTimeUnits << "\"" << std::endl;
	ofAC << "files:" << std::endl;
	for (size_t f = 0; f < m_vecFiles.size(); f++) {
		ofAC << "  - \"" << m_vecFiles[f] << "\"" << std::endl;
	}
	ofAC << "times:" << std::endl;
	for (auto it = m_mapTimeToTimeFileIx.begin(); it != m_mapTimeToTimeFileIx.end(); it++) {
		ofAC << "  - \"" << it->first.ToString() << "\"" << std::endl;
		for (size_t ft = 0; ft < it->second.size(); ft++) {
			ofAC << "    - [" << it->second[ft].first << "," << it->second[ft].second << "]" << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void AutoCuratorDataset::FromYAMLFile(
	const std::string & strYAMLFile
) {
	size_t sLine = 0;

	m_vecFiles.clear();
	m_mapTimeToTimeFileIx.clear();

	std::ifstream ifAC(strYAMLFile);
	if (!ifAC.is_open()) {
		_EXCEPTION1("Unable to open file \"%s\" for reading", strYAMLFile.c_str());
	}

	std::string strBuf;

	std::getline(ifAC, strBuf); sLine++;
	if (strBuf != "---") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}

	// Calendar
	std::getline(ifAC, strBuf); sLine++;
	if (strBuf.substr(0,9) != "calendar:") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	strBuf = strBuf.substr(10);
	STLStringHelper::RemoveWhitespaceInPlace(strBuf);
	m_eCalendarType = Time::CalendarTypeFromString(strBuf);

	// Time nctype
	std::getline(ifAC, strBuf); sLine++;
	if (strBuf.substr(0,12) != "time_nctype:") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	strBuf = strBuf.substr(12);
	STLStringHelper::RemoveWhitespaceInPlace(strBuf);
	m_eNcTimeType = (NcType)std::stoi(strBuf);

	// Time units
	std::getline(ifAC, strBuf); sLine++;
	if (strBuf.substr(0,11) != "time_units:") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	strBuf = strBuf.substr(11);
	STLStringHelper::RemoveWhitespaceInPlace(strBuf);
	if (strBuf.length() < 2) {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	if ((strBuf[0] != '\"') || (strBuf[strBuf.length()-1] != '\"')) {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	m_strNcTimeUnits = strBuf.substr(1,strBuf.length()-2);

	// File names
	std::getline(ifAC, strBuf); sLine++;
	if (strBuf != "files:") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	std::getline(ifAC, strBuf); sLine++;
	while(strBuf != "times:") {
		if (ifAC.eof()) {
			_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
		}
		if (strBuf.substr(0,4) != "  - ") {
			_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
		}
		strBuf = strBuf.substr(4);
		STLStringHelper::RemoveWhitespaceInPlace(strBuf);
		if ((strBuf[0] != '\"') || (strBuf[strBuf.length()-1] != '\"')) {
			_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
		}
		m_vecFiles.push_back(strBuf.substr(1,strBuf.length()-2));
		std::getline(ifAC, strBuf); sLine++;
	}

	// Times
	if (strBuf != "times:") {
		_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
	}
	std::getline(ifAC, strBuf); sLine++;
	while(!ifAC.eof()) {
		if (strBuf.substr(0,4) != "  - ") {
			_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
		}
		strBuf = strBuf.substr(4);
		STLStringHelper::RemoveWhitespaceInPlace(strBuf);
		if ((strBuf[0] != '\"') || (strBuf[strBuf.length()-1] != '\"')) {
			_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
		}
		Time time(m_eCalendarType);
		time.FromFormattedString(strBuf.substr(1,strBuf.length()-2));
		std::pair<TimeToFileTimeIxMap::iterator, bool> map_insert_res =
			m_mapTimeToTimeFileIx.insert(
				std::pair<Time, FileTimeIxVector>(
					time, FileTimeIxVector()));
		if (!map_insert_res.second) {
			_EXCEPTIONT("Error inserting object into map");
		}

		// Load file time pairs
		std::getline(ifAC, strBuf); sLine++;
		while(!ifAC.eof()) {
			if (strBuf.substr(0,4) == "  - ") {
				break;
			}
			if (strBuf.substr(0,6) != "    - ") {
				_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
			}
			strBuf = strBuf.substr(6);
			if (strBuf.length() < 2) {
				_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
			}
			if ((strBuf[0] != '[') || (strBuf[strBuf.length()-1] != ']')) {
				_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
			}
			std::string strBufFirst;
			std::string strBufSecond;
			for (int ix = 1; ix < strBuf.length()-1; ix++) {
				if (strBuf[ix] == ',') {
					strBufFirst = strBuf.substr(1,ix-1);
					strBufSecond = strBuf.substr(ix+1,strBuf.length()-ix-2);
					break;
				}
			}
			if (!STLStringHelper::IsInteger(strBufFirst)) {
				_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
			}
			if (!STLStringHelper::IsInteger(strBufSecond)) {
				_EXCEPTION2("Invalid format of \"%s\" on line %lu", strYAMLFile.c_str(), sLine);
			}

			map_insert_res.first->second.push_back(
				FileTimeIx(
					std::stoi(strBufFirst),
					std::stoi(strBufSecond)));

			std::getline(ifAC, strBuf); sLine++;
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
// AutoCurator
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

	// Open file
	NcFile ncFile(strFile.c_str());
	if (!ncFile.is_valid()) {
		_EXCEPTION1("Unable to open input file \"%s\"", strFile.c_str());
	}

	// Check for time dimension and variable
	NcDim * dimTime = NcGetTimeDimension(ncFile);
	NcVar * varTime = NcGetTimeVariable(ncFile);

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

		// Add file to list of files indexed
		int iFileIx = static_cast<int>(m_vecFiles.size());
		m_vecFiles.push_back(strFile);

		// Index file
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
			_EXCEPTION1("\"time\" variable may have only one dimension in file \"%s\"",
				strFile.c_str());
		}
		if (!NcIsTimeDimension(varTime->get_dim(0)) != 0) {
			_EXCEPTION1("First dimension of \"time\" is not a time dimension in file \"%s\"",
				strFile.c_str());
		}

		// Check if this file is a daily mean climatology
		bool fDailyMeanClimatology = false;
		NcAtt * attType = varTime->get_att("type");
		if (attType != NULL) {
			std::string strType = attType->as_string(0);
			if (strType == "daily mean climatology") {
				fDailyMeanClimatology = true;
			}
		}

		AutoCuratorDataset * pacd = this;
		if (fDailyMeanClimatology) {
			pacd = &m_acdDailyMean;
		}

		// Get calendar and units information
		Time::CalendarType eCalendarType = Time::CalendarUnknown;

		NcAtt * attCalendar = varTime->get_att("calendar");
		NcAtt * attUnits = varTime->get_att("units");
		if (attCalendar == NULL) {
			if (pacd->m_eCalendarType != Time::CalendarUnknown) {
				eCalendarType = pacd->m_eCalendarType;
				Announce("WARNING: \"time\" variable in \"%s\" missing "
					"\"calendar\" attribute. Assuming calendar is \"%s\"",
					strFile.c_str(),
					Time::StringFromCalendarType(eCalendarType).c_str());
			} else {
				_EXCEPTION1("\"time\" variable in \"%s\" missing "
					"\"calendar\" attribute", strFile.c_str());
			}
		} else {
			eCalendarType = Time::CalendarTypeFromString(attCalendar->as_string(0));
		}

		if (attUnits == NULL) {
			_EXCEPTION1("\"time\" variable in \"%s\" missing "
				"\"units\" attribute", strFile.c_str());
		}
		std::string strCurrentNcTimeUnits = attUnits->as_string(0);

		// Index file
		int iFileIx = static_cast<int>(pacd->m_vecFiles.size());
		pacd->m_vecFiles.push_back(strFile);

		if (pacd->m_strNcTimeUnits == "") {
			pacd->m_eNcTimeType = varTime->type();
			pacd->m_strNcTimeUnits = strCurrentNcTimeUnits;
		}

		if (pacd->m_eCalendarType == Time::CalendarUnknown) {
			pacd->m_eCalendarType = eCalendarType;
		} else if (pacd->m_eCalendarType != eCalendarType) {
			_EXCEPTION2("CalendarType mismatch in \"%s\", found \"%s\"",
				strFile.c_str(),
				attCalendar->as_string(0));
		}

		DataArray1D<int> vecTimeInt;
		DataArray1D<float> vecTimeFloat;
		DataArray1D<double> vecTimeDouble;
		if (varTime->type() == ncInt) {
			vecTimeInt.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeInt[0]), dimTime->size());

		} else if (varTime->type() == ncInt64) {
			DataArray1D<ncint64> vecTimeInt64;
			vecTimeInt64.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeInt64[0]), dimTime->size());

			vecTimeInt.Allocate(dimTime->size());
			for (int t = 0; t < dimTime->size(); t++) {
				vecTimeInt[t] = static_cast<int>(vecTimeInt64[t]);
			}

		} else if (varTime->type() == ncFloat) {
			vecTimeFloat.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeFloat[0]), dimTime->size());

		} else if (varTime->type() == ncDouble) {
			vecTimeDouble.Allocate(dimTime->size());
			varTime->set_cur((long)0);
			varTime->get(&(vecTimeDouble[0]), dimTime->size());

		} else {
			_EXCEPTION1("Variable \"time\" has invalid type in file \"%s\" (require int, float or double)",
				strFile.c_str());
		}

		for (int t = 0; t < dimTime->size(); t++) {
			Time time(pacd->m_eCalendarType);
			if ((varTime->type() == ncInt) || (varTime->type() == ncInt64)) {
				time.FromCFCompliantUnitsOffsetInt(
					strCurrentNcTimeUnits,
					vecTimeInt[t]);

			} else if (varTime->type() == ncFloat) {
				time.FromCFCompliantUnitsOffsetDouble(
					strCurrentNcTimeUnits,
					static_cast<double>(vecTimeFloat[t]));

			} else if (varTime->type() == ncDouble) {
				time.FromCFCompliantUnitsOffsetDouble(
					strCurrentNcTimeUnits,
					vecTimeDouble[t]);
			}

#if defined(ROUND_TIMES_TO_NEAREST_MINUTE)
			time.RoundToNearestMinute();
#endif

			pacd->InsertTimeToFileTimeIx(time, iFileIx, t);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

bool AutoCurator::FindFilesAtTime(
	const Time & time,
	NcFileVector & vecncDataFiles
) const {

	// Get normal files at this time
	{
		FilenameTimePairVector vec = Find(time);

		vecncDataFiles.clear();
		vecncDataFiles.SetTime(time);
		for (int i = 0; i < vec.size(); i++) {
			vecncDataFiles.InsertFile(
				vec[i].first.c_str(),
				vec[i].second);
		}
	}

	// Get daily mean climatology at this time
	{
		Time timeDailyMean(m_acdDailyMean.m_eCalendarType);
		timeDailyMean.SetYear(1);
		timeDailyMean.SetMonth(time.GetMonth());
		timeDailyMean.SetDay(time.GetDay());
		if ((time.IsLeapDay()) &&
		    (m_acdDailyMean.m_eCalendarType == Time::CalendarNoLeap)
		) {
			timeDailyMean.SetDay(28);
		}

		FilenameTimePairVector vec = m_acdDailyMean.Find(timeDailyMean);

		for (int i = 0; i < vec.size(); i++) {
			Announce("Using daily mean climatology from %s (index %i)",
				time.ToDateString().c_str(),
				vec[i].second);

			vecncDataFiles.InsertFile(
				vec[i].first.c_str(),
				vec[i].second);
		}
	}

	if (vecncDataFiles.size() == 0) {
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

bool AutoCurator::FindFilesNearTime(
	const Time & time,
	NcFileVector & vecncDataFiles,
	const Time & timeMaxDelta,
	bool fVerbose
) const {
	if (m_mapTimeToTimeFileIx.size() == 0) {
		return false;
	}

	// Find first element before time and store in iter0
	// Find first element after time and store in iter1
	TimeToFileTimeIxMap::const_iterator iter1 =
		m_mapTimeToTimeFileIx.lower_bound(time);

	TimeToFileTimeIxMap::const_iterator iter0 = iter1;

	if (iter0 != m_mapTimeToTimeFileIx.begin()) {
		iter0--;
	}
	if (iter1 == m_mapTimeToTimeFileIx.end()) {
		iter1--;
	}

	// TODO: Use exact arithmetic
	double dMaxDeltaSeconds = timeMaxDelta.AsSeconds();
	_ASSERT(dMaxDeltaSeconds >= 0.0);

	// Find the closest time within timeMaxDelta and use that instead
	double dTimeDelta0 = fabs(time.DeltaSeconds(iter0->first));
	double dTimeDelta1 = fabs(time.DeltaSeconds(iter1->first));

	if (dTimeDelta0 < dTimeDelta1) {
		if (dTimeDelta0 <= dMaxDeltaSeconds) {
			if (dTimeDelta0 != 0.0) {
				Announce("Substituting time %s", iter0->first.ToString().c_str());
			}
			return FindFilesAtTime(iter0->first, vecncDataFiles);
		}
	} else {
		if (dTimeDelta1 <= dMaxDeltaSeconds) {
			if (dTimeDelta1 != 0.0) {
				Announce("Substituting time %s", iter1->first.ToString().c_str());
			}
			return FindFilesAtTime(iter1->first, vecncDataFiles);
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

