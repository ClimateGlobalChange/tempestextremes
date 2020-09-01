///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchNodes.cpp
///	\author  Paul Ullrich
///	\version August 14, 2018
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "FilenameList.h"
#include "NodeFileUtilities.h"
#include "CoordTransforms.h"

#include "kdtree.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////

void ParseVariableList(
	const std::string & strVariables,
	std::vector< std::string > & vecVariableStrings
) {
	int iVarBegin = 0;
	int iVarCurrent = 0;

	// Parse variable name
	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
			(strVariables[iVarCurrent] == ',') ||
			(strVariables[iVarCurrent] == ' ') ||
			(strVariables[iVarCurrent] == '\t') ||
			(strVariables[iVarCurrent] == '\n') ||
			(strVariables[iVarCurrent] == '\r')
		) {
			if (iVarCurrent == iVarBegin) {
				if (iVarCurrent >= strVariables.length()) {
					break;
				}

				iVarCurrent++;
				iVarBegin++;
				continue;
			}

			vecVariableStrings.push_back(
				strVariables.substr(iVarBegin, iVarCurrent - iVarBegin));

			iVarBegin = iVarCurrent + 1;
		}

		iVarCurrent++;
	}
}

///////////////////////////////////////////////////////////////////////////////

typedef std::vector< std::vector<std::string> > TimesliceCandidateInformation;

typedef std::pair<Time, TimesliceCandidateInformation> TimeToCandidateInfoPair;

typedef std::map<Time, TimesliceCandidateInformation> TimeToCandidateInfoMap;

typedef TimeToCandidateInfoMap::iterator TimeToCandidateInfoMapIterator;

///////////////////////////////////////////////////////////////////////////////

void ParseDetectNodesFile(
	const std::string & strInputFile,
	const std::vector< std::string > & vecFormatStrings,
	TimeToCandidateInfoMap & mapCandidateInfo,
	Time::CalendarType caltype,
	bool fAllowRepeatedTimes
) {
	// Open file for reading
	FILE * fp = fopen(strInputFile.c_str(), "r");

	if (fp == NULL) {
		_EXCEPTION1("Unable to open input file \"%s\"", strInputFile.c_str());
	}

	// Buffer storage
	std::string strLine;
	char szBuffer[1024];

	// Insufficient candidate information warning
	bool fWarnInsufficientCandidateInfo = false;

	// Current read state
	enum ReadState {
		ReadState_Time,
		ReadState_Candidate
	} eReadState = ReadState_Time;

	// Number of entries per candidate
	int nFormatEntries = vecFormatStrings.size();

	// Current candidate at this time
	int iCandidate = 0;

	// Current candidate information
	TimeToCandidateInfoMapIterator iterCurrentTime = mapCandidateInfo.end();

	// Total candidates at this time
	int nCandidates = 0;

	for (;;) {

		// Load in one line
		strLine.clear();
		for (;;) {
			fgets(szBuffer, 1024, fp);
			strLine += szBuffer;
			if (strlen(szBuffer) != 1023) {
				break;
			}
		}

		// Check for eof
		if (feof(fp)) {
			break;
		}

		// Check for blank line
		if (strLine.size() == 0) {
			continue;
		}

		// Check for comment
		if (strLine[0] == '#') {
			continue;
		}

		// Ignore comments
		if (strlen(szBuffer) > 0) {
			if (szBuffer[0] == '#') {
				continue;
			}
		}

		// Parse the time
		if (eReadState == ReadState_Time) {

			// Parse the time data
			std::vector<std::string> vecTimeString;
			ParseVariableList(strLine, vecTimeString);

			if (vecTimeString.size() != 5) {
				_EXCEPTION1("Malformed time string:\n%s", strLine.c_str());
			}

			iCandidate = 0;
			nCandidates = std::stoi(vecTimeString[3]);

			Time time(caltype);
			time.SetYear(std::stoi(vecTimeString[0]));
			time.SetMonth(std::stoi(vecTimeString[1]));
			time.SetDay(std::stoi(vecTimeString[2]));
			time.SetSecond(std::stoi(vecTimeString[4]) * 3600);

			auto it = mapCandidateInfo.find(time);
			if (it != mapCandidateInfo.end()) {
				if (fAllowRepeatedTimes) {
					iterCurrentTime = it;
					iCandidate = iterCurrentTime->second.size();
					iterCurrentTime->second.resize(iCandidate + nCandidates);

				} else {
					_EXCEPTION1("Repated time \"%s\" found in candidate files",
						time.ToString().c_str());
				}

			} else {
				auto ins =
					mapCandidateInfo.insert(
						TimeToCandidateInfoPair(
							time, TimesliceCandidateInformation(nCandidates)));

				if (!ins.second) {
					_EXCEPTION1("Insertion of time \"%s\" into candidate info map failed",
						time.ToString().c_str());
				}

				iterCurrentTime = ins.first;
			}

			// Prepare to parse candidate data
			if (nCandidates != 0) {
				eReadState = ReadState_Candidate;
			}

		// Parse candidate information
		} else if (eReadState == ReadState_Candidate) {

			// Parse candidates
			_ASSERT(iterCurrentTime != mapCandidateInfo.end());

			TimesliceCandidateInformation & tscinfo = iterCurrentTime->second;

			_ASSERT(iCandidate < tscinfo.size());

			ParseVariableList(strLine, tscinfo[iCandidate]);

			if (tscinfo[iCandidate].size() != nFormatEntries) {
				fWarnInsufficientCandidateInfo = true;
			}

			iCandidate++;
			if (iCandidate == tscinfo.size()) {
				eReadState = ReadState_Time;
				iCandidate = 0;
				iterCurrentTime = mapCandidateInfo.end();
			}
		}
	}

	fclose(fp);	

	// Insufficient candidate information
	if (fWarnInsufficientCandidateInfo) {
		Announce("WARNING: One or more candidates do not have matching"
				" --in_fmt entries");
	}
}

///////////////////////////////////////////////////////////////////////////////

struct Node {
	double x;
	double y;
	double z;

	double lonRad;
	double latRad;
};

///////////////////////////////////////////////////////////////////////////////

class SimplePathSegment {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SimplePathSegment(
		int iTime0,
		int iCandidate0,
		int iTime1,
		int iCandidate1
	) {
		m_iTime[0] = iTime0;
		m_iTime[1] = iTime1;

		m_iCandidate[0] = iCandidate0;
		m_iCandidate[1] = iCandidate1;
	}

public:
	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator< (const SimplePathSegment & seg) const {
		std::pair<int,int> pr0(m_iTime[0], m_iCandidate[0]);
		std::pair<int,int> pr1(seg.m_iTime[0], seg.m_iCandidate[0]);
		return (pr0 < pr1);
	}

public:
	///	<summary>
	///		Begin and end time.
	///	</summary>
	int m_iTime[2];

	///	<summary>
	///		Begin and end candidate.
	///	</summary>
	int m_iCandidate[2];
};

typedef std::set<SimplePathSegment> SimplePathSegmentSet;

///////////////////////////////////////////////////////////////////////////////

class SimplePath {

public:
	///	<summary>
	///		Array of times.
	///	</summary>
	std::vector<int> m_iTimes;

	///	<summary>
	///		Array of candidates.
	///	</summary>
	std::vector<int> m_iCandidates;
};

typedef std::vector<SimplePath> SimplePathVector;

///////////////////////////////////////////////////////////////////////////////

class PathThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		GreaterThan,
		LessThan,
		GreaterThanEqualTo,
		LessThanEqualTo,
		EqualTo,
		NotEqualTo,
		AbsGreaterThanEqualTo,
		AbsLessThanEqualTo
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		const std::string & strOp,
		const std::vector< std::string > & vecFormatStrings
	) {
		// Read mode
		enum {
			ReadMode_Column,
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_MinCount,
			ReadMode_Invalid
		} eReadMode = ReadMode_Column;

		// Loop through string
		int iLast = 0;
		for (int i = 0; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in column name
				if (eReadMode == ReadMode_Column) {

					int j = 0;
					for (; j < vecFormatStrings.size(); j++) {
						if (strSubStr == vecFormatStrings[j]) {
							m_iColumn = j;
							break;
						}
					}
					if (j == vecFormatStrings.size()) {
						_EXCEPTION1("Threshold column name \"%s\" "
							"not found in --format", strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Op;

				// Read in operation
				} else if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = NotEqualTo;
					} else if (strSubStr == "|>=") {
						m_eOp = AbsGreaterThanEqualTo;
					} else if (strSubStr == "|<=") {
						m_eOp = AbsLessThanEqualTo;
					} else {
						_EXCEPTION1("Threshold invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Value;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_MinCount;

				// Read in minimum count
				} else if (eReadMode == ReadMode_MinCount) {
					if (strSubStr == "all") {
						m_nMinimumCount = (-1);
					} else {
						m_nMinimumCount = atoi(strSubStr.c_str());
					}

					if (m_nMinimumCount < -1) {
						_EXCEPTION1("Invalid minimum count \"%i\"",
							m_nMinimumCount);
					}

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("Too many entries in threshold string \"%s\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("Insufficient entries in threshold string \"%s\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription;
		strDescription += vecFormatStrings[m_iColumn];
		if (m_eOp == GreaterThan) {
			strDescription += " greater than ";
		} else if (m_eOp == LessThan) {
			strDescription += " less than ";
		} else if (m_eOp == GreaterThanEqualTo) {
			strDescription += " greater than or equal to ";
		} else if (m_eOp == LessThanEqualTo) {
			strDescription += " less than or equal to ";
		} else if (m_eOp == EqualTo) {
			strDescription += " equal to ";
		} else if (m_eOp == NotEqualTo) {
			strDescription += " not equal to ";
		} else if (m_eOp == AbsGreaterThanEqualTo) {
			strDescription += " magnitude greater than or equal to ";
		} else if (m_eOp == AbsLessThanEqualTo) {
			strDescription += " magnitude less than or equal to ";
		}

		char szValue[128];
		sprintf(szValue, "%f", m_dValue);
		strDescription += szValue;

		char szMinCount[160];
		if (m_nMinimumCount == -1) {
			strDescription += " at all times";
		} else {
			sprintf(szMinCount, " at least %i time(s)", m_nMinimumCount);
			strDescription += szMinCount;
		}

		Announce("%s", strDescription.c_str());
	}

	///	<summary>
	///		Verify that the specified path satisfies the threshold op.
	///	</summary>
	bool Apply(
		const SimplePath & path,
		const std::vector<TimeToCandidateInfoMapIterator> & vecCandidates
	) {
		int nCount = 0;
		for (int s = 0; s < path.m_iTimes.size(); s++) {
			int t = path.m_iTimes[s];
			int i = path.m_iCandidates[s];

			_ASSERT((t >= 0) && (t < vecCandidates.size()));

			double dCandidateValue =
				std::stod(vecCandidates[t]->second[i][m_iColumn]);

			if ((m_eOp == GreaterThan) &&
				(dCandidateValue > m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == LessThan) &&
				(dCandidateValue < m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == GreaterThanEqualTo) &&
				(dCandidateValue >= m_dValue)
			) {
				nCount++;
			
			} else if (
				(m_eOp == LessThanEqualTo) &&
				(dCandidateValue <= m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == EqualTo) &&
				(dCandidateValue == m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == NotEqualTo) &&
				(dCandidateValue != m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == AbsGreaterThanEqualTo) &&
				(fabs(dCandidateValue) >= m_dValue)
			) {
				nCount++;

			} else if (
				(m_eOp == AbsLessThanEqualTo) &&
				(fabs(dCandidateValue) <= m_dValue)
			) {
				nCount++;
			}
		}

		// Check that the criteria is satisfied for all segments
		if (m_nMinimumCount == (-1)) {
			if (nCount == (int)(path.m_iTimes.size())) {
				return true;
			} else {
				return false;
			}
		}

		// Check total count against min count
		if (nCount >= m_nMinimumCount) {
			return true;
		} else {
			return false;
		}
	}

protected:
	///	<summary>
	///		Active column.
	///	</summary>
	int m_iColumn;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Minimum number of segments that need to satisfy the op.
	///	</summary>
	int m_nMinimumCount;
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Input file
	std::string strInputFile;

	// List of input files
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivityFile;

	// Output file
	std::string strOutputFile;

	// Format string
	std::string strFormatOld;

	// Format string
	std::string strFormat;

	// Range (in degrees)
	double dRange;

	// Minimum duration of path
	std::string strMinTime;

	// Minimum path length
	int nMinPathLength;

	// Minimum distance between endpoints of path
	double dMinEndpointDistance;

	// Minimum path length
	double dMinPathDistance;

	// Maximum time gap (in time steps or duration)
	std::string strMaxGapSize;

	// Time step stride
	int nTimeStride;

	// Calendar type
	std::string strCalendar;

	// Allow repeated times
	bool fAllowRepeatedTimes;

	// Output format
	std::string strOutputFileFormat;

	// Thresholds
	std::string strThreshold;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineString(strConnectivityFile, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strFormatOld, "*format", "");
		CommandLineString(strFormat, "in_fmt", "lon,lat");
		CommandLineDoubleD(dRange, "range", 5.0, "(degrees)");
		CommandLineString(strMinTime, "mintime", "3");
		CommandLineInt(nMinPathLength, "*minlength", 1);
		CommandLineDoubleD(dMinEndpointDistance, "min_endpoint_dist", 0.0, "(degrees)");
		CommandLineDoubleD(dMinPathDistance, "min_path_dist", 0.0, "(degrees)");
		CommandLineString(strMaxGapSize, "maxgap", "0");
		CommandLineStringD(strThreshold, "threshold", "",
			"[col,op,value,count;...]");
		CommandLineStringD(strCalendar, "caltype", "standard", "(none|standard|noleap|360_day)");
		CommandLineBool(fAllowRepeatedTimes, "allow_repeated_times");
		//CommandLineInt(nTimeStride, "timestride", 1);
		CommandLineStringD(strOutputFileFormat, "out_file_format", "gfdl", "(gfdl|csv|csvnohead)");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("At least one of --in or --in_list must be specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("At most one of --in or --in_list must be specified");
	}

	// Check format
/*
	std::vector<size_t> nGridDim;
	if (strFormatOld != "") {
		Announce("WARNING: --format is deprecated.  Consider using --in_fmt and --in_connect instead");

		if (strFormatOld.substr(0,4) == "i,j,") {
			strFormat = strFormatOld.substr(4);

			nGridDim.resize(2);
			nGridDim[0] = static_cast<size_t>(-1);
			nGridDim[1] = static_cast<size_t>(-1);

		} else if (strFormatOld.substr(0,2) == "i,") {
			strFormat = strFormatOld.substr(2);

			nGridDim.resize(1);
			nGridDim[0] = static_cast<size_t>(-1);

		} else {
			_EXCEPTIONT("Invalid --format string; expected first entries to be \"i,\" or \"i,j,\"");
		}

	} else {
		if (strConnectivityFile != "") {
			nGridDim.resize(1);
			nGridDim[0] = static_cast<size_t>(-1);
		} else {
			nGridDim.resize(2);
			nGridDim[0] = static_cast<size_t>(-1);
			nGridDim[1] = static_cast<size_t>(-1);
		}
	}
*/
	if (strFormatOld != "") {
		Announce("WARNING: --format is deprecated.  Consider using --in_fmt and --in_connect instead");
		strFormat = strFormatOld;

	} else {
		if (strConnectivityFile != "") {
			strFormat = std::string("i,") + strFormat;
		} else {
			strFormat = std::string("i,j,") + strFormat;
		}
	}

	// Notify that --minlength is deprecated
	if (nMinPathLength != 1) {
		Announce("WARNING: --minlength is deprecated.  Consider using --mintime instead");
	}

	// Parse --mintime
	int nMinTimeSteps = nMinPathLength;
	double dMinTimeSeconds = 0.0;
	if (STLStringHelper::IsIntegerIndex(strMinTime)) {
		nMinTimeSteps = std::stoi(strMinTime);
		if (nMinTimeSteps < 1) {
			_EXCEPTIONT("Invalid value of --mintime; expected positive integer or time delta");
		}

	} else {
		Time timeMinTime;
		timeMinTime.FromFormattedString(strMinTime);
		dMinTimeSeconds = timeMinTime.AsSeconds();
		if (dMinTimeSeconds <= 0.0) {
			_EXCEPTIONT("Invalid value for --mintime; expected positive integer or positive time delta");
		}
	}

	// Parse --maxgap
	int nMaxGapSteps = (-1);
	double dMaxGapSeconds = -1.0;
	if (STLStringHelper::IsIntegerIndex(strMaxGapSize)) {
		nMaxGapSteps = std::stoi(strMaxGapSize);
		if (nMaxGapSteps < 0) {
			_EXCEPTIONT("Invalid value of --maxgap; expected nonnegative integer or time delta");
		}

	} else {
		Time timeMaxGap;
		timeMaxGap.FromFormattedString(strMaxGapSize);
		dMaxGapSeconds = timeMaxGap.AsSeconds();
		if (dMaxGapSeconds < 0.0) {
			_EXCEPTIONT("Invalid value for --maxgap; expected nonnegative integer or positive time delta");
		}
	}

	// Parse calendar type
	Time::CalendarType caltype = Time::CalendarTypeFromString(strCalendar);

	if (caltype == Time::CalendarNone) {
		if (dMinTimeSeconds > 0.0) {
			_EXCEPTIONT("A calendar type (--caltype) must be specified if using time deltas for --mintime");
		}
		if (dMaxGapSeconds >= 0.0) {
			_EXCEPTIONT("A calendar type (--caltype) must be specified if using time deltas for --maxgap");
		}
	}

	// Output format
	if ((strOutputFileFormat != "gfdl") &&
		(strOutputFileFormat != "csv") &&
		(strOutputFileFormat != "csvnohead")
	) {
		_EXCEPTIONT("Output format must be either \"gfdl\", \"csv\", or \"csvnohead\"");
	}

	// Input file list
	FilenameList vecInputFiles;

	if (strInputFile != "") {
		vecInputFiles.push_back(strInputFile);
	}
	if (strInputFileList != "") {
		vecInputFiles.FromFile(strInputFileList, false);
	}

	// Parse format string
	std::vector< std::string > vecFormatStrings;
	ParseVariableList(strFormat, vecFormatStrings);

	if (vecFormatStrings.size() == 0) {
		_EXCEPTIONT("No format specified");
	}

	// Check format string for lat/lon
	int iLatIndex = (-1);
	int iLonIndex = (-1);

	for (int i = 0; i < vecFormatStrings.size(); i++) {
		if (vecFormatStrings[i] == "lat") {
			iLatIndex = i;
		}
		if (vecFormatStrings[i] == "lon") {
			iLonIndex = i;
		}
	}

	if (iLatIndex == (-1)) {
		_EXCEPTIONT("Latitude \"lat\" must be specified in format");
	}
	if (iLonIndex == (-1)) {
		_EXCEPTIONT("Longitude \"lon\" must be specified in format");
	}

	// Parse the threshold string
	std::vector<PathThresholdOp> vecThresholdOp;

	if (strThreshold != "") {
		AnnounceStartBlock("Parsing thresholds");

		int iLast = 0;
		for (int i = 0; i <= strThreshold.length(); i++) {

			if ((i == strThreshold.length()) || (strThreshold[i] == ';')) {
				std::string strSubStr = strThreshold.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecThresholdOp.size());
				vecThresholdOp.resize(iNextOp + 1);
				vecThresholdOp[iNextOp].Parse(strSubStr, vecFormatStrings);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the input
	std::vector<TimeToCandidateInfoMapIterator> vecIterCandidates;
	TimeToCandidateInfoMap mapCandidateInfo;
	{
		AnnounceStartBlock("Loading candidate data");

		// Parse all input files
		for (size_t f = 0; f < vecInputFiles.size(); f++) {
			if (vecInputFiles.size() > 1) {
				Announce("File (%lu/%lu) \"%s\"", f+1, vecInputFiles.size(), vecInputFiles[f].c_str());
			}
			ParseDetectNodesFile(
				vecInputFiles[f],
				vecFormatStrings,
				mapCandidateInfo,
				caltype,
				fAllowRepeatedTimes);
		}

		if (mapCandidateInfo.size() == 0) {
			_EXCEPTIONT("Zero candidate nodes found within input files; nothing to stitch");
		}

		vecIterCandidates.reserve(mapCandidateInfo.size());
		for (auto it = mapCandidateInfo.begin(); it != mapCandidateInfo.end(); it++) {
			vecIterCandidates.push_back(it);
		}

		_ASSERT(vecIterCandidates.size() == mapCandidateInfo.size());

		const Time & timeFirst = vecIterCandidates[0]->first;
		const Time & timeLast = vecIterCandidates[vecIterCandidates.size()-1]->first;

		Announce("Discrete times: %i (%s to %s)",
			vecIterCandidates.size(),
			timeFirst.ToString().c_str(),
			timeLast.ToString().c_str());

		// Verify total time is larger than mintime
		if (dMinTimeSeconds > 0.0) {
			double dDeltaSeconds = timeFirst.DeltaSeconds(timeLast);

			if (dMinTimeSeconds > dDeltaSeconds) {
				_EXCEPTION3("Duration of record (%id%ih%i) is shorter than --mintime; no paths possible",
					static_cast<int>(dDeltaSeconds) / 86400,
					(static_cast<int>(dDeltaSeconds) % 86400) / 3600,
					static_cast<int>(dDeltaSeconds) % 3600);
			}
		}

		// Verify adjacent times are smaller than maxgap
		if (dMaxGapSeconds >= 0.0) {
			for (size_t t = 0; t < vecIterCandidates.size()-1; t++) {
				const Time & timeCurrent = vecIterCandidates[t]->first;
				const Time & timeNext = vecIterCandidates[t+1]->first;

				double dDeltaSeconds = timeCurrent.DeltaSeconds(timeNext);
				if (dDeltaSeconds > dMaxGapSeconds) {
					Announce("WARNING: Discrete times %i (%s) and %i (%s) differ by more than --maxgap (%s)",
						t,
						timeCurrent.ToString().c_str(),
						t+1,
						timeNext.ToString().c_str(),
						strMaxGapSize.c_str());
				}
			}
		}

		AnnounceEndBlock("Done");
	}

	/////////////////////////////////////////////
	// 	Build KD trees at each time slice

	AnnounceStartBlock("Creating KD trees at each time level");

	// Null pointer
	int * noptr = NULL;

	// Vector of lat/lon values
	std::vector< std::vector<Node> > vecNodes;
	vecNodes.resize(mapCandidateInfo.size());

	// Vector of KD trees
	std::vector<kdtree *> vecKDTrees;
	vecKDTrees.resize(mapCandidateInfo.size());

	for (size_t t = 0; t < vecIterCandidates.size(); t++) {

		const TimesliceCandidateInformation & tscinfo = vecIterCandidates[t]->second;

		// Create a new kdtree
		if (tscinfo.size() == 0) {
			vecKDTrees[t] = NULL;
			continue;
		}

		vecKDTrees[t] = kd_create(3);

		vecNodes[t].resize(tscinfo.size());

		// Insert all points at this time level
		for (int i = 0; i < tscinfo.size(); i++) {
			vecNodes[t][i].lonRad = DegToRad(std::stod(tscinfo[i][iLonIndex].c_str()));
			vecNodes[t][i].latRad = DegToRad(std::stod(tscinfo[i][iLatIndex].c_str()));

			RLLtoXYZ_Rad(
				vecNodes[t][i].lonRad,
				vecNodes[t][i].latRad,
				vecNodes[t][i].x,
				vecNodes[t][i].y,
				vecNodes[t][i].z);

			kd_insert3(
				vecKDTrees[t],
				vecNodes[t][i].x,
				vecNodes[t][i].y,
				vecNodes[t][i].z,
				reinterpret_cast<void*>(noptr+i));
		}
	}

	AnnounceEndBlock("Done");

	/////////////////////////////////////////////
	// 	Populate the set of path segments

	AnnounceStartBlock("Populating set of path segments");

	std::vector<SimplePathSegmentSet> vecPathSegmentsSet;
	vecPathSegmentsSet.resize(mapCandidateInfo.size()-1);

	// Search nodes at current time level
	for (size_t t = 0; t < vecIterCandidates.size()-1; t++) {

		const Time & timeCurrent = vecIterCandidates[t]->first;
		const TimesliceCandidateInformation & tscinfo = vecIterCandidates[t]->second;

		// Loop through all candidates at current time level
		for (int i = 0; i < tscinfo.size(); i++) {

			// Check future timesteps for next candidate in path
			for (int g = 1; ; g++) {

				if ((nMaxGapSteps != (-1)) && (g > nMaxGapSteps+1)) {
					break;
				}
				if (t+g >= vecIterCandidates.size()) {
					break;
				}
		
				if (dMaxGapSeconds >= 0.0) {
					const Time & timeNext = vecIterCandidates[t+g]->first;
					double dDeltaSeconds = timeCurrent.DeltaSeconds(timeNext);
					if (dDeltaSeconds > dMaxGapSeconds) {
						break;
					}
				}

				if (vecKDTrees[t+g] == NULL) {
					continue;
				}

				kdres * set = kd_nearest3(vecKDTrees[t+g], vecNodes[t][i].x, vecNodes[t][i].y, vecNodes[t][i].z);

				if (kd_res_size(set) == 0) {
					kd_res_free(set);
					break;
				}

				int iRes =
					  reinterpret_cast<int*>(kd_res_item_data(set))
					- reinterpret_cast<int*>(noptr);

				kd_res_free(set);

				// Great circle distance between points
				double dR =
					GreatCircleDistance_Deg(
						vecNodes[t][i].lonRad,
						vecNodes[t][i].latRad,
						vecNodes[t+g][iRes].lonRad,
						vecNodes[t+g][iRes].latRad);

				// Verify great circle distance satisfies range requirement
				if (dR <= dRange) {

					// Insert new path segment into set of path segments
					vecPathSegmentsSet[t].insert(
						SimplePathSegment(t, i, t+g, iRes));

					break;
				}
			}
		}
	}

	AnnounceEndBlock("Done");

	/////////////////////////////////////////////
	// 	Find all paths

	AnnounceStartBlock("Constructing paths");

	std::vector<SimplePath> vecPaths;

	int nRejectedMinTimePaths = 0;
	int nRejectedMinEndpointDistPaths = 0;
	int nRejectedMinPathDistPaths = 0;
	int nRejectedThresholdPaths = 0;

	// Remove path
	for (size_t t = 0; t < vecIterCandidates.size()-1; t++) {

		// Loop through all remaining segments
		while (vecPathSegmentsSet[t].size() > 0) {

			// Create a new path
			SimplePath path;

			SimplePathSegmentSet::iterator iterSeg
				= vecPathSegmentsSet[t].begin();

			path.m_iTimes.push_back(iterSeg->m_iTime[0]);
			path.m_iCandidates.push_back(iterSeg->m_iCandidate[0]);

			int tx = t;

			for (;;) {
				path.m_iTimes.push_back(iterSeg->m_iTime[1]);
				path.m_iCandidates.push_back(iterSeg->m_iCandidate[1]);

				int txnext = iterSeg->m_iTime[1];

				if (txnext >= vecIterCandidates.size()-1) {
					vecPathSegmentsSet[tx].erase(iterSeg);
					break;
				}

				SimplePathSegment segFind(
					iterSeg->m_iTime[1], iterSeg->m_iCandidate[1], 0, 0);

				vecPathSegmentsSet[tx].erase(iterSeg);

				iterSeg = vecPathSegmentsSet[txnext].find(segFind);

				if (iterSeg == vecPathSegmentsSet[txnext].end()) {
					break;
				}

				tx = txnext;
			}

			// Reject path due to minimum time
			if (path.m_iTimes.size() < nMinTimeSteps) {
				nRejectedMinTimePaths++;
				continue;
			}
			if (dMinTimeSeconds > 0.0) {
				int nT = path.m_iTimes.size();
				int iTime0 = path.m_iTimes[0];
				int iTime1 = path.m_iTimes[nT-1];

				_ASSERT((iTime0 >= 0) && (iTime0 < vecIterCandidates.size()));
				_ASSERT((iTime1 >= 0) && (iTime1 < vecIterCandidates.size()));

				const Time & timeFirst = vecIterCandidates[iTime0]->first;
				const Time & timeLast = vecIterCandidates[iTime1]->first;

				double dDeltaSeconds = timeFirst.DeltaSeconds(timeLast);
				if (dDeltaSeconds < dMinTimeSeconds) {
					nRejectedMinTimePaths++;
					continue;
				}
			}

			// Reject path due to minimum endpoint distance
			if (dMinEndpointDistance > 0.0) {
				int nT = path.m_iTimes.size();

				int iTime0 = path.m_iTimes[0];
				int iRes0  = path.m_iCandidates[0];

				int iTime1 = path.m_iTimes[nT-1];
				int iRes1  = path.m_iCandidates[nT-1];

				double dLonRad0 = vecNodes[iTime0][iRes0].lonRad;
				double dLatRad0 = vecNodes[iTime0][iRes0].latRad;

				double dLonRad1 = vecNodes[iTime1][iRes1].lonRad;
				double dLatRad1 = vecNodes[iTime1][iRes1].latRad;

				double dR = GreatCircleDistance_Deg(dLonRad0, dLatRad0, dLonRad1, dLatRad1);

				if (dR < dMinEndpointDistance) {
					nRejectedMinEndpointDistPaths++;
					continue;
				}
			}

			// Reject path due to minimum total path distance
			if (dMinPathDistance > 0.0) {
				double dTotalPathDistance = 0.0;
				for (int i = 0; i < path.m_iTimes.size() - 1; i++) {
					int iTime0 = path.m_iTimes[i];
					int iRes0 = path.m_iCandidates[i];

					int iTime1 = path.m_iTimes[i+1];
					int iRes1 = path.m_iCandidates[i+1];

					double dLonRad0 = vecNodes[iTime0][iRes0].lonRad;
					double dLatRad0 = vecNodes[iTime0][iRes0].latRad;

					double dLonRad1 = vecNodes[iTime1][iRes1].lonRad;
					double dLatRad1 = vecNodes[iTime1][iRes1].latRad;

					double dR = GreatCircleDistance_Deg(dLonRad0, dLatRad0, dLonRad1, dLatRad1);

					dTotalPathDistance += dR;
				}

				if (dTotalPathDistance < dMinPathDistance) {
					nRejectedMinPathDistPaths++;
					continue;
				}
			}

			// Reject path due to threshold
			bool fOpResult = true;
			for (int x = 0; x < vecThresholdOp.size(); x++) {
				fOpResult =
					vecThresholdOp[x].Apply(
						path,
						vecIterCandidates);

				if (!fOpResult) {
					break;
				}
			}
			if (!fOpResult) {
				nRejectedThresholdPaths++;
				continue;
			}

			// Add path to array of paths
			vecPaths.push_back(path);
		}
	}

	Announce("Paths rejected (mintime): %i", nRejectedMinTimePaths);
	Announce("Paths rejected (minendpointdist): %i", nRejectedMinEndpointDistPaths);
	Announce("Paths rejected (minpathdist): %i", nRejectedMinPathDistPaths);
	Announce("Paths rejected (threshold): %i", nRejectedThresholdPaths);
	Announce("Total paths found: %i", vecPaths.size());
	AnnounceEndBlock("Done");

	/////////////////////////////////////////////
	// 	Write results

	AnnounceStartBlock("Writing results");

	// Convert the result into a NodeFile
	{
		NodeFile nodefile;
		nodefile.m_ePathType = NodeFile::PathTypeSN;
		nodefile.m_pathvec.resize(vecPaths.size());

		for (int i = 0; i < vecFormatStrings.size(); i++) {
			nodefile.m_cdh.push_back(vecFormatStrings[i]);
		}

		for (int i = 0; i < vecPaths.size(); i++) {
			Path & path = nodefile.m_pathvec[i];
			path.resize(vecPaths[i].m_iTimes.size());

			if (vecPaths[i].m_iTimes.size() == 0) {
				_EXCEPTIONT("Zero length Path found");
			}

			for (int t = 0; t < vecPaths[i].m_iTimes.size(); t++) {
				PathNode & pathnode = path[t];

				int iTime = vecPaths[i].m_iTimes[t];

				int iCandidate = vecPaths[i].m_iCandidates[t];

				path[t].m_time = vecIterCandidates[iTime]->first;

				for (int j = 0; j < vecIterCandidates[iTime]->second[iCandidate].size(); j++) {
					pathnode.m_vecColumnData.push_back(
						new ColumnDataString(
							vecIterCandidates[iTime]->second[iCandidate][j]));
				}
			}

			path.m_timeStart = path[0].m_time;
		}

		if (strOutputFileFormat == "gfdl") {
			nodefile.Write(strOutputFile);
		} else if (strOutputFileFormat == "csv") {
			nodefile.Write(strOutputFile, NULL, NULL, NodeFile::FileFormatCSV, true);
		} else if (strOutputFileFormat == "csvnohead") {
			nodefile.Write(strOutputFile, NULL, NULL, NodeFile::FileFormatCSV, false);
		}
	}
/*
	if (strOutputFileFormat == "std") {
		FILE * fp = fopen(strOutputFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Failed to open output file \"%s\"",
				strOutputFile.c_str());
		}

		for (int i = 0; i < vecPaths.size(); i++) {
			int iStartTime = vecPaths[i].m_iTimes[0];

			fprintf(fp, "start");
			fprintf(fp, "\t%li", vecPaths[i].m_iTimes.size());

			int jEnd = vecTimes[iStartTime].size();
			if (jEnd > 5) {
				jEnd = 5;
			}
			for (int j = 0; j < jEnd; j++) {
				if (j == 3) {
					continue;
				}
				fprintf(fp, "\t%s", vecTimes[iStartTime][j].c_str());
			}
			fprintf(fp, "\n");

			for (int t = 0; t < vecPaths[i].m_iTimes.size(); t++) {
				int iTime = vecPaths[i].m_iTimes[t];
				int iCandidate = vecPaths[i].m_iCandidates[t];

				for (int j = 0; j < vecCandidates[iTime][iCandidate].size(); j++) {
					fprintf(fp, "\t%s",
						vecCandidates[iTime][iCandidate][j].c_str());
				}
				for (int j = 0; j < jEnd; j++) {
					if (j == 3) {
						continue;
					}
					fprintf(fp, "\t%s", vecTimes[iTime][j].c_str());
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);

	} else if (strOutputFileFormat == "visit") {
		FILE * fp = fopen(strOutputFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Failed to open output file \"%s\"",
				strOutputFile.c_str());
		}

		// Write output format
		fprintf(fp, "#id,time_id,year,month,day,hour,");
		fprintf(fp, "%s", strFormat.c_str());
		fprintf(fp, "\n");

		for (int i = 0; i < vecPaths.size(); i++) {
			for (int t = 0; t < vecPaths[i].m_iTimes.size(); t++) {
				int iTime = vecPaths[i].m_iTimes[t];
				int iCandidate = vecPaths[i].m_iCandidates[t];

				fprintf(fp, "%i,\t%i,\t%s,\t%s,\t%s,\t%s,\t",
					i+1, t+1,
					vecTimes[iTime][2].c_str(),
					vecTimes[iTime][1].c_str(),
					vecTimes[iTime][0].c_str(),
					vecTimes[iTime][4].c_str());

				fprintf(fp, "\t");
				for (int j = 0; j < vecCandidates[iTime][iCandidate].size(); j++) {
					fprintf(fp, "%s",
						vecCandidates[iTime][iCandidate][j].c_str());
					if (j != vecCandidates[iTime][iCandidate].size()-1) {
						fprintf(fp, ",\t");
					}
				}
				fprintf(fp, "\n");
			}
		}
	}
*/
	AnnounceEndBlock("Done");

	// Cleanup
	AnnounceStartBlock("Cleanup");

	for (int t = 0; t < vecKDTrees.size(); t++) {
		kd_free(vecKDTrees[t]);
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

