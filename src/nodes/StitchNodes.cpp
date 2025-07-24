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
#include "TimeMatch.h"

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

///	<summary>
///		Read a DetectNodes output file into mapCandidateInfo.
///	</summary>
void ParseDetectNodesFile(
	const std::string & strInputFile,
	const std::vector< std::string > & vecFormatStrings,
	TimeToCandidateInfoMap & mapCandidateInfo,
	Time::CalendarType caltype,
	bool fAllowRepeatedTimes,
	const Time & timeBegin,
	const Time & timeEnd,
	size_t sGridDims,
	const std::regex & reTimeSubset
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
	size_t sFormatEntries = vecFormatStrings.size();

	if (sFormatEntries < sGridDims) {
		_EXCEPTION2("Number of format entries (%lu) must be greater than number of grid dimensions (%lu)",
			sFormatEntries, sGridDims);
	}

	// Current candidate at this time
	int iCandidate = 0;

	// Current candidate information
	TimeToCandidateInfoMapIterator iterCurrentTime = mapCandidateInfo.end();

	// Line number
	size_t sLineNumber = 0;

	// Total candidates at this time
	int nCandidates = 0;

	// Do not include this time slice
	bool fIncludeTimeSlice = true;

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
		sLineNumber++;

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
				_EXCEPTION2("Malformed time string on line (%lu):\n%s", sLineNumber, strLine.c_str());
			}

			iCandidate = 0;
			nCandidates = std::stoi(vecTimeString[3]);

			Time time(caltype);
			time.SetYear(std::stoi(vecTimeString[0]));
			time.SetMonth(std::stoi(vecTimeString[1]));
			time.SetDay(std::stoi(vecTimeString[2]));

			if (vecTimeString[4].length() == 5) {
				time.SetSecond(std::stoi(vecTimeString[4]));
			} else {
				time.SetSecond(std::stoi(vecTimeString[4]) * 3600);
			}

			fIncludeTimeSlice = true;
			if (timeBegin.GetCalendarType() != Time::CalendarUnknown) {
				if (time < timeBegin) {
					fIncludeTimeSlice = false;
				}
			}
			if (timeEnd.GetCalendarType() != Time::CalendarUnknown) {
				if (time > timeEnd) {
					fIncludeTimeSlice = false;
				}
			}
			
#ifndef TEMPEST_NOREGEX
			std::smatch match;
			if (!std::regex_search(vecTimeString[4], match, reTimeSubset)) {
				fIncludeTimeSlice = false;
			}
#endif
			if (fIncludeTimeSlice) {
				auto it = mapCandidateInfo.find(time);
				if (it != mapCandidateInfo.end()) {
					if (fAllowRepeatedTimes) {
						iterCurrentTime = it;
						iCandidate = iterCurrentTime->second.size();
						nCandidates += iCandidate;
						iterCurrentTime->second.resize(iCandidate + nCandidates);

					} else {
						_EXCEPTION2("Repeated time \"%s\" found in candidate files on line (%lu)",
							time.ToString().c_str(), sLineNumber);
					}

				} else {
					auto ins =
						mapCandidateInfo.insert(
							TimeToCandidateInfoPair(
								time, TimesliceCandidateInformation(nCandidates)));

					if (!ins.second) {
						_EXCEPTION2("Insertion of time \"%s\" into candidate info map failed on line (%lu)",
							time.ToString().c_str(), sLineNumber);
					}

					iterCurrentTime = ins.first;
				}
			}

			// Prepare to parse candidate data
			if (nCandidates != 0) {
				eReadState = ReadState_Candidate;
			}

		// Parse candidate information
		} else if (eReadState == ReadState_Candidate) {

			if (fIncludeTimeSlice) {

				// Parse candidates
				_ASSERT(iterCurrentTime != mapCandidateInfo.end());

				TimesliceCandidateInformation & tscinfo = iterCurrentTime->second;

				_ASSERT(iCandidate < tscinfo.size());

				ParseVariableList(strLine, tscinfo[iCandidate]);

				if (tscinfo[iCandidate].size() != sFormatEntries) {
					fWarnInsufficientCandidateInfo = true;
				}

				if (tscinfo[iCandidate].size() < sGridDims) {
					_EXCEPTION2("On line %lu there are too few columns in data; first %lu columns must be integer indices", sLineNumber, sGridDims);
				}
				for (size_t d = 0; d < sGridDims; d++) {
					if (!STLStringHelper::IsIntegerIndex(tscinfo[iCandidate][d])) {
						_EXCEPTION2("On line %lu first %lu columns are not integer indices (equal to number of grid dimensions); "
							"if an unstructured grid is being used make sure --in_connect is specified.",
							sLineNumber, sGridDims);
					}
				}
			}

			// Advance candidate number
			iCandidate++;
			if (iCandidate == nCandidates) {
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

///	<summary>
///		A (time,candidate index) pair
///	</summary>
struct TimeCandidatePair : public std::pair<int,int> {
	TimeCandidatePair(int _timeix, int _candidate) :
		std::pair<int,int>(_timeix, _candidate)
	{ }

	inline const int & timeix() const {
		return first;
	}

	inline const int & candidate() const {
		return second;
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A segment in a simple path connecting candidate indices at two times.
///	</summary>
class SimplePathSegment :
	public std::pair<TimeCandidatePair, TimeCandidatePair>
{

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SimplePathSegment(
		int iTime0,
		int iCandidate0,
		int iTime1,
		int iCandidate1
	) :
		std::pair<TimeCandidatePair, TimeCandidatePair>(
			TimeCandidatePair(iTime0, iCandidate0),
			TimeCandidatePair(iTime1, iCandidate1))
	{ }

public:
	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator< (const SimplePathSegment & seg) const {
		return (first < seg.first);
	}
};

typedef std::set<SimplePathSegment> SimplePathSegmentSet;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of times and candidate indices at those times defining a path.
///	</summary>
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

///	<summary>
///		A vector of SimplePaths.
///	</summary>
typedef std::vector<SimplePath> SimplePathVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An operator which is used to enforce some criteria on paths.
///	</summary>
class PathThresholdOp {

public:
	///	<summary>
	///		Index denoting all points along path.
	///	</summary>
	static const int Count_All = (-1);

	///	<summary>
	///		Index denoting only the first point along the path.
	///	</summary>
	static const int Count_First = (-2);

	///	<summary>
	///		Index denoting only the last point along the path.
	///	</summary>
	static const int Count_Last = (-3);

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
						m_nMinimumCount = Count_All;
					} else if (strSubStr == "first") {
						m_nMinimumCount = Count_First;
					} else if (strSubStr == "last") {
						m_nMinimumCount = Count_Last;
					} else {
						m_nMinimumCount = atoi(strSubStr.c_str());
					}

					if (m_nMinimumCount < Count_Last) {
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
		snprintf(szValue, 128, "%f", m_dValue);
		strDescription += szValue;

		char szMinCount[160];
		if (m_nMinimumCount == Count_All) {
			strDescription += " at all times";
		} else if (m_nMinimumCount == Count_First) {
			strDescription += " at the first time";
		} else if (m_nMinimumCount == Count_Last) {
			strDescription += " at the last time";
		} else {
			snprintf(szMinCount, 160, " at least %i time(s)", m_nMinimumCount);
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

			if ((m_nMinimumCount == Count_First) && (s > 0)) {
				continue;
			}
			if ((m_nMinimumCount == Count_Last) && (s < path.m_iTimes.size()-1)) {
				continue;
			}

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
		if (m_nMinimumCount == Count_All) {
			if (nCount == (int)(path.m_iTimes.size())) {
				return true;
			} else {
				return false;
			}
		}
	
		// Check that the criteria are satisfied for the first or last segment
		if ((m_nMinimumCount == Count_First) || (m_nMinimumCount == Count_Last)) {
			if (nCount == 1) {
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

void GeneratePathSegmentsSetBasic(
	const std::vector<TimeToCandidateInfoMapIterator> & vecIterCandidates,
	const std::vector< std::vector<Node> > & vecNodes,
	const std::vector<kdtree *> & vecKDTrees,
	int nMaxGapSteps,
	double dMaxGapSeconds,
	double dRangeDeg,
	std::vector<SimplePathSegmentSet> & vecPathSegmentsSet
) {

	vecPathSegmentsSet.resize(vecIterCandidates.size()-1);

	// Null pointer
	int * noptr = NULL;

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
				double dRDeg =
					GreatCircleDistance_Deg(
						vecNodes[t][i].lonRad,
						vecNodes[t][i].latRad,
						vecNodes[t+g][iRes].lonRad,
						vecNodes[t+g][iRes].latRad);
                // printf("%1.15e %1.15e\n", dRDeg, dRangeDeg);
				// Verify great circle distance satisfies range requirement

				if (dRDeg <= dRangeDeg) {
					// Insert new path segment into set of path segments
					vecPathSegmentsSet[t].insert(
						SimplePathSegment(t, i, t+g, iRes));

					break;
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		First-come-first-serve method for generating paths (default).
///	</summary>
void GeneratePathsBasic(
	const std::vector<TimeToCandidateInfoMapIterator> & vecIterCandidates,
	std::vector<SimplePathSegmentSet> & vecPathSegmentsSet,
	std::vector<SimplePath> & vecPaths
) {
	// Loop through all times
	for (size_t t = 0; t < vecIterCandidates.size()-1; t++) {

		// Loop through all remaining segments at this time
		while (vecPathSegmentsSet[t].size() > 0) {

			// Create a new path starting with this segment
			SimplePath path;

			SimplePathSegmentSet::iterator iterSeg
				= vecPathSegmentsSet[t].begin();

			path.m_iTimes.push_back(iterSeg->first.timeix());
			path.m_iCandidates.push_back(iterSeg->first.candidate());

			int tx = t;

			for (;;) {
				path.m_iTimes.push_back(iterSeg->second.timeix());
				path.m_iCandidates.push_back(iterSeg->second.candidate());

				int txnext = iterSeg->second.timeix();

				if (txnext >= vecIterCandidates.size()-1) {
					vecPathSegmentsSet[tx].erase(iterSeg);
					break;
				}

				SimplePathSegment segFind(
					iterSeg->second.timeix(), iterSeg->second.candidate(), 0, 0);

				vecPathSegmentsSet[tx].erase(iterSeg);

				iterSeg = vecPathSegmentsSet[txnext].find(segFind);

				if (iterSeg == vecPathSegmentsSet[txnext].end()) {
					break;
				}

				tx = txnext;
			}

			// Add path to array of paths
			vecPaths.push_back(path);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Priority method for generating path segments.
///	</summary>
void GeneratePathSegmentsWithPriority(
	const std::vector<TimeToCandidateInfoMapIterator> & vecIterCandidates,
	const std::vector< std::vector<Node> > & vecNodes,
	const std::vector<kdtree *> & vecKDTrees,
	const int ixPriorityCol,
	int nMaxGapSteps,
	double dMaxGapSeconds,
	double dRangeDeg,
	std::vector<SimplePathSegmentSet> & vecPathSegmentsSet
) {
	// For target candidates, priority includes both timestep delta and distance
	typedef std::pair<int, double> PriorityPair;

	// Null pointer
	int * noptr = NULL;

	// Chord distance
	const double dChordLength = ChordLengthFromGreatCircleDistance_Deg(dRangeDeg);

	// Initialize the path segment set
	vecPathSegmentsSet.resize(vecIterCandidates.size()-1);

	// Get priority of all candidates
	if (ixPriorityCol < 0) {
		_EXCEPTIONT("Invalid priority column index");
	}
	std::multimap<double, TimeCandidatePair> mapPriority;
	for (size_t t = 0; t < vecIterCandidates.size()-1; t++) {
		const TimesliceCandidateInformation & tscinfo = vecIterCandidates[t]->second;

		for (int i = 0; i < tscinfo.size(); i++) {
			if (tscinfo[i].size() <= ixPriorityCol) {
				_EXCEPTION2("Priority column index out of range (%i <= %i)", tscinfo[i].size(), ixPriorityCol);
			}

			double dPriority = std::stod(tscinfo[i][ixPriorityCol]);
			mapPriority.insert(std::pair<double, TimeCandidatePair>(dPriority, TimeCandidatePair(t,i)));
		}
	}

	Announce("%lu total candidates found", mapPriority.size());

	// Set of candidates that are already targets of another node
	std::set<TimeCandidatePair> setUsedCandidates;

	// Loop through all candidates in increasing priority
	for (auto & it : mapPriority) {
		const int t = it.second.timeix();
		const int i = it.second.candidate();

		const Time & timeCurrent = vecIterCandidates[t]->first;

		// Find all candidates within prescribed distance
		std::multimap<std::pair<int, double>, TimeCandidatePair> mmapTargets;
		for (int g = 1; ; g++) {

			if ((nMaxGapSteps != (-1)) && (g > nMaxGapSteps+1)) {
				break;
			}
			if (t+g >= vecKDTrees.size()) {
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

			kdres * set =
				kd_nearest_range3(
					vecKDTrees[t+g],
					vecNodes[t][i].x,
					vecNodes[t][i].y,
					vecNodes[t][i].z,
					dChordLength);

			if (set == NULL) {
				_EXCEPTIONT("Fatal exception in kd_nearest_range3");
			}
			if (kd_res_size(set) == 0) {
				kd_res_free(set);
				continue;
			}

			do {
				int iRes =
					  reinterpret_cast<int*>(kd_res_item_data(set))
					- reinterpret_cast<int*>(noptr);

				//printf("%i %i %i %i - %lu %lu\n", t, i, t+g, iRes);

				double dRDeg =
					GreatCircleDistance_Deg(
						vecNodes[t][i].lonRad,
						vecNodes[t][i].latRad,
						vecNodes[t+g][iRes].lonRad,
						vecNodes[t+g][iRes].latRad);

				mmapTargets.insert(
					std::pair<PriorityPair, TimeCandidatePair>(
						PriorityPair(g,dRDeg),
						TimeCandidatePair(t+g,iRes)));

			} while (kd_res_next(set));

			kd_res_free(set);
		}

		// Find the shortest target node and insert a segment
		for (auto & itTarget : mmapTargets) {
			auto itUsed = setUsedCandidates.find(itTarget.second);
			if (itUsed == setUsedCandidates.end()) {
				setUsedCandidates.insert(itTarget.second);

				_ASSERT((t >= 0) && (t < vecPathSegmentsSet.size()));
				vecPathSegmentsSet[t].insert(
					SimplePathSegment(t, i, itTarget.second.timeix(), itTarget.second.candidate()));

				break;
			}
		}
	}

/*
	// Build the reverse candidate map, which identifies, for each candidate,
	// all candidates pointing to it.
	typedef std::vector<TimeCandidatePair> TimeCandidateVector;
	typedef std::map<TimeCandidatePair, TimeCandidateVector> ReverseCandidateMap;

	ReverseCandidateMap mapReverseCandidate;

	for (auto & setPathSegmentsAtTime : vecPathSegmentsSet) {
		for (auto & aSimplePathSegment : setPathSegmentsAtTime) {
			auto it = mapReverseCandidate.find(aSimplePathSegment.second);
			if (it == mapReverseCandidate.end()) {
				std::pair<ReverseCandidateMap::iterator, bool> prResult =
					mapReverseCandidate.insert(
						ReverseCandidateMap::value_type(
							aSimplePathSegment.second,
							TimeCandidateVector()));

				prResult.first->second.push_back(
					aSimplePathSegment.first);

			} else {
				it->second.push_back(
					aSimplePathSegment.first);
			}
		}
	}

	// Loop through map and identify priority candidate
	std::set<TimeCandidatePair> setConsumedCandidates;

	// Only retain that segment in array
*/
}

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
	double dRangeDeg;

	// Minimum duration of path
	std::string strMinTime;

	// Begin time for analysis
	std::string strTimeBegin;

	// End time for analysis
	std::string strTimeEnd;

	// Variable to use when prioritizing paths
	std::string strPrioritize;

	// Minimum path length
	int nMinPathLength;

	// Minimum distance between endpoints of path
	double dMinEndpointDistance;

	// Minimum path length
	double dMinPathDistance;

	// Maximum time gap (in time steps or duration)
	std::string strMaxGapSize;

	// Time subset
	std::string strTimeFilter;

	// Calendar type
	std::string strCalendar;

	// Allow repeated times
	bool fAllowRepeatedTimes;

	// Add relative velocity to the output
	bool fAddVelocity;

	// Output format
	std::string strOutputFileFormat;

	// Thresholds
	std::string strThreshold;

	// Output in seconds
	bool fOutputSeconds;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineString(strConnectivityFile, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strFormatOld, "*format", "");
		CommandLineString(strFormat, "in_fmt", "lon,lat");
		CommandLineDoubleD(dRangeDeg, "range", 5.0, "(degrees)");
		//TODO: CommandLineDoubleD(dSpeed, "maxspeed", 0.0, "(degrees/hour)");
		CommandLineString(strMinTime, "mintime", "3");
		CommandLineInt(nMinPathLength, "*minlength", 1);
		CommandLineString(strTimeBegin, "time_begin", "");
		CommandLineString(strTimeEnd, "time_end", "");
		CommandLineString(strPrioritize, "prioritize", "");
		CommandLineDoubleD(dMinEndpointDistance, "min_endpoint_dist", 0.0, "(degrees)");
		CommandLineDoubleD(dMinPathDistance, "min_path_dist", 0.0, "(degrees)");
		CommandLineString(strMaxGapSize, "maxgap", "0");
		CommandLineStringD(strThreshold, "threshold", "",
			"[col,op,value,count;...]");
		CommandLineStringD(strCalendar, "caltype", "standard", "(none|standard|noleap|360_day)");
		CommandLineString(strTimeFilter, "timefilter", "");
		CommandLineBool(fAllowRepeatedTimes, "allow_repeated_times");
		CommandLineBool(fAddVelocity, "add_velocity");
		//CommandLineInt(nTimeStride, "timestride", 1);
		CommandLineStringD(strOutputFileFormat, "out_file_format", "gfdl", "(gfdl|csv|csvnohead)");
		CommandLineBool(fOutputSeconds, "out_seconds");

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
	size_t sGridDims = 0;
	if (strFormatOld != "") {
		Announce("WARNING: --format is deprecated.  Consider using --in_fmt and --in_connect instead");
		strFormat = strFormatOld;

	} else {
		if (strConnectivityFile != "") {
			strFormat = std::string("i,") + strFormat;
			sGridDims = 1;
		} else {
			strFormat = std::string("i,j,") + strFormat;
			sGridDims = 2;
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

	// Parse --time_begin and --time_end
	Time timeBegin(Time::CalendarUnknown);
	Time timeEnd(Time::CalendarUnknown);
	if (strTimeBegin != "") {
		timeBegin = Time(caltype);
		timeBegin.FromFormattedString(strTimeBegin);
		if (timeBegin.GetTimeType() != Time::TypeFixed) {
			_EXCEPTIONT("Time specified by --time_begin must correspond to a specific date");
		}
	}
	if (strTimeEnd != "") {
		timeEnd = Time(caltype);
		timeEnd.FromFormattedString(strTimeEnd);
		if (timeEnd.GetTimeType() != Time::TypeFixed) {
			_EXCEPTIONT("Time specified by --time_end must correspond to a specific date");
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

	// Parse --timefilter
	if (strTimeFilter == "3hr") {
		strTimeFilter = "(0|3|6|9|12|15|18|21)";
	}
	if (strTimeFilter == "6hr") {
		strTimeFilter = "(0|6|12|18)";
	}
	if (strTimeFilter == "daily") {
		strTimeFilter = "0";
	}
	if (strTimeFilter == "") {
		strTimeFilter = "\\d+";
	}

#ifdef TEMPEST_NOREGEX
	if (strTimeFilter != "") {
		_EXCEPTIONT("Cannot use --timefilter with -DTEMPEST_NOREGEX compiler flag");
	}
#endif
#ifndef TEMPEST_NOREGEX
	std::regex reTimeSubset;
	if (strTimeFilter != "") {
		// Test regex support
		TestRegex();

		try {
			reTimeSubset.assign(strTimeFilter);
		} catch(std::regex_error & reerr) {
			_EXCEPTION2("Parse error in --timefilter regular expression \"%s\" (code %i)",
				strTimeFilter.c_str(), reerr.code());
		}
	}
#endif

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

	// Parse the input into mapCandidateInfo, then store iterators to all elements
	// of mapCandidateInfo in vecIterCandidates to enable sequential access.
	TimeToCandidateInfoMap mapCandidateInfo;
	std::vector<TimeToCandidateInfoMapIterator> vecIterCandidates;
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
				fAllowRepeatedTimes,
				timeBegin,
				timeEnd,
				sGridDims,
				reTimeSubset);
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
	// 	Remove segments 

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
			vecNodes[t][i].lonRad = DegToRad(std::stod(tscinfo[i][iLonIndex]));
			vecNodes[t][i].latRad = DegToRad(std::stod(tscinfo[i][iLatIndex]));

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
	// 	Find all paths
	AnnounceStartBlock("Constructing paths");

	std::vector<SimplePath> vecPaths;

	std::vector<SimplePathSegmentSet> vecPathSegmentsSet;

	AnnounceStartBlock("Populating set of path segments");

	if (strPrioritize == "") {

		// Vector over time of sets of all path segments that start at that time

		GeneratePathSegmentsSetBasic(
			vecIterCandidates,
			vecNodes,
			vecKDTrees,
			nMaxGapSteps,
			dMaxGapSeconds,
			dRangeDeg,
			vecPathSegmentsSet);

	} else {
		int ixPriorityCol = (-1);
		for (int i = 0; i < vecFormatStrings.size(); i++) {
			if (vecFormatStrings[i] == strPrioritize) {
				ixPriorityCol = i;
				break;
			}
		}

		if (ixPriorityCol == (-1)) {
			_EXCEPTION1("Format string (--in_fmt) does not contain priority column \"%s\"", strPrioritize.c_str());
		}

		GeneratePathSegmentsWithPriority(
			vecIterCandidates,
			vecNodes,
			vecKDTrees,
			ixPriorityCol,
			nMaxGapSteps,
			dMaxGapSeconds,
			dRangeDeg,
			vecPathSegmentsSet);
	}

	AnnounceEndBlock("Done");

	AnnounceStartBlock("Connecting path segments");

	GeneratePathsBasic(
		vecIterCandidates,
		vecPathSegmentsSet,
		vecPaths);

	AnnounceEndBlock("Done");

	AnnounceEndBlock("Done");

	/////////////////////////////////////////////
	// 	Filter paths
	AnnounceStartBlock("Filtering paths");

	int nRejectedMinTimePaths = 0;
	int nRejectedMinEndpointDistPaths = 0;
	int nRejectedMinPathDistPaths = 0;
	int nRejectedThresholdPaths = 0;

	std::vector<SimplePath> vecPathsUnfiltered = vecPaths;
	vecPaths.clear();

	// Loop through all paths
	for (const SimplePath & path : vecPathsUnfiltered) {

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

		// Additional columns for zonal and meridional velocity of track
		if (fAddVelocity) {
			nodefile.m_cdh.push_back("uvel");
			nodefile.m_cdh.push_back("vvel");
		}

		// Copy over column data from candidate file to track file for output
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
					pathnode.PushColumnDataString(
						vecIterCandidates[iTime]->second[iCandidate][j]);
				}
			}

			// Add velocity of the system in m/s
			if (fAddVelocity) {

				char szUvel[32];
				char szVvel[32];

				if (vecPaths[i].m_iTimes.size() == 1) {
					snprintf(szUvel, 32, "%3.6f", 0.0);
					snprintf(szVvel, 32, "%3.6f", 0.0);
				}

				for (int t = 0; t < path.size()-1; t++) {

					PathNode & pathnode_curr = path[t];
					PathNode & pathnode_next = path[t+1];

					// Get longitude and latitude of current and next points.
					double dLon0Deg = pathnode_curr.GetColumnDataAsDouble(iLonIndex);
					double dLat0Deg = pathnode_curr.GetColumnDataAsDouble(iLatIndex);

					double dLon1Deg = pathnode_next.GetColumnDataAsDouble(iLonIndex);
					double dLat1Deg = pathnode_next.GetColumnDataAsDouble(iLatIndex);

					double dDeltaSeconds = pathnode_next.m_time - pathnode_curr.m_time;

					// Calculate great circle direction using stereographic projection
					// At poles use a consistent stereographic plane.
					double dUvelRad;
					double dVvelRad;

					GreatCircleDirection_Rad(
						DegToRad(dLon0Deg),
						DegToRad(dLat0Deg),
						DegToRad(dLon1Deg),
						DegToRad(dLat1Deg),
						dUvelRad,
						dVvelRad);

					// Calculate velocities
					snprintf(szUvel, 32, "%3.6e", EarthRadius * dUvelRad / dDeltaSeconds);
					snprintf(szVvel, 32, "%3.6e", EarthRadius * dVvelRad / dDeltaSeconds);

					pathnode_curr.PushColumnDataString(szUvel);
					pathnode_curr.PushColumnDataString(szVvel);
				}

				PathNode & pathnode_last = path[path.size()-1];

				pathnode_last.PushColumnDataString(szUvel);
				pathnode_last.PushColumnDataString(szVvel);
			}

			path.m_timeStart = path[0].m_time;
		}

		if (strOutputFileFormat == "gfdl") {
			nodefile.Write(strOutputFile, NULL, NULL, NodeFile::FileFormatGFDL, false, fOutputSeconds);
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


