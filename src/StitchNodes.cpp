///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchNodes.cpp
///	\author  Paul Ullrich
///	\version October 1st, 2014
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

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"

#include "kdtree.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////

struct Time {
	int nDay;
	int nMonth;
	int nYear;
	int nCount;
	int nHour;
};

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

void ParseInput(
	const std::string & strInputFile,
	const std::vector< std::string > & vecFormatStrings,
	std::vector< std::vector<std::string> > & vecTimes,
	std::vector< std::vector< std::vector<std::string> > > & vecCandidates,
	int nTimeStride = 1
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

	int iAllTime = 0;
	int iTime = 0;

	int iCandidate = 0;
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
			vecTimes.resize(iTime + 1);

			ParseVariableList(strLine, vecTimes[iTime]);

			if (vecTimes[iTime].size() != 5) {
				_EXCEPTION1("Malformed time string:\n%s", strLine.c_str());
			}

			iCandidate = 0;
			nCandidates = atoi(vecTimes[iTime][3].c_str());

			// Verify that this time is on stride
			if (iAllTime % nTimeStride != 0) {
				if (nCandidates != 0) {
					eReadState = ReadState_Candidate;
				} else {
					iAllTime++;
				}
				vecTimes.resize(iTime);
				continue;
			}

			// Prepare to parse candidate data
			vecCandidates.resize(iTime + 1);
			vecCandidates[iTime].resize(nCandidates);

			if (nCandidates != 0) {
				eReadState = ReadState_Candidate;
			} else {
				iTime++;
				iAllTime++;
			}

		// Parse candidate information
		} else if (eReadState == ReadState_Candidate) {

			// Ignore candidates that are not on stride
			if (iAllTime % nTimeStride != 0) {
				iCandidate++;
				if (iCandidate == nCandidates) {
					eReadState = ReadState_Time;
					iAllTime++;
					iCandidate = 0;
				}
				continue;
			}

			// Parse candidates
			ParseVariableList(strLine, vecCandidates[iTime][iCandidate]);

			if (vecCandidates[iTime][iCandidate].size() != nFormatEntries) {
				fWarnInsufficientCandidateInfo = true;
			}

			iCandidate++;
			if (iCandidate == nCandidates) {
				eReadState = ReadState_Time;
				iTime++;
				iAllTime++;
				iCandidate = 0;
			}
		}
	}

	fclose(fp);	

	// Insufficient candidate information
	if (fWarnInsufficientCandidateInfo) {
		Announce("WARNING: One or more candidates do not have match"
				" --format entries");
	}
}

///////////////////////////////////////////////////////////////////////////////

struct Node {
	double x;
	double y;
	double z;

	double lat;
	double lon;
};

///////////////////////////////////////////////////////////////////////////////

typedef std::vector< std::vector<std::string> > TimesVector;

typedef std::vector< std::vector< std::vector<std::string> > > CandidateVector;

///////////////////////////////////////////////////////////////////////////////

class PathSegment {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PathSegment(
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
	bool operator< (const PathSegment & seg) const {
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

typedef std::set<PathSegment> PathSegmentSet;

///////////////////////////////////////////////////////////////////////////////

class Path {

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

typedef std::vector<Path> PathVector;

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
		const Path & path,
		const CandidateVector & vecCandidates
	) {
		int nCount = 0;
		for (int s = 0; s < path.m_iTimes.size(); s++) {
			int t = path.m_iTimes[s];
			int i = path.m_iCandidates[s];

			double dCandidateValue =
				atof(vecCandidates[t][i][m_iColumn].c_str());

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

	// Output file
	std::string strOutputFile;

	// Format string
	std::string strFormat;

	// Range (in degrees)
	double dRange;

	// Minimum path length
	int nMinPathLength;

	// Minimum distance between endpoints of path
	double dMinEndpointDistance;

	// Minimum path length
	double dMinPathDistance;

	// Maximum time gap (in time steps)
	int nMaxGapSize;

	// Time step stride
	int nTimeStride;

	// Output format
	std::string strOutputFormat;

	// Thresholds
	std::string strThreshold;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strFormat, "format", "no,i,j,lon,lat");
		CommandLineDoubleD(dRange, "range", 5.0, "(degrees)");
		CommandLineInt(nMinPathLength, "minlength", 3);
		CommandLineDoubleD(dMinEndpointDistance, "min_endpoint_dist", 0.0, "(degrees)");
		CommandLineDoubleD(dMinPathDistance, "min_path_dist", 0.0, "(degrees)");
		CommandLineInt(nMaxGapSize, "maxgap", 0);
		CommandLineStringD(strThreshold, "threshold", "",
			"[col,op,value,count;...]");
		CommandLineInt(nTimeStride, "timestride", 1);
		CommandLineStringD(strOutputFormat, "out_format", "std", "(std|visit)");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if (strInputFile.size() == 0) {
		_EXCEPTIONT("No input file specified");
	}

	// Output format
	if ((strOutputFormat != "std") &&
		(strOutputFormat != "visit")
	) {
		_EXCEPTIONT("Output format must be either \"std\" or \"visit\"");
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
	TimesVector vecTimes;
	CandidateVector vecCandidates;

	{
		AnnounceStartBlock("Loading candidate data");

		ParseInput(
			strInputFile,
			vecFormatStrings,
			vecTimes,
			vecCandidates,
			nTimeStride);

		Announce("Discrete times: %i", vecTimes.size());

		AnnounceEndBlock("Done");
	}

	// Create kdtree at each time
	AnnounceStartBlock("Creating KD trees at each time level");

	// Null pointer
	int * noptr = NULL;

	// Vector of lat/lon values
	std::vector< std::vector<Node> > vecNodes;
	vecNodes.resize(vecTimes.size());

	// Vector of KD trees
	std::vector<kdtree *> vecKDTrees;
	vecKDTrees.resize(vecTimes.size());

	for (int t = 0; t < vecTimes.size(); t++) {

		// Create a new kdtree
		if (vecCandidates[t].size() == 0) {
			vecKDTrees[t] = NULL;
			continue;
		}

		vecKDTrees[t] = kd_create(3);

		vecNodes[t].resize(vecCandidates[t].size());

		// Insert all points at this time level
		for (int i = 0; i < vecCandidates[t].size(); i++) {
			double dLat = atof(vecCandidates[t][i][iLatIndex].c_str());
			double dLon = atof(vecCandidates[t][i][iLonIndex].c_str());

			dLat *= M_PI / 180.0;
			dLon *= M_PI / 180.0;

			double dX = sin(dLon) * cos(dLat);
			double dY = cos(dLon) * cos(dLat);
			double dZ = sin(dLat);

			vecNodes[t][i].lat = dLat;
			vecNodes[t][i].lon = dLon;

			vecNodes[t][i].x = dX;
			vecNodes[t][i].y = dY;
			vecNodes[t][i].z = dZ;

			kd_insert3(vecKDTrees[t], dX, dY, dZ, reinterpret_cast<void*>(noptr+i));
		}
	}

	AnnounceEndBlock("Done");

	// Create set of path segments
	AnnounceStartBlock("Populating set of path segments");

	std::vector<PathSegmentSet> vecPathSegmentsSet;
	vecPathSegmentsSet.resize(vecTimes.size()-1);

	// Insert nodes from this time level
	for (int t = 0; t < vecTimes.size()-1; t++) {

		// Loop through all points at the current time level
		for (int i = 0; i < vecCandidates[t].size(); i++) {

			if (vecCandidates[t+1].size() == 0) {
				break;
			}

			double dX = vecNodes[t][i].x;
			double dY = vecNodes[t][i].y;
			double dZ = vecNodes[t][i].z;

			double dLat = vecNodes[t][i].lat;
			double dLon = vecNodes[t][i].lon;

			for (int g = 1; g <= nMaxGapSize+1; g++) {
				if (t+g >= vecTimes.size()) {
					break;
				}

				if (vecKDTrees[t+g] == NULL) {
					continue;
				}

				kdres * set = kd_nearest3(vecKDTrees[t+g], dX, dY, dZ);

				if (kd_res_size(set) == 0) {
					kd_res_free(set);
					break;
				}

				int iRes =
					  reinterpret_cast<int*>(kd_res_item_data(set))
					- reinterpret_cast<int*>(noptr);

				kd_res_free(set);

				// Great circle distance between points
				double dLonC = vecNodes[t+g][iRes].lon;
				double dLatC = vecNodes[t+g][iRes].lat;

				double dR = 180.0 / M_PI * acos(sin(dLatC) * sin(dLat)
						+ cos(dLatC) * cos(dLat) * cos(dLon - dLonC));

				// Verify great circle distance satisfies range requirement
				if (dR <= dRange) {

					// Insert new path segment into set of path segments
					vecPathSegmentsSet[t].insert(
						PathSegment(t, i, t+g, iRes));

					break;
				}
			}
		}
	}

	AnnounceEndBlock("Done");

	// Work forwards to find all paths
	AnnounceStartBlock("Constructing paths");

	std::vector< Path > vecPaths;

	int nRejectedMinLengthPaths = 0;
	int nRejectedMinEndpointDistPaths = 0;
	int nRejectedMinPathDistPaths = 0;
	int nRejectedThresholdPaths = 0;

	// Loop through all times
	for (int t = 0; t < vecTimes.size()-1; t++) {

		// Loop through all remaining segments
		while (vecPathSegmentsSet[t].size() > 0) {

			// Create a new path
			Path path;

			PathSegmentSet::iterator iterSeg
				= vecPathSegmentsSet[t].begin();

			path.m_iTimes.push_back(iterSeg->m_iTime[0]);
			path.m_iCandidates.push_back(iterSeg->m_iCandidate[0]);

			int tx = t;

			for (;;) {
				path.m_iTimes.push_back(iterSeg->m_iTime[1]);
				path.m_iCandidates.push_back(iterSeg->m_iCandidate[1]);

				int txnext = iterSeg->m_iTime[1];

				if (txnext >= vecTimes.size()-1) {
					vecPathSegmentsSet[tx].erase(iterSeg);
					break;
				}

				PathSegment segFind(
					iterSeg->m_iTime[1], iterSeg->m_iCandidate[1], 0, 0);

				vecPathSegmentsSet[tx].erase(iterSeg);

				iterSeg = vecPathSegmentsSet[txnext].find(segFind);

				if (iterSeg == vecPathSegmentsSet[txnext].end()) {
					break;
				}

				tx = txnext;
			}

			// Reject path due to minimum length
			if (path.m_iTimes.size() < nMinPathLength) {
				nRejectedMinLengthPaths++;
				continue;
			}

			// Reject path due to minimum endpoint distance
			if (dMinEndpointDistance > 0.0) {
				int nT = path.m_iTimes.size();

				int iTime0 = path.m_iTimes[0];
				int iRes0  = path.m_iCandidates[0];

				int iTime1 = path.m_iTimes[nT-1];
				int iRes1  = path.m_iCandidates[nT-1];

				double dLon0 = vecNodes[iTime0][iRes0].lon;
				double dLat0 = vecNodes[iTime0][iRes0].lat;

				double dLon1 = vecNodes[iTime1][iRes1].lon;
				double dLat1 = vecNodes[iTime1][iRes1].lat;

				double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLat1)
						+ cos(dLat0) * cos(dLat1) * cos(dLon0 - dLon1));

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

					double dLon0 = vecNodes[iTime0][iRes0].lon;
					double dLat0 = vecNodes[iTime0][iRes0].lat;

					double dLon1 = vecNodes[iTime1][iRes1].lon;
					double dLat1 = vecNodes[iTime1][iRes1].lat;

					double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLat1)
						+ cos(dLat0) * cos(dLat1) * cos(dLon0 - dLon1));

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
						vecCandidates);

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

	Announce("Paths rejected (minlength): %i", nRejectedMinLengthPaths);
	Announce("Paths rejected (minendpointdist): %i", nRejectedMinEndpointDistPaths);
	Announce("Paths rejected (minpathdist): %i", nRejectedMinPathDistPaths);
	Announce("Paths rejected (threshold): %i", nRejectedThresholdPaths);
	Announce("Total paths found: %i", vecPaths.size());
	AnnounceEndBlock("Done");

	// Write results out
	AnnounceStartBlock("Writing results");
	if (strOutputFormat == "std") {
		FILE * fp = fopen(strOutputFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Failed to open output file \"%s\"",
				strOutputFile.c_str());
		}

		for (int i = 0; i < vecPaths.size(); i++) {
			int iStartTime = vecPaths[i].m_iTimes[0];

			fprintf(fp, "start\t");
			fprintf(fp, "%li\t", vecPaths[i].m_iTimes.size());
			for (int j = 0; j < vecTimes[iStartTime].size(); j++) {
				if (j == 3) {
					continue;
				}
				fprintf(fp, "%s\t", vecTimes[iStartTime][j].c_str());
			}
			fprintf(fp, "\n");

			for (int t = 0; t < vecPaths[i].m_iTimes.size(); t++) {
				int iTime = vecPaths[i].m_iTimes[t];
				int iCandidate = vecPaths[i].m_iCandidates[t];

				fprintf(fp, "\t");
				for (int j = 0; j < vecCandidates[iTime][iCandidate].size(); j++) {
					fprintf(fp, "%s\t",
						vecCandidates[iTime][iCandidate][j].c_str());
				}
				for (int j = 0; j < vecTimes[iTime].size(); j++) {
					if (j == 3) {
						continue;
					}
					fprintf(fp, "%s\t", vecTimes[iTime][j].c_str());
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);

	} else if (strOutputFormat == "visit") {
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

