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
			(strVariables[iVarCurrent] == '\n')
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
	std::vector< std::vector< std::vector<std::string> > > & vecCandidates
) {
	// Open file for reading
	FILE * fp = fopen(strInputFile.c_str(), "r");

	// Buffer storage
	std::string strLine;
	char szBuffer[1024];

	// Current read state
	enum ReadState {
		ReadState_Time,
		ReadState_Candidate
	} eReadState = ReadState_Time;

	// Number of entries per candidate
	int nFormatEntries = vecFormatStrings.size();

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
			vecTimes.resize(iTime + 1);

			ParseVariableList(strLine, vecTimes[iTime]);

			if (vecTimes[iTime].size() != 5) {
				_EXCEPTION1("Malformed time string:\n%s", szBuffer);
			}

			iCandidate = 0;
			nCandidates = atoi(vecTimes[iTime][3].c_str());

			vecCandidates.resize(iTime + 1);
			vecCandidates[iTime].resize(nCandidates);

			eReadState = ReadState_Candidate;

		// Parse candidate information
		} else if (eReadState == ReadState_Candidate) {
			ParseVariableList(strLine, vecCandidates[iTime][iCandidate]);

			iCandidate++;
			if (iCandidate == nCandidates) {
				eReadState = ReadState_Time;
				iTime++;
				iCandidate = 0;
			}
		}
	}

	fclose(fp);	
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

	// Maximum time gap (in time steps)
	int nMaxGapSize;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strFormat, "format", "no,i,j,lon,lat");
		CommandLineDouble(dRange, "range", 5.0);
		CommandLineInt(nMinPathLength, "minlength", 3);
		CommandLineInt(nMaxGapSize, "maxgap", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

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

	// Parse the input
	std::vector< std::vector<std::string> > vecTimes;
	std::vector< std::vector< std::vector<std::string> > > vecCandidates;

	{
		AnnounceStartBlock("Loading candidate data");

		ParseInput(
			strInputFile,
			vecFormatStrings,
			vecTimes,
			vecCandidates);

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

			double dX = vecNodes[t][i].x;
			double dY = vecNodes[t][i].y;
			double dZ = vecNodes[t][i].z;

			double dLat = vecNodes[t][i].lat;
			double dLon = vecNodes[t][i].lon;

			kdres * set = kd_nearest3(vecKDTrees[t+1], dX, dY, dZ);

			int iRes =
				  reinterpret_cast<int*>(kd_res_item_data(set))
				- reinterpret_cast<int*>(noptr);

			kd_res_free(set);

			// Great circle distance between points
			double dLonC = vecNodes[t+1][iRes].lon;
			double dLatC = vecNodes[t+1][iRes].lat;

			double dR = 180.0 / M_PI * acos(sin(dLatC) * sin(dLat)
					+ cos(dLatC) * cos(dLat) * cos(dLon - dLonC));

			// Verify great circle distance satisfies range requirement
			if (dR <= dRange) {

				// Insert new path segment into set of path segments
				vecPathSegmentsSet[t].insert(
					PathSegment(t, i, t+1, iRes));
			}
		}
	}

	AnnounceEndBlock("Done");

	// Work forwards to find all paths
	AnnounceStartBlock("Constructing paths");

	std::vector< Path > vecPaths;

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

				vecPathSegmentsSet[tx].erase(iterSeg);

				if (tx == vecTimes.size()-1) {
					break;
				}

				PathSegment segFind(
					iterSeg->m_iTime[1], iterSeg->m_iCandidate[1], 0, 0);

				iterSeg = vecPathSegmentsSet[txnext].find(segFind);

				if (iterSeg == vecPathSegmentsSet[txnext].end()) {
					break;
				}

				tx = txnext;
			}

			if (path.m_iTimes.size() > nMinPathLength) {
/*
				printf("Start Time: %i\n", path.m_iTimes[0]);
				for (int i = 0; i < path.m_iTimes.size(); i++) {
					printf("%i, ", path.m_iCandidates[i]);
				}
*/
				vecPaths.push_back(path);
			}
		}
	}

	Announce("Total paths found: %i", vecPaths.size());
	AnnounceEndBlock("Done");

	// Write results out
	AnnounceStartBlock("Writing results");
	{
		FILE * fp = fopen(strOutputFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Failed to open output file \"%s\"",
				strOutputFile.c_str());
		}

		for (int i = 0; i < vecPaths.size(); i++) {
			int iStartTime = vecPaths[i].m_iTimes[0];

			fprintf(fp, "start\t");
			for (int j = 0; j < vecTimes[iStartTime].size(); j++) {
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
					fprintf(fp, "%s\t", vecTimes[iTime][j].c_str());
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
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

