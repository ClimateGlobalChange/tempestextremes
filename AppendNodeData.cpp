///////////////////////////////////////////////////////////////////////////////
///
///	\file    AppendNodeData.cpp
///	\author  Paul Ullrich
///	\version February 20, 2015
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

#include "DataVector.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <queue>
#include <set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load in the contents of a text file containing one filename per
///		line and store in a vector of strings.
///	</summary>
void GetFileList(
	const std::string & strFileList,
	std::vector<std::string> & vecFiles
) {
	FILE * fp = fopen(strFileList.c_str(), "r");

	char szBuffer[1024];
	for (;;) {
		fgets(szBuffer, 1024, fp);

		if (feof(fp)) {
			break;
		}

		// Remove end-of-line characters
		for (;;) {
			int nLen = strlen(szBuffer);
			if ((szBuffer[nLen-1] == '\n') ||
				(szBuffer[nLen-1] == '\r') ||
				(szBuffer[nLen-1] == ' ')
			) {
				szBuffer[nLen-1] = '\0';
				continue;
			}
			break;
		}

		vecFiles.push_back(szBuffer);
	}

	if (vecFiles.size() == 0) {
		_EXCEPTION1("No files found in file \"%s\"", strFileList.c_str());
	}

	fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get a DataVector containing the time variable across a list of
///		input files.
///	</summary>
void GetTimeIndices(
	const std::vector<std::string> & vecInputFiles,
	std::vector<int> & vecTimes
) {
	vecTimes.resize(vecInputFiles.size() + 1);

	int nTime = 0;

	for (int f = 0; f < vecInputFiles.size(); f++) {
		NcFile ncFile(vecInputFiles[f].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		NcDim * dimTime = ncFile.get_dim("time");
		if (dimTime == NULL) {
			_EXCEPTIONT("No dimension \"time\" found in input file");
		}

		vecTimes[f] = nTime;

		nTime += static_cast<int>(dimTime->size());
	}
	vecTimes[vecInputFiles.size()] = nTime;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
void FindLocalAverage(
	const DataMatrix<float> & data,
	const DataVector<double> & dataLat,
	const DataVector<double> & dataLon,
	int iLat,
	int iLon,
	double dMaxDist,
	float & dAverage,
	float & dMaxValue
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Number of latitudes and longitudes
	const int nLat = dataLat.GetRows();
	const int nLon = dataLon.GetRows();

	// Queue of nodes that remain to be visited
	std::queue< std::pair<int, int> > queueNodes;
	queueNodes.push( std::pair<int, int>(iLat, iLon) );

	// Set of nodes that have already been visited
	std::set< std::pair<int, int> > setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = dataLat[iLat];
	double dLon0 = dataLon[iLon];

	// Sum
	float dSum = 0.0f;
	float dArea = 0.0f;

	// Reset average
	dAverage = 0.0f;

	// Reset maximum value
	dMaxValue = data[iLat][iLon];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		std::pair<int, int> pr = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(pr) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(pr);

		double dLatThis = dataLat[pr.first];
		double dLonThis = dataLon[pr.second];

		// Great circle distance to this element
		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if (dR > dMaxDist) {
			continue;
		}

		// Add value to sum
		float dLocalArea = cos(dLatThis);

		dArea += dLocalArea;
		dSum += data[pr.first][pr.second] * dLocalArea;

		if (data[pr.first][pr.second] > dMaxValue) {
			dMaxValue = data[pr.first][pr.second];
		}

		// Add all neighbors of this point
		std::pair<int,int> prWest(pr.first, (pr.second + nLon - 1) % nLon);
		if (setNodesVisited.find(prWest) == setNodesVisited.end()) {
			queueNodes.push(prWest);
		}

		std::pair<int,int> prEast(pr.first, (pr.second + 1) % nLon);
		if (setNodesVisited.find(prEast) == setNodesVisited.end()) {
			queueNodes.push(prEast);
		}

		std::pair<int,int> prNorth(pr.first + 1, pr.second);
		if ((prNorth.first < nLat) &&
			(setNodesVisited.find(prNorth) == setNodesVisited.end())
		) {
			queueNodes.push(prNorth);
		}

		std::pair<int,int> prSouth(pr.first - 1, pr.second);
		if ((prSouth.first >= 0) &&
			(setNodesVisited.find(prSouth) == setNodesVisited.end())
		) {
			queueNodes.push(prSouth);
		}
	}

	// Set average
	dAverage = dSum / dArea;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

	// All data files
	std::vector<NcFile *> vecDataNcFiles;

try {

	// Input file
	std::string strInputFile;

	// Data file
	std::string strDataFile;

	// Data file list
	std::string strDataFileList;

	// Maximum distance from PSL min to perform computation (in degrees) 
	double dMaxDist;

	// Input file format
	std::string strInputFormat;

	// Output file (NetCDF)
	std::string strOutputFile;

	// Time column
	int iTimeIxCol;

	// Column in which the longitude index appears
	int iLonIxCol;

	// Column in which the latitude index appears
	int iLatIxCol;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strDataFile, "data", "");
		CommandLineString(strDataFileList, "datalist", "");
		CommandLineDoubleD(dMaxDist, "maxdist", 5.0, "(degrees)");
		CommandLineStringD(strInputFormat, "in_format", "visit", "(std|visit)");
		CommandLineString(strOutputFile, "out", "");
		CommandLineInt(iTimeIxCol, "itimecol", 2);
		CommandLineInt(iLonIxCol, "iloncol", 8);
		CommandLineInt(iLatIxCol, "ilatcol", 9);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if (strInputFile == "") {
		_EXCEPTIONT("No input file (--in) specified");
	}
	if ((strDataFile == "") && (strDataFileList == "")) {
		_EXCEPTIONT("No data files (--data) or (--datalist) specified");
	}
	if ((strDataFile != "") && (strDataFileList != "")) {
		_EXCEPTIONT("Only one (--data) or (--datalist) allowed");
	}
	if (strInputFormat != "visit") {
		_EXCEPTIONT("UNIMPLEMENTED:  Only \"--in_format visit\" supported");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}
	
	// Data file list
	std::vector<std::string> vecDataFiles;

	if (strDataFile != "") {
		vecDataFiles.push_back(strDataFile);
	}
	if (strDataFileList != "") {
		GetFileList(strDataFileList, vecDataFiles);
	}

	// Open all data files
	int nLat;
	int nLon;

	DataVector<double> dataLat;
	DataVector<double> dataLon;

	for (int f = 0; f < vecDataFiles.size(); f++) {
		vecDataNcFiles.push_back(
			new NcFile(vecDataFiles[f].c_str(), NcFile::ReadOnly));

		if (!vecDataNcFiles[f]->is_valid()) {
			_EXCEPTION1("Unable to open data file \"%s\"",
				vecDataFiles[f].c_str());
		}

		// Determine number of longitudes and latitudes
		NcDim * dimLat = vecDataNcFiles[f]->get_dim("lat");
		NcDim * dimLon = vecDataNcFiles[f]->get_dim("lon");

		if (dimLat == NULL) {
			_EXCEPTION1("File \"%s\" does not contain dimension \"lat\"",
				vecDataFiles[f].c_str());
		}
		if (dimLon == NULL) {
			_EXCEPTION1("File \"%s\" does not contain dimension \"lon\"",
				vecDataFiles[f].c_str());
		}

		// Initialize latitude / longitude data
		if (f == 0) {
			nLat = dimLat->size();
			nLon = dimLon->size();

			NcVar * varLat = vecDataNcFiles[f]->get_var("lat");
			NcVar * varLon = vecDataNcFiles[f]->get_var("lon");

			if (varLat == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"lat\"",
					vecDataFiles[f].c_str());
			}
			if (varLon == NULL) {
				_EXCEPTION1("File \"%s\" does not contain variable \"lon\"",
					vecDataFiles[f].c_str());
			}

			dataLat.Initialize(nLat);
			dataLon.Initialize(nLon);

			varLat->get(dataLat, nLat);
			varLon->get(dataLon, nLon);

			for (int j = 0; j < nLat; j++) {
				dataLat[j] *= M_PI / 180.0;
			}

			for (int i = 0; i < nLon; i++) {
				dataLon[i] *= M_PI / 180.0;
			}

		} else {
			if (nLat != dimLat->size()) {
				_EXCEPTION1("File \"%s\" latitude dimension size mismatch",
					vecDataFiles[f].c_str());
			}
			if (nLon != dimLon->size()) {
				_EXCEPTION1("File \"%s\" longitude dimension size mismatch",
					vecDataFiles[f].c_str());
			}
		}
	}

	// Get time indices
	std::vector<int> vecTimes;

	GetTimeIndices(vecDataFiles, vecTimes);

	int nTime = vecTimes[vecTimes.size()-1];

	// String buffer
	char szBuffer[1024];
	char szSecondBuffer[1024];

	// Open the input file
	FILE * fp = fopen(strInputFile.c_str(), "r");
	if (fp == NULL) {
		_EXCEPTION1("Unable to open input file \"%s\"",
			strInputFile.c_str());
	}

	// Open the output file
	FILE * fpout = fopen(strOutputFile.c_str(), "w");
	if (fpout == NULL) {
		_EXCEPTION1("Unable to open output file \"%s\"",
			strOutputFile.c_str());
	}

	// PRECT data matrix
	DataMatrix<float> dPRECT(nLat, nLon);

	// Loop through all lines of input file
	for (;;) {

		// Read in the next line
		fgets(szBuffer, 1024, fp);

		int nLength = strlen(szBuffer);

		// Check for end of file
		if (feof(fp)) {
			break;
		}

		// Check for comment line
		if (szBuffer[0] == '#') {
			continue;
		}
		if (nLength == 0) {
			continue;
		}

		strcpy(szSecondBuffer, szBuffer);
		szSecondBuffer[nLength-1] = '\0';

		// Parse line
		int iTime = 0;
		int iLon = 0;
		int iLat = 0;

		int iCol = 0;
		int iLast = 0;

		bool fWhitespace = true;

		for (int i = 0; i <= nLength; i++) {
			if ((szBuffer[i] == ' ') ||
				(szBuffer[i] == ',') ||
				(szBuffer[i] == '\t') ||
				(szBuffer[i] == '\0')
			) {
				if (!fWhitespace) {
					if (iCol == iTimeIxCol) {
						szBuffer[i] = '\0';
						iTime = atoi(&(szBuffer[iLast]));
					}
					if (iCol == iLonIxCol) {
						szBuffer[i] = '\0';
						iLon = atoi(&(szBuffer[iLast]));
					}
					if (iCol == iLatIxCol) {
						szBuffer[i] = '\0';
						iLat = atoi(&(szBuffer[iLast]));
					}
				}

				fWhitespace = true;

			} else {
				if (fWhitespace) {
					iLast = i;
					iCol++;
				}
				fWhitespace = false;
			}
		}

		if ((iLat < 0) || (iLat >= nLat)) {
			_EXCEPTION1("Latitude index (%i) out of range", iLat);
		}
		if ((iLon < 0) || (iLon >= nLon)) {
			_EXCEPTION1("Longitude index (%i) out of range", iLon);
		}
		if ((iTime < 0) || (iTime >= nTime)) {
			_EXCEPTION1("Time index (%i) out of range", iTime);
		}

		// Find the correct file
		int iFile = (-1);
		for (int f = 0; f < vecTimes.size()-1; f++) {
			if ((iTime >= vecTimes[f]) && (iTime <= vecTimes[f+1])) {
				iFile = f;
			}
		}
		if (iFile == (-1)) {
			_EXCEPTION1("Time index (%i) out of range", iTime);
		}

		// Load in PRECT from file
		NcVar * varPRECT = vecDataNcFiles[iFile]->get_var("PRECT");
		if (varPRECT == NULL) {
			_EXCEPTION1("File \"%s\" does not contain variable \"PRECT\"",
				vecDataFiles[iFile].c_str());
		}
		
		varPRECT->set_cur(iTime - vecTimes[iFile], 0, 0);
		varPRECT->get(&(dPRECT[0][0]), 1, nLat, nLon);

		// Compute average of PRECT from file
		float dAverage;
		float dMaxValue;

		FindLocalAverage(
			dPRECT,
			dataLat,
			dataLon,
			iLat,
			iLon,
			dMaxDist,
			dAverage,
			dMaxValue);

		// Write to file
		fprintf(fpout, "%s,\t%1.5e,\t%1.5e\n", szSecondBuffer, dAverage, dMaxValue);
	}

	fclose(fp);
	fclose(fpout);

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

	// Close all open files
	for (int f = 0; f < vecDataNcFiles.size(); f++) {
		delete vecDataNcFiles[f];
	}
}

///////////////////////////////////////////////////////////////////////////////

