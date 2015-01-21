///////////////////////////////////////////////////////////////////////////////
///
///	\file    DensityNodes.cpp
///	\author  Paul Ullrich
///	\version January 20, 2015
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

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Load in the contents of a text file containing one filename per
///		line and store in a vector of strings.
///	</summary>
void GetInputFileList(
	const std::string & strInputFileList,
	std::vector<std::string> & vecInputFiles
) {
	FILE * fp = fopen(strInputFileList.c_str(), "r");

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

		vecInputFiles.push_back(szBuffer);
	}

	if (vecInputFiles.size() == 0) {
		_EXCEPTION1("No files found in file \"%s\"", strInputFileList.c_str());
	}

	fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {

	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Input file format
	std::string strInputFormat;

	// NetCDF file containing latitude and longitude arrays
	std::string strLatLonFile;

	// Output file (NetCDF)
	std::string strOutputFile;

	// Output variable name
	std::string strOutputVariable;

	// Column in which the longitude index appears
	int iLonIxCol;

	// Column in which the latitude index appears
	int iLatIxCol;

	// Number of latitudes in output
	int nLat;

	// Number of longitudes in output
	int nLon;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "inlist", "");
		CommandLineStringD(strInputFormat, "in_format", "visit", "(std|visit)");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputVariable, "outvar", "density");
		CommandLineInt(iLonIxCol, "iloncol", 8);
		CommandLineInt(iLatIxCol, "ilatcol", 9);
		CommandLineInt(nLat, "nlat", 0);
		CommandLineInt(nLon, "nlon", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if ((strInputFile == "") && (strInputFileList == "")) {
		_EXCEPTIONT("No input file (--in) or (--inlist) specified");
	}
	if ((strInputFile != "") && (strInputFileList != "")) {
		_EXCEPTIONT("Only one input file (--in) or (--inlist) allowed");
	}
	if (strInputFormat != "visit") {
		_EXCEPTIONT("UNIMPLEMENTED:  Only \"--in_format visit\" supported");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check output variable
	if (strOutputVariable == "") {
		_EXCEPTIONT("No output variable name (--outvar) specified");
	}

	// Number of latitudes and longitudes
	if (nLat == 0) {
		_EXCEPTIONT("UNIMPLEMENTED: --nlat must be specified currently");
	}
	if (nLon == 0) {
		_EXCEPTIONT("UNIMPLEMENTED: --nlon must be specified currently");
	}

	// Input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile != "") {
		vecInputFiles.push_back(strInputFile);
	}
	if (strInputFileList != "") {
		GetInputFileList(strInputFileList, vecInputFiles);
	}

	int nFiles = vecInputFiles.size();

	// Density
	DataMatrix<int> nCounts;
	nCounts.Initialize(nLat, nLon);

	// Loop through all files in list
	AnnounceStartBlock("Processing files");

	std::string strBuffer;
	strBuffer.reserve(1024);

	for (int f = 0; f < nFiles; f++) {
		Announce("File \"%s\"", vecInputFiles[f].c_str());

		FILE * fp = fopen(vecInputFiles[f].c_str(), "r");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		for (;;) {

			// Read in the next line
			fgets(&(strBuffer[0]), 1024, fp);

			int nLength = strlen(&(strBuffer[0]));

			// Check for end of file
			if (feof(fp)) {
				break;
			}

			// Check for comment line
			if (strBuffer[0] == '#') {
				continue;
			}

			// Parse line
			int iLon = 0;
			int iLat = 0;

			int iCol = 0;
			int iLast = 0;

			bool fWhitespace = true;

			for (int i = 0; i <= nLength; i++) {
				if ((strBuffer[i] == ' ') ||
					(strBuffer[i] == ',') ||
					(strBuffer[i] == '\t') ||
					(strBuffer[i] == '\0')
				) {
					if (!fWhitespace) {
						if (iCol == iLonIxCol) {
							strBuffer[i] = '\0';
							iLon = atoi(&(strBuffer[iLast]));
						}
						if (iCol == iLatIxCol) {
							strBuffer[i] = '\0';
							iLat = atoi(&(strBuffer[iLast]));
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

			nCounts[iLat][iLon]++;
		}

		fclose(fp);
	}

	AnnounceEndBlock("Done");

	// Output results
	AnnounceStartBlock("Output results");

	// Load the netcdf output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"",
			strOutputFile.c_str());
	}

	// Create output
	NcDim * dimLat = ncOutput.add_dim("lat", nLat);
	NcDim * dimLon = ncOutput.add_dim("lon", nLon);

	NcVar * varCount =
		ncOutput.add_var(
			strOutputVariable.c_str(),
			ncInt,
			dimLat,
			dimLon);

	varCount->put(&(nCounts[0][0]), nLat, nLon);

	ncOutput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

