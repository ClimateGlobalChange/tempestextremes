///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculatePosthocOutput.cpp
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
#include "Variable.h"
#include "AutoCurator.h"

#include "netcdfcpp.h"

#include <fstream>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct VelocityVector {
	double dU;
	double dV;
};

///////////////////////////////////////////////////////////////////////////////

struct PathNode {
	std::vector<int> m_coord;
	std::vector<std::string> m_vecData;
	std::vector<VelocityVector> m_dVelocityLat;
	std::vector<double> m_dVelocityLon;
};

///////////////////////////////////////////////////////////////////////////////

struct Path {
	Time m_timeStart;
	std::vector<PathNode> m_vecPathNodes;
};

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<Path> PathVector;

///////////////////////////////////////////////////////////////////////////////

void CalculateRadialWindProfile(
) {
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {

	// Input text file
	std::string strInputFile;

	// Input list of text files
	std::string strInputFileList;

	// Input file type
	std::string strInputFileType;

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Connectivity file
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Output file
	std::string strOutputFile;

	// Append output to input file
	bool fOutputAppend;

	// Variables for calculating radial wind profile
	std::string strRadialWindProfileVars;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_file", "");
		//CommandLineString(strInputFileList, "in_file_list", "");
		CommandLineStringD(strInputFileType, "in_file_type", "SN", "[DCU|SN]");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strOutputFile, "out_file", "");
		CommandLineBool(fOutputAppend, "out_append");

		CommandLineStringD(strRadialWindProfileVars, "radial_wind_profile", "", "U,V,bins,radius");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Create autocurator
	AutoCurator autocurator;

	// Check arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_file) or (--in_file_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_file) or (--in_file_list)"
			" may be specified");
	}

	// Input file type
	enum InputFileType {
		InputFileTypeDCU,
		InputFileTypeSN
	} iftype;
	if (strInputFileType == "DCU") {
		iftype = InputFileTypeDCU;
	} else if (strInputFileType == "SN") {
		iftype = InputFileTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_file_type, expected \"SN\" or \"DCU\"");
	}

	// Radial wind profile can only be calculated with StitchNodes output
	if ((strRadialWindProfileVars != "") && (iftype != InputFileTypeSN)) {
		_EXCEPTIONT("--radial_wind_profile can only be used with --in_file_type SN");
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Curate input data
	if (strInputData.length() != 0) {
		autocurator.IndexFiles(strInputData);

	} else {
		std::ifstream ifInputDataList(strInputDataList.c_str());
		if (!ifInputDataList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputDataList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			autocurator.IndexFiles(strFileLine);
		}
	}

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

	// No connectivity file; check for latitude/longitude dimension
	} else {
		const std::vector<std::string> & vecFiles = autocurator.GetFilenames();

		if (vecFiles.size() < 1) {
			_EXCEPTIONT("No data files specified");
		}

		NcFile ncFile(vecFiles[0].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", vecFiles[0].c_str());
		}

		grid.GenerateLatitudeLongitude(&ncFile, fRegional);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputFiles.push_back(strFileLine);
		}
	}

	// Loop through all files
	std::string strBuffer;

	std::vector<int> coord;
	coord.resize(grid.m_nGridDim.size());

	// Vector of Path information loaded from StitchNodes file
	PathVector vecPaths;

	// A map from Times to file lines, used for StitchNodes formatted output
	std::map<Time, std::vector<int> > mapTimeToFileLines;

	for (int f = 0; f < vecInputFiles.size(); f++) {

		std::ifstream ifInput(vecInputFiles[f]);
		if (!ifInput.is_open()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[f].c_str());
		}

		int iLine = 1;
		for (;;) {

			int nCount = 0;
			Time time;

			// Read header lines
			{
				getline(ifInput, strBuffer);
				if (ifInput.eof()) {
					break;
				}

				std::istringstream iss(strBuffer);

				// DetectCyclonesUnstructured output
				if (iftype == InputFileTypeDCU) {
					int iYear;
					int iMonth;
					int iDay;
					int iHour;

					iss >> iYear;
					iss >> iMonth;
					iss >> iDay;
					iss >> nCount;
					iss >> iHour;

					if (iss.eof()) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					time = Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						autocurator.GetCalendarType());

				// StitchNodes output
				} else if (iftype == InputFileTypeSN) {
					std::string strStart;
					int iYear;
					int iMonth;
					int iDay;
					int iHour;

					iss >> strStart;
					iss >> nCount;
					iss >> iYear;
					iss >> iMonth;
					iss >> iDay;
					iss >> iHour;

					if (iss.eof()) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					vecPaths.resize(vecPaths.size() + 1);
					vecPaths[vecPaths.size()-1].m_timeStart =
						Time(
							iYear,
							iMonth-1,
							iDay-1,
							3600 * iHour,
							0,
							autocurator.GetCalendarType());

					vecPaths[vecPaths.size()-1].m_vecPathNodes.resize(nCount);
				}

				iLine++;
			}

			// Read contents under each header line
			for (int i = 0; i < nCount; i++) {

				getline(ifInput, strBuffer);
				if (ifInput.eof()) {
					break;
				}

				std::istringstream iss(strBuffer);

				for (int n = 0; n < grid.m_nGridDim.size(); n++) {
					iss >> coord[n];
					if ((coord[n] < 0) || (coord[n] >= grid.m_nGridDim[n])) {
						_EXCEPTION1("Coordinate index out of range on line %i", coord[n]);
					}
				}
				if (iss.eof()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, vecInputFiles[f].c_str());
				}

				std::string strBuf;
				std::vector<std::string> vecDelimitedOutput;
				for (;;) {
					iss >> strBuf;
					if (iss.eof()) {
						break;
					}
					vecDelimitedOutput.push_back(strBuf);
				}

				int nOutputSize = vecDelimitedOutput.size();

				// StitchNodes requires us to reorder the data 
				if (iftype == InputFileTypeSN) {
					PathNode & pathnode =
						vecPaths[vecPaths.size()-1].m_vecPathNodes[i];

					if (nOutputSize < 4) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					int iYear = std::stoi(vecDelimitedOutput[nOutputSize-4]);
					int iMonth = std::stoi(vecDelimitedOutput[nOutputSize-3]);
					int iDay = std::stoi(vecDelimitedOutput[nOutputSize-2]);
					int iHour = std::stoi(vecDelimitedOutput[nOutputSize-1]);

					time = Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						autocurator.GetCalendarType());

					std::map<Time, std::vector<int> >::iterator iter =
						mapTimeToFileLines.find(time);
					if (iter == mapTimeToFileLines.end()) {
						iter = mapTimeToFileLines.insert(
							std::map<Time, std::vector<int> >::value_type(
								time, std::vector<int>())).first;
					}
					iter->second.push_back(iLine);

					pathnode.m_coord = coord;
					pathnode.m_vecData = vecDelimitedOutput;
				}

				iLine++;
			}

			// Calculate velocity at each point
		}
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

