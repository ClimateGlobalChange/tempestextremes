///////////////////////////////////////////////////////////////////////////////
///
///	\file    SpineARs.cpp
///	\author  Paul Ullrich
///	\version December 1st, 2016
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

#include "Variable.h"
#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <set>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse the list of input files.
///	</summary>
void ParseInputFiles(
	const std::string & strInputFile,
	std::vector<NcFile *> & vecFiles
) {
	int iLast = 0;
	for (int i = 0; i <= strInputFile.length(); i++) {
		if ((i == strInputFile.length()) ||
		    (strInputFile[i] == ';')
		) {
			std::string strFile =
				strInputFile.substr(iLast, i - iLast);

			NcFile * pNewFile = new NcFile(strFile.c_str());

			if (!pNewFile->is_valid()) {
				_EXCEPTION1("Cannot open input file \"%s\"",
					strFile.c_str());
			}

			vecFiles.push_back(pNewFile);
			iLast = i+1;
		}
	}

	if (vecFiles.size() == 0) {
		_EXCEPTION1("No input files found in \"%s\"",
			strInputFile.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

class SpineARsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SpineARsParam() :
		fpLog(NULL),
		ixSearchByVar(0),
		strOutputVariable(""),
		iLaplacianSize(0),
		dMinLaplacian(0.0),
		dMinAbsLat(0.0),
		dEqBandMaxLat(0.0),
		dMinIWV(0.0),
		dMinArea(0.0),
		dMaxArealFraction(0.0),
		dZonalMeanWeight(0.0),
		dZonalMaxWeight(0.0),
		dMeridMeanWeight(0.0),
		dMeridMaxWeight(0.0),
		nMinArea(0),
		nAddTimeDim(0),
		strAddTimeDimUnits(""),
		fOutputLaplacian(false),
		fRegional(false),
		strLongitudeName("lon"),
		strLatitudeName("lat")
	{ }

public:
	// Log
	FILE * fpLog;

	// Variable index to search on
	VariableIndex ixSearchByVar;

	// Output variable name
	std::string strOutputVariable;

	// Size of the Laplacian
	int iLaplacianSize;

	// Minimum Laplacian
	double dMinLaplacian;

	// Minimum absolute latitude
	double dMinAbsLat;

	// Maximum absolute latitude of equatorial band
	double dEqBandMaxLat;

	// Minimum pointwise integrated water vapor
	double dMinIWV;

	// Minimum area (grid cells)
	double dMinArea;

	// Maximal areal fraction
	double dMaxArealFraction;

	// Zonal mean weight
	double dZonalMeanWeight;

	// Zonal max weight
	double dZonalMaxWeight;

	// Meridional mean weight
	double dMeridMeanWeight;

	// Meridional max weight
	double dMeridMaxWeight;

	// Minimum area
	int nMinArea;

	// Add a time dimension if absent
	int nAddTimeDim;

	// Time dimension units
	std::string strAddTimeDimUnits;

	// Output Laplacian
	bool fOutputLaplacian;

	// Regional data
	bool fRegional;

	// Name of longitude variabe
	std::string strLongitudeName;

	// Name of latitude variable
	std::string strLatitudeName;

};

///////////////////////////////////////////////////////////////////////////////

void SpineARs(
	int iFile,
	const std::string & strInputFiles,
	const std::string & strOutputFile,
	VariableRegistry & varreg,
	const SpineARsParam & param
) {

	// Check if zonal weights are specified
	bool fHasZonalWeight = false;
	if ((param.dZonalMeanWeight != 0.0) || (param.dZonalMaxWeight != 0.0)) {
		fHasZonalWeight = true;
	}

	// Check if meridional weights are specified
	bool fHasMeridWeight = false;
	if ((param.dMeridMeanWeight != 0.0) || (param.dMeridMaxWeight != 0.0)) {
		fHasMeridWeight = true;
	}

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Load in the benchmark file
	NcFileVector vecFiles;

	ParseInputFiles(strInputFiles, vecFiles);

	// Define the SimpleGrid
	SimpleGrid grid;
	AnnounceStartBlock("Generating RLL grid data");
	grid.GenerateLatitudeLongitude(
		vecFiles[0],
		param.fRegional,
		param.strLatitudeName,
		param.strLongitudeName);
	AnnounceEndBlock("Done");

	// Get the time dimension
	NcDim * dimTime = vecFiles[0]->get_dim("time");
	//if (dimTime == NULL) {
	//	_EXCEPTIONT("Error accessing dimension \"time\"");
	//}

	// Get the longitude dimension
	NcDim * dimLon = vecFiles[0]->get_dim(param.strLongitudeName.c_str());
	if (dimLon == NULL) {
		_EXCEPTION1("Error accessing dimension \"%s\"", param.strLongitudeName.c_str());
	}

	// Get the longitude variable
	NcVar * varLon = vecFiles[0]->get_var(param.strLongitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"", param.strLongitudeName.c_str());
	}

	DataArray1D<double> dLonDeg(dimLon->size());
	varLon->get(&(dLonDeg[0]), dimLon->size());

	// Get the latitude dimension
	NcDim * dimLat = vecFiles[0]->get_dim(param.strLatitudeName.c_str());
	if (dimLat == NULL) {
		_EXCEPTION1("Error accessing dimension \"%s\"", param.strLatitudeName.c_str());
	}

	// Get the latitude variable
	NcVar * varLat = vecFiles[0]->get_var(param.strLatitudeName.c_str());
	if (varLat == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"", param.strLatitudeName.c_str());
	}

	DataArray1D<double> dLatDeg(dimLat->size());
	varLat->get(&(dLatDeg[0]), dimLat->size());

	// Open the NetCDF output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Copy over latitude, longitude and time variables to output file
	NcDim * dimTimeOut = NULL;
	if (dimTime != NULL) {
		CopyNcVar(*(vecFiles[0]), ncOutput, "time", true);
		dimTimeOut = ncOutput.get_dim("time");
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

	} else if (param.nAddTimeDim != -1) {
		dimTimeOut = ncOutput.add_dim("time", 0);
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error creating dimension \"time\" in output file");
		}

		NcVar * varTimeOut = ncOutput.add_var("time", ncDouble, dimTimeOut);
		if (varTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

		double dTime = static_cast<double>(param.nAddTimeDim);
		varTimeOut->set_cur((long)0);
		varTimeOut->put(&dTime, 1);

		if (param.strAddTimeDimUnits != "") {
			varTimeOut->add_att("units", param.strAddTimeDimUnits.c_str());
		}
		varTimeOut->add_att("long_name", "time");
		varTimeOut->add_att("calendar", "standard");
		varTimeOut->add_att("standard_name", "time");
	}

	CopyNcVar(*(vecFiles[0]), ncOutput, param.strLatitudeName.c_str(), true);
	CopyNcVar(*(vecFiles[0]), ncOutput, param.strLongitudeName.c_str(), true);

	NcDim * dimLonOut = ncOutput.get_dim(param.strLongitudeName.c_str());
	if (dimLonOut == NULL) {
		_EXCEPTION1("Error copying variable \"%s\" to output file",
			param.strLongitudeName.c_str());
	}
	NcDim * dimLatOut = ncOutput.get_dim(param.strLatitudeName.c_str());
	if (dimLatOut == NULL) {
		_EXCEPTION1("Error copying variable \"%s\" to output file",
			param.strLatitudeName.c_str());
	}

	// Tagged cell array
	DataArray2D<int> dVarDatatag(dimLat->size(), dimLon->size());

	NcVar * varIWVtag = NULL;
	if (dimTimeOut != NULL) {
		varIWVtag = ncOutput.add_var(
			"ar_binary_tag",
			ncByte,
			dimTimeOut,
			dimLatOut,
			dimLonOut);

	} else {
		varIWVtag = ncOutput.add_var(
			"ar_binary_tag",
			ncByte,
			dimLatOut,
			dimLonOut);
	}

	// Laplacian
	DataArray2D<double> dLaplacian(dimLat->size(), dimLon->size());

	NcVar * varLaplacian = NULL;
	if (param.fOutputLaplacian) {
		if (dimTimeOut != NULL) {
			varLaplacian = ncOutput.add_var(
				"ar_dx2",
				ncDouble,
				dimTimeOut,
				dimLatOut,
				dimLonOut);

		} else {
			varLaplacian = ncOutput.add_var(
				"ar_dx2",
				ncDouble,
				dimLatOut,
				dimLonOut);
		}
	}

	AnnounceEndBlock("Done");

	// Delta longitude
	double dDeltaLon = (dLonDeg[1] - dLonDeg[0]) / 180.0 * M_PI;
	double dDeltaLat = (dLatDeg[1] - dLatDeg[0]) / 180.0 * M_PI;

	double dX = dDeltaLon * static_cast<double>(param.iLaplacianSize);
	double dY = dDeltaLat * static_cast<double>(param.iLaplacianSize);

	double dX2 = dX * dX;
	double dY2 = dY * dY;

	// Loop through all times
	int nTimes = 1;
	if (dimTime != NULL) {
		nTimes = dimTime->size();
	}

	for (int t = 0; t < nTimes; t++) {

		char szBuffer[20];
		sprintf(szBuffer, "Time %i", t);
		AnnounceStartBlock(szBuffer);

		// Get the search-by variable array
		Variable & varSearchBy = varreg.Get(param.ixSearchByVar);
		varSearchBy.LoadGridData(varreg, vecFiles, grid, t);
		DataArray1D<float> & dataSearchBy = varSearchBy.GetData();

		DataArray2D<float> dVarData(dimLat->size(), dimLon->size(), false);
		dVarData.AttachToData(&(dataSearchBy[0]));

		dVarDatatag.Zero();

		// Calculate Laplacian of field at each point
		AnnounceStartBlock("Compute Laplacian");

		double dA = 1.0 / 12.0 * (1.0/dX2 + 1.0/dY2);
		double dB = 5.0 / (6.0 * dX2) - 1.0 / (6.0 * dY2);
		double dC = -1.0 / (6.0 * dX2) + 5.0 / (6.0 * dY2);
		double dD = -5.0 / 3.0 * (1.0/dX2 + 1.0/dY2);

		int j_begin = param.iLaplacianSize;
		int j_end = dimLat->size() - param.iLaplacianSize;

		int i_begin = 0;
		int i_end = dimLon->size();

		if (param.fRegional) {
			i_begin = param.iLaplacianSize;
			i_end = dimLon->size() - param.iLaplacianSize;
		}

		for (int j = j_begin; j < j_end; j++) {
		for (int i = i_begin; i < i_end; i++) {
			int i0 = (i + dimLon->size() - param.iLaplacianSize) % (dimLon->size());
			int i2 = (i + param.iLaplacianSize) % (dimLon->size());

			int j0 = j - param.iLaplacianSize;
			int j2 = j + param.iLaplacianSize;

			dLaplacian[j][i] =
				  dA * dVarData[j0][i0]
				+ dB * dVarData[j ][i0]
				+ dA * dVarData[j2][i0]
				+ dC * dVarData[j0][i ]
				+ dD * dVarData[j ][i ]
				+ dC * dVarData[j2][i ]
				+ dA * dVarData[j0][i2]
				+ dB * dVarData[j ][i2]
				+ dA * dVarData[j2][i2];
		}
		}

		AnnounceEndBlock("Done");

		AnnounceStartBlock("Build tagged cell array");

		// Build tagged cell array
		for (int j = 0; j < dimLat->size(); j++) {
		for (int i = 0; i < dimLon->size(); i++) {
			if (fabs(dLatDeg[j]) < param.dMinAbsLat) {
				continue;
			}
			if (dVarData[j][i] < param.dMinIWV) {
				continue;
			}
			if (dLaplacian[j][i] > -param.dMinLaplacian) {
				continue;
			}

			dVarDatatag[j][i] = 1;
		}
		}

		// Has zonal weights
		if (fHasZonalWeight) {

			// Compute zonal threshold
			DataArray1D<float> dZonalThreshold(dimLat->size());
			for (int j = 0; j < dimLat->size(); j++) {
				float dMaxZonalIWV = dVarData[j][0];
				for (int i = 0; i < dimLon->size(); i++) {
					dZonalThreshold[j] += dVarData[j][i];
					if (dVarData[j][i] > dMaxZonalIWV) {
						dMaxZonalIWV = dVarData[j][i];
					}
				}
				dZonalThreshold[j] /= static_cast<float>(dimLon->size());

				dZonalThreshold[j] =
					param.dZonalMeanWeight * dZonalThreshold[j]
					+ param.dZonalMaxWeight * dMaxZonalIWV;
			}

			// Adjust the tag
			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dVarData[j][i] < dZonalThreshold[j]) {
					dVarDatatag[j][i] = 0;
				}
			}
			}
		}

		// Has meridional weights
		if (fHasMeridWeight) {

			// Compute meridional threshold
			DataArray1D<float> dMeridThreshold(dimLon->size());
			for (int i = 0; i < dimLon->size(); i++) {
				float dMaxMeridVarData = dVarData[0][i];
				for (int j = 0; j < dimLat->size(); j++) {
					dMeridThreshold[i] += dVarData[j][i];
					if (dVarData[j][i] > dMaxMeridVarData) {
						dMaxMeridVarData = dVarData[j][i];
					}
				}
				dMeridThreshold[i] /= static_cast<float>(dimLon->size());

				dMeridThreshold[i] =
					param.dMeridMeanWeight * dMeridThreshold[i]
					+ param.dMeridMaxWeight * dMaxMeridVarData;
			}

			// Adjust the tag
			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dVarData[j][i] < dMeridThreshold[i]) {
					dVarDatatag[j][i] = 0;
				}
			}
			}
		}

		// Only retain blobs with minimum area
		if (param.nMinArea > 0) {
			std::set<Point> setBlobs;

			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dVarDatatag[j][i] != 0) {
					setBlobs.insert(std::pair<int,int>(j,i));
				}
			}
			}

			for (;;) {
				if (setBlobs.size() == 0) {
					break;
				}

				std::set<Point> setThisBlob;
				std::set<Point> setPointsToExplore;

				Point pt = *(setBlobs.begin());
				setBlobs.erase(setBlobs.begin());
				setPointsToExplore.insert(pt);

				for (;;) {
					if (setPointsToExplore.size() == 0) {
						break;
					}
					pt = *(setPointsToExplore.begin());
					setThisBlob.insert(pt);
					setPointsToExplore.erase(setPointsToExplore.begin());

					for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						Point pt2(pt.first + j, pt.second + i);
						if (pt2.first < 0) {
							continue;
						}
						if (pt2.first >= dimLat->size()) {
							continue;
						}
						if (pt2.second < 0) {
							pt2.second += dimLon->size();
						}
						if (pt2.second >= dimLon->size()) {
							pt2.second -= dimLon->size();
						}

						std::set<Point>::iterator iter = setBlobs.find(pt2);
						if (iter != setBlobs.end()) {
							setPointsToExplore.insert(*iter);
							setBlobs.erase(iter);
						}
					}
					}
				}

				if (setThisBlob.size() < param.nMinArea) {
					std::set<Point>::iterator iter = setThisBlob.begin();
					for (; iter != setThisBlob.end(); iter++) {
						dVarDatatag[iter->first][iter->second] = 0;
					}
				}
			}
		}
		AnnounceEndBlock("Done");

		AnnounceStartBlock("Writing results");

		// Output tagged cell array
		if (dimTimeOut != NULL) {
			if (varLaplacian != NULL) {
				varLaplacian->set_cur(t, 0, 0);
				varLaplacian->put(&(dLaplacian[0][0]), 1, dimLatOut->size(), dimLonOut->size());
			}

			varIWVtag->set_cur(t, 0, 0);
			varIWVtag->put(&(dVarDatatag[0][0]), 1, dimLatOut->size(), dimLonOut->size());

		} else {
			if (varLaplacian != NULL) {
				varLaplacian->set_cur(0, 0);
				varLaplacian->put(&(dLaplacian[0][0]), dimLatOut->size(), dimLonOut->size());
			}

			varIWVtag->set_cur(0, 0);
			varIWVtag->put(&(dVarDatatag[0][0]), dimLatOut->size(), dimLonOut->size());
		}

		AnnounceEndBlock("Done");

		AnnounceEndBlock(NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {
	// Input file
	std::string strInputFile;

	// Input file list
	std::string strInputFileList;

	// Output file
	std::string strOutputFile;

	// Output file list
	std::string strOutputFileList;

	// Input integrated water vapor variable name
	std::string strSearchByVariable;

	// SpineARs parameter
	SpineARsParam arparam;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strInputFileList, "in_list", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFileList, "out_list", "");
		CommandLineString(strSearchByVariable, "var", "");
		CommandLineString(arparam.strOutputVariable, "outvar", "");
		CommandLineInt(arparam.iLaplacianSize, "laplaciansize", 5);
		CommandLineDouble(arparam.dMinLaplacian, "minlaplacian", 0.5e4);
		CommandLineDouble(arparam.dMinAbsLat, "minabslat", 15.0);
		CommandLineDouble(arparam.dEqBandMaxLat, "eqbandmaxlat", 15.0);
		CommandLineDouble(arparam.dMinIWV, "minval", 20.0);
		CommandLineDoubleD(arparam.dZonalMeanWeight, "zonalmeanwt", 0.0, "(suggested 0.7)");
		CommandLineDoubleD(arparam.dZonalMaxWeight, "zonalmaxwt", 0.0, "(suggested 0.3)");
		CommandLineDoubleD(arparam.dMeridMeanWeight, "meridmeanwt", 0.0, "(suggested 0.9)");
		CommandLineDoubleD(arparam.dMeridMaxWeight, "meridmaxwt", 0.0, "(suggested 0.1)");
		CommandLineInt(arparam.nMinArea, "minarea", 0);
		CommandLineInt(arparam.nAddTimeDim, "addtimedim", -1);
		CommandLineString(arparam.strAddTimeDimUnits, "addtimedimunits", "");
		CommandLineBool(arparam.fOutputLaplacian, "laplacianout");
		CommandLineBool(arparam.fRegional, "regional");
		CommandLineString(arparam.strLongitudeName, "lonname", "lon");
		CommandLineString(arparam.strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	AnnounceStartBlock("Loading data");

	// Check input
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in) or (--in_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in) or (--in_list)"
			" may be specified");
	}

	// Check output
	if ((strOutputFile.length() == 0) && (strOutputFileList.length() == 0)) {
		_EXCEPTIONT("No output data file (--out) or (--out_list)"
			" specified");
	}
	if ((strOutputFile.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_list)"
			" may be specified");
	}

	// Check input/output
	if ((strInputFileList.length() != 0) && (strOutputFileList.length() == 0)) {
		_EXCEPTIONT("Arguments (--in_list) and (--out_list) must be specified together");
	}
	if ((strInputFile.length() != 0) && (strOutputFile.length() == 0)) {
		_EXCEPTIONT("Arguments (--in) and (--out) must be specified together");
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

	// Load output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFile.length() != 0) {
		vecOutputFiles.push_back(strOutputFile);

	} else {
		std::ifstream ifOutputFileList(strOutputFileList.c_str());
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputFiles.push_back(strFileLine);
		}
	}

	// Check length
	if (vecInputFiles.size() != vecOutputFiles.size()) {
		_EXCEPTIONT("Input and output file list length mismatch");
	}

	// Check variable
	if (strSearchByVariable == "") {
		_EXCEPTIONT("No search-by variable name (--var) specified");
	}

	// Check output variable
	if (arparam.strOutputVariable.length() == 0) {
		arparam.strOutputVariable = strSearchByVariable + "tag";
	}

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
#endif

	// Create Variable registry
	VariableRegistry varreg;

	// Get the integrated water vapor variable
	Variable varSearchByArg;
	varSearchByArg.ParseFromString(varreg, strSearchByVariable);
	arparam.ixSearchByVar = varreg.FindOrRegister(varSearchByArg);

	// Loop over all files to be processed
	for (int f = 0; f < vecInputFiles.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
#endif

		// Detect atmospheric rivers
		SpineARs(
			f,
			vecInputFiles[f],
			vecOutputFiles[f],
			varreg,
			arparam);
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

