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

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse the list of input files.
///	</summary>
void ParseInputFiles(
	const std::string & strInputFiles,
	std::vector<NcFile *> & vecFiles
) {
	int iLast = 0;
	for (int i = 0; i <= strInputFiles.length(); i++) {
		if ((i == strInputFiles.length()) ||
		    (strInputFiles[i] == ';')
		) {
			std::string strFile =
				strInputFiles.substr(iLast, i - iLast);

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
			strInputFiles.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input file
	std::string strInputFiles;

	// Output file
	std::string strOutputFile;

	// Input integrated water vapor variable name
	std::string strSearchByVariable;

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

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFiles, "in", "");
		//CommandLineString(strInputFilesList, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strSearchByVariable, "var", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineInt(iLaplacianSize, "laplaciansize", 5);
		CommandLineDouble(dMinLaplacian, "minlaplacian", 0.5e4);
		CommandLineDouble(dMinAbsLat, "minabslat", 15.0);
		CommandLineDouble(dEqBandMaxLat, "eqbandmaxlat", 15.0);
		CommandLineDouble(dMinIWV, "minval", 20.0);
		CommandLineDouble(dZonalMeanWeight, "zonalmeanwt", 0.7);
		CommandLineDouble(dZonalMaxWeight, "zonalmaxwt", 0.3);
		CommandLineDouble(dMeridMeanWeight, "meridmeanwt", 0.9);
		CommandLineDouble(dMeridMaxWeight, "meridmaxwt", 0.1);
		CommandLineInt(nMinArea, "minarea", 0);
		CommandLineInt(nAddTimeDim, "addtimedim", -1);
		CommandLineString(strAddTimeDimUnits, "addtimedimunits", "");
		CommandLineBool(fOutputLaplacian, "laplacianout");
		CommandLineBool(fRegional, "regional");
		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	AnnounceStartBlock("Loading data");

	// Create Variable registry
	VariableRegistry varreg;

	// Check input
	if (strInputFiles == "") {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check variable
	if (strSearchByVariable == "") {
		_EXCEPTIONT("No IWV variable name (--iwvvar) specified");
	}

	// Check output variable
	if (strOutputVariable.length() == 0) {
		strOutputVariable = strSearchByVariable + "tag";
	}

	// Load in the benchmark file
	NcFileVector vecFiles;

	ParseInputFiles(strInputFiles, vecFiles);

	// Define the SimpleGrid
	SimpleGrid grid;
	AnnounceStartBlock("Generating RLL grid data");
	grid.GenerateLatitudeLongitude(
		vecFiles[0],
		fRegional,
		strLatitudeName,
		strLongitudeName);
	AnnounceEndBlock("Done");

	// Get the time dimension
	NcDim * dimTime = vecFiles[0]->get_dim("time");
	//if (dimTime == NULL) {
	//	_EXCEPTIONT("Error accessing dimension \"time\"");
	//}

	// Get the longitude dimension
	NcDim * dimLon = vecFiles[0]->get_dim(strLongitudeName.c_str());
	if (dimLon == NULL) {
		_EXCEPTION1("Error accessing dimension \"%s\"", strLongitudeName.c_str());
	}

	// Get the longitude variable
	NcVar * varLon = vecFiles[0]->get_var(strLongitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"", strLongitudeName.c_str());
	}

	DataArray1D<double> dLonDeg(dimLon->size());
	varLon->get(&(dLonDeg[0]), dimLon->size());

	// Get the latitude dimension
	NcDim * dimLat = vecFiles[0]->get_dim(strLatitudeName.c_str());
	if (dimLat == NULL) {
		_EXCEPTION1("Error accessing dimension \"%s\"", strLatitudeName.c_str());
	}

	// Get the latitude variable
	NcVar * varLat = vecFiles[0]->get_var(strLatitudeName.c_str());
	if (varLat == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"", strLatitudeName.c_str());
	}

	DataArray1D<double> dLatDeg(dimLat->size());
	varLat->get(&(dLatDeg[0]), dimLat->size());

	// Get the integrated water vapor variable
	Variable varSearchByArg;
	varSearchByArg.ParseFromString(varreg, strSearchByVariable);
	int ixSearchByVar = varreg.FindOrRegister(varSearchByArg);

	DataArray2D<float> dVarData(dimLat->size(), dimLon->size());

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

	} else if (nAddTimeDim != -1) {
		dimTimeOut = ncOutput.add_dim("time", 0);
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error creating dimension \"time\" in output file");
		}

		NcVar * varTimeOut = ncOutput.add_var("time", ncDouble, dimTimeOut);
		if (varTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

		double dTime = static_cast<double>(nAddTimeDim);
		varTimeOut->set_cur((long)0);
		varTimeOut->put(&dTime, 1);

		if (strAddTimeDimUnits != "") {
			varTimeOut->add_att("units", strAddTimeDimUnits.c_str());
		}
		varTimeOut->add_att("long_name", "time");
		varTimeOut->add_att("calendar", "standard");
		varTimeOut->add_att("standard_name", "time");
	}

	CopyNcVar(*(vecFiles[0]), ncOutput, strLatitudeName.c_str(), true);
	CopyNcVar(*(vecFiles[0]), ncOutput, strLongitudeName.c_str(), true);

	NcDim * dimLonOut = ncOutput.get_dim(strLongitudeName.c_str());
	if (dimLonOut == NULL) {
		_EXCEPTION1("Error copying variable \"%s\" to output file", strLongitudeName.c_str());
	}
	NcDim * dimLatOut = ncOutput.get_dim(strLatitudeName.c_str());
	if (dimLatOut == NULL) {
		_EXCEPTION1("Error copying variable \"%s\" to output file", strLatitudeName.c_str());
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
	if (fOutputLaplacian) {
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

	double dX = dDeltaLon * static_cast<double>(iLaplacianSize);
	double dY = dDeltaLat * static_cast<double>(iLaplacianSize);

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
		Variable & varSearchBy = varreg.Get(ixSearchByVar);
		varSearchBy.LoadGridData(varreg, vecFiles, grid, t);
		DataArray1D<float> & dataSearchBy = varSearchBy.GetData();

		DataArray2D<float> dataSearchBy2D(dimLat->size(), dimLon->size(), false);
		dataSearchBy2D.AttachToData(&(dataSearchBy[0]));

		dVarDatatag.Zero();
/*
		// Identify ridges
		AnnounceStartBlock("Identify ridges");
		for (int j = iLaplacianSize; j < dimLat->size() - iLaplacianSize; j++) {
		for (int i = 0; i < dimLon->size(); i++) {
			const int nPoints = 32;
			for (int n = 0; n < nPoints; n++) {
				double dLonPt = 
			}
		}
		}
		AnnounceEndBlock("Done");
*/

		AnnounceStartBlock("Compute Laplacian");

		// Compute Laplacian
		double dA = 1.0 / 12.0 * (1.0/dX2 + 1.0/dY2);
		double dB = 5.0 / (6.0 * dX2) - 1.0 / (6.0 * dY2);
		double dC = -1.0 / (6.0 * dX2) + 5.0 / (6.0 * dY2);
		double dD = -5.0 / 3.0 * (1.0/dX2 + 1.0/dY2);

		int j_begin = iLaplacianSize;
		int j_end = dimLat->size() - iLaplacianSize;

		int i_begin = 0;
		int i_end = dimLon->size();

		if (fRegional) {
			i_begin = iLaplacianSize;
			i_end = dimLon->size() - iLaplacianSize;
		}

		for (int j = j_begin; j < j_end; j++) {
		for (int i = i_begin; i < i_end; i++) {
			int i0 = (i + dimLon->size() - iLaplacianSize) % (dimLon->size());
			int i2 = (i + iLaplacianSize) % (dimLon->size());

			int j0 = j - iLaplacianSize;
			int j2 = j + iLaplacianSize;

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
/*
			if (dLaplacian[j][i] > -dMinLaplacian) {
				dLaplacian[j][i] = 0.0;
			} else {
				dLaplacian[j][i] = 1.0;
			}
*/
		}
		}

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
				dZonalMeanWeight * dZonalThreshold[j]
				+ dZonalMaxWeight * dMaxZonalIWV;
		}

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
				dMeridMeanWeight * dMeridThreshold[i]
				+ dMeridMaxWeight * dMaxMeridVarData;
		}

		// Announce
		AnnounceEndBlock("Done");

		AnnounceStartBlock("Build tagged cell array");

		// Build tagged cell array
		for (int j = 0; j < dimLat->size(); j++) {
		for (int i = 0; i < dimLon->size(); i++) {
			if (fabs(dLatDeg[j]) < dMinAbsLat) {
				continue;
			}
			if (dVarData[j][i] < dMinIWV) {
				continue;
			}
			if (dVarData[j][i] < dZonalThreshold[j]) {
				continue;
			}
			if (dVarData[j][i] < dMeridThreshold[i]) {
				continue;
			}

			//dVarDatatag[j][i] = 1 + static_cast<int>(dLaplacian[j][i]);
			if (dLaplacian[j][i] > -dMinLaplacian) {
				dVarDatatag[j][i] = 0;
			} else {
				dVarDatatag[j][i] = 1;
			}
		}
		}

		// Only retain blobs with minimum area
		if (nMinArea > 0) {
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

				if (setThisBlob.size() < nMinArea) {
					std::set<Point>::iterator iter = setThisBlob.begin();
					for (; iter != setThisBlob.end(); iter++) {
						dVarDatatag[iter->first][iter->second] = 0;
					}
				}
			}
		}

/*
		// Remove points connected with equatorial moisture band
		for (int i = 0; i < dimLon->size(); i++) {
			bool fSouthDone = false;
			bool fNorthDone = false;

			for (int j = 0; j < dimLat->size(); j++) {
				if ((!fSouthDone) && (dLatDeg[j] > -dMinAbsLat)) {
					for (int k = j-1; k > 0; k--) {
						if (dVarDatatag[k][i] == 0) {
							break;
						}
						if (dLatDeg[k] < -dEqBandMaxLat) {
							break;
						}
						dVarDatatag[k][i] = 0;
					}
					fSouthDone = true;
				}
				if ((!fNorthDone) && (dLatDeg[j] > dMinAbsLat)) {
					for (int k = j-1; k < dimLat->size(); k++) {
						if (dVarDatatag[k][i] == 0) {
							break;
						}
						if (dLatDeg[k] > dEqBandMaxLat) {
							break;
						}
						dVarDatatag[k][i] = 0;
					}
					fNorthDone = true;
				}
			}
		}
*/
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

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

