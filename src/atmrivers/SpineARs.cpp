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

#include <algorithm>
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

float SampleCurve(
	double dLon,
	double dLat,
	double dLon0,
	double dLat0,
	const DataArray1D<float> & dParam
) {
	double dYp = (dLat - dLat0);
	//double dYp = (dLat - dLat0) * cos(dParam[2]) - (dLon - dLon0) * sin(dParam[2]);
	//return static_cast<float>(dParam[0] + dParam[1] * exp(-dYp * dYp / (dParam[3] * dParam[3])));
	return static_cast<float>(600.0 + 100.0 * exp(-dYp * dYp / (dParam[3] * dParam[3])));
}

///////////////////////////////////////////////////////////////////////////////

float SampleCurveRMSE(
	const DataArray2D<float> & dVar,
	const DataArray1D<double> & dLon,
	const DataArray1D<double> & dLat,
	int iLon,
	int iLat,
	int iFitRadius,
	double dFalloff,
	const DataArray1D<float> & dParam
) {
	if (iLat < iFitRadius) {
		_EXCEPTION();
	}
	if (iLat >= dLat.GetRows() - iFitRadius) {
		_EXCEPTION();
	}

	const double dLon0 = dLon[iLon];
	const double dLat0 = dLat[iLat];

	float dRMSE = 0.0;
	for (int jp = -iFitRadius; jp <= iFitRadius; jp++) {
		double dLatX = dLat[iLat + jp];
		for (int ip = -iFitRadius; ip <= iFitRadius; ip++) {
			int ix = iLon + ip;
			if (ix < 0) {
				ix += dLon.GetRows();
			}
			if (ix >= dLat.GetRows()) {
				ix -= dLon.GetRows();
			}
			double dLonX = dLon[ip];

			float dSample = SampleCurve(dLonX, dLatX, dLon0, dLat0, dParam) - dVar[iLat][ix];
/*
			float dDist = 2.0;
			if ((ip != 0) || (jp != 0)) {
				dDist = exp(-0.5 * dFalloff * log(
					static_cast<float>(jp) * static_cast<float>(jp)
					+ static_cast<float>(ip) * static_cast<float>(ip)));
			}

			dRMSE += dSample * dSample * dDist;
*/
			dRMSE += dSample * dSample;
		}
	}

	return dRMSE;
}

///////////////////////////////////////////////////////////////////////////////

void CurveFit(
	const DataArray2D<float> & dVar,
	const DataArray1D<double> & dLon,
	const DataArray1D<double> & dLat,
	int iLon,
	int iLat,
	int iFitRadius,
	double dFalloff,
	DataArray1D<float> & dParam
) {
	// Maximum number of iterations
	static const int MaxIterations = 100;

	// Allocate 4 entries
	dParam.Allocate(4);
	dParam[0] = 0.8 * dVar[iLat][iLon];
	dParam[1] = 0.2 * dVar[iLat][iLon];
	dParam[2] = 0.0;
	dParam[3] = 5.0 * M_PI / 180.0;

	// Temp param
	DataArray1D<float> dTempParam(4);

	// Approximate delta
	DataArray1D<float> dDelta(4);
	dDelta[0] = 1.0;
	dDelta[1] = 1.0;
	dDelta[2] = 0.0001;
	dDelta[3] = 0.0001;

	// Gamma in each direction
	DataArray1D<float> dGamma(4);
	dGamma[0] = 0.001;
	dGamma[1] = 0.001;
	dGamma[2] = 0.00000001;
	dGamma[3] = 0.0000001;

	// Gradient in each direction
	DataArray1D<float> dGrad(4);

	// Residual at each iteration
	for (int i = 0; i < MaxIterations; i++) {

		// Value of the function at this point
		float dVal =
			SampleCurveRMSE(
				dVar, dLon, dLat, iLon, iLat,
				iFitRadius,
				dFalloff,
				dParam);

		printf("%i %1.6e %i %i - %1.6e %1.6e %1.6e %1.6e\n", i, dVal, iLon, iLat, dParam[0], dParam[1], dParam[2], dParam[3]);

		// Calculate gradient in each direction
		for (int p = 0; p < dParam.GetRows(); p++) {
			dTempParam = dParam;
			dTempParam[p] += dDelta[p];

			dGrad[p] =
				SampleCurveRMSE(
					dVar, dLon, dLat, iLon, iLat,
					iFitRadius,
					dFalloff,
					dTempParam);

			printf("%1.6e %1.6e %1.6e %1.6e -- %1.6e\n", dTempParam[0], dTempParam[1], dTempParam[2], dTempParam[3], dGrad[p]);

			dGrad[p] = (dGrad[p] - dVal) / dDelta[p];
		}

		// Gradient descent
		for (int p = 0; p < dParam.GetRows(); p++) {
			printf("%1.6e %1.6e %1.6e %1.6e\n", dGrad[0], dGrad[1], dGrad[2], dGrad[3]);
			dParam[p] -= dGamma[p] * dGrad[p];
		}
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
		iSpineWidth(-1),
		nSpineOrientationDist(3),
		fOutputSpineOrientation(false),
		dMinLaplacian(0.0),
		fOutputLaplacian(false),
		dMinAbsGrad(0.0),
		fOutputAbsGrad(false),
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

	// Spine width
	int iSpineWidth;

	// Spine orientation distance
	int nSpineOrientationDist;

	// Output spine orientation
	bool fOutputSpineOrientation;

	// Minimum Laplacian
	double dMinLaplacian;

	// Output Laplacian
	bool fOutputLaplacian;

	// Minimum absolute gradient
	double dMinAbsGrad;

	// Output minimum absolute gradient
	bool fOutputAbsGrad;

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
	const SpineARsParam & param,
	bool fVerbose
) {
	if (fVerbose) AnnounceStartBlock("Processing files \"%s\"", strInputFiles.c_str());

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
	if (fVerbose) AnnounceStartBlock("Generating RLL grid data");
	grid.GenerateLatitudeLongitude(
		vecFiles[0],
		param.fRegional,
		param.strLatitudeName,
		param.strLongitudeName);
	if (fVerbose) AnnounceEndBlock("Done");

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

	DataArray1D<double> dLonRad(dimLon->size());
	for (int i = 0; i < dLonDeg.GetRows(); i++) {
		dLonRad[i] = dLonDeg[i] * M_PI / 180.0;
	}

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

	DataArray1D<double> dLatRad(dimLat->size());
	for (int i = 0; i < dLatDeg.GetRows(); i++) {
		dLatRad[i] = dLatDeg[i] * M_PI / 180.0;
	}

	// Check iSpineWidth
	if (param.iSpineWidth != (-1)) {
		if (param.iSpineWidth >= dimLon->size()) {
			_EXCEPTIONT("spinewidth must be < longitude dimension size");
		}
		if (param.iSpineWidth >= dimLat->size()) {
			_EXCEPTIONT("spinewidth must be < latitude dimension size");
		}
	}

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

	// Absolute gradient
	DataArray2D<double> dAbsGrad;
	if (param.dMinAbsGrad != 0.0) {
		dAbsGrad.Allocate(dimLat->size(), dimLon->size());
	}

	NcVar * varAbsGrad = NULL;
	if (param.fOutputAbsGrad) {
		if (dimTimeOut != NULL) {
			varAbsGrad = ncOutput.add_var(
				"ar_absgrad",
				ncDouble,
				dimTimeOut,
				dimLatOut,
				dimLonOut);

		} else {
			varAbsGrad = ncOutput.add_var(
				"ar_absgrad",
				ncDouble,
				dimLatOut,
				dimLonOut);
		}
	}

	// Orientation variable
	DataArray2D<double> dOrientation;

	NcVar * varOrientation = NULL;
	if (param.fOutputSpineOrientation) {
		dOrientation.Allocate(dimLat->size(), dimLon->size());

		if (dimTimeOut != NULL) {
			varOrientation =
				ncOutput.add_var(
					"ar_orientation",
					ncDouble,
					dimTimeOut,
					dimLatOut,
					dimLonOut);

		} else {
			varOrientation =
				ncOutput.add_var(
					"ar_orientation",
					ncDouble,
					dimLatOut,
					dimLonOut);
		}
		varOrientation->add_att("_FillValue", -1.0);
	}

	// Delta longitude
	double dDeltaLon = fabs(dLonDeg[1] - dLonDeg[0]) / 180.0 * M_PI;
	double dDeltaLat = fabs(dLatDeg[1] - dLatDeg[0]) / 180.0 * M_PI;

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

		if (fVerbose) AnnounceStartBlock("Time %i", t);

		// Get the search-by variable array
		Variable & varSearchBy = varreg.Get(param.ixSearchByVar);
		varSearchBy.LoadGridData(varreg, vecFiles, grid, t);
		DataArray1D<float> & dataSearchBy = varSearchBy.GetData();

		DataArray2D<float> dVarData(dimLat->size(), dimLon->size(), false);
		dVarData.AttachToData(&(dataSearchBy[0]));

		dVarDatatag.Zero();

		// Calculate Laplacian of field at each point
		if (fVerbose) AnnounceStartBlock("Compute Laplacian");

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

		if (fVerbose) AnnounceEndBlock("Done");

		if (param.dMinAbsGrad != 0.0) {
			if (fVerbose) AnnounceStartBlock("Compute AbsGrad");

			for (int j = j_begin; j < j_end; j++) {
			for (int i = i_begin; i < i_end; i++) {
				int i0 = (i + dimLon->size() - param.iLaplacianSize) % (dimLon->size());
				int i2 = (i + param.iLaplacianSize) % (dimLon->size());

				int j0 = j - param.iLaplacianSize;
				int j2 = j + param.iLaplacianSize;

				dAbsGrad[j][i] =
				    fabs(dVarData[j2][i2] - dVarData[j0][i0]) / 2.0 / sqrt(dX2 + dY2)
				  + fabs(dVarData[j2][i0] - dVarData[j0][i2]) / 2.0 / sqrt(dX2 + dY2)
				  + fabs(dVarData[j ][i2] - dVarData[j ][i0]) / 2.0 / dX
				  + fabs(dVarData[j2][i ] - dVarData[j0][i ]) / 2.0 / dY;
			}
			}

			if (fVerbose) AnnounceEndBlock("Done");
		}

		if (fVerbose) AnnounceStartBlock("Build tagged cell array");

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

			// Tag spines
			if (param.iSpineWidth != (-1)) {

				// Check for local maximum in each of the four directions
				bool fLocalMax[4] = {true, true, true, true};

				for (int k = -param.iSpineWidth; k <= param.iSpineWidth; k++) {

					if (k == 0) {
						continue;
					}

					int ix = i + k;
					if (ix < 0) {
						ix += dimLon->size();
					}
					if (ix >= dimLon->size()) {
						ix -= dimLon->size();
					}

					int jx1 = j + k;
					if (jx1 < 0) {
						jx1 = 0;
					}
					if (jx1 >= dimLat->size()) {
						jx1 = dimLat->size()-1;
					}

					int jx2 = j - k;
					if (jx2 < 0) {
						jx2 = 0;
					}
					if (jx2 >= dimLat->size()) {
						jx2 = dimLat->size()-1;
					}
/*
					if ((i == 175) && (j == 238)) {
						printf("\n");
						printf("(%1.10e)\n", dVarData[j][i]);
						printf("%i %i (%1.10e)\n", jx1, i, dVarData[jx1][i]);
						printf("%i %i (%1.10e)\n", j, ix, dVarData[j][ix]);
						printf("%i %i (%1.10e)\n", jx1, ix, dVarData[jx1][ix]);
						printf("%i %i (%1.10e)\n", jx2, ix, dVarData[jx2][ix]);
					}
*/
					if (dVarData[jx1][i] >= dVarData[j][i]) {
						fLocalMax[0] = false;
					}
					if (dVarData[j][ix] >= dVarData[j][i]) {
						fLocalMax[1] = false;
					}
					//if (dVarData[jx1][ix] >= dVarData[j][i]) {
					//	fLocalMax[2] = false;
					//}
					//if (dVarData[jx2][ix] >= dVarData[j][i]) {
					//	fLocalMax[3] = false;
					//}
				}

				if (fLocalMax[0] || fLocalMax[1]) {
					//if ((i == 175) && (j == 238)) {
					//	std::cout << "TEST " << j << std::endl;
					//}
					dVarDatatag[j][i] = 2;
/*
					if ((i == 225) && (j == 257)) {
						DataArray1D<float> dParam(4);
						CurveFit(
							dVarData,
							dLonRad,
							dLatRad,
							i,
							j,
							4,
							1.0,
							dParam);
					}
*/
					//_EXCEPTION();

				}

				//if (fLocalMax[0] || fLocalMax[1] || fLocalMax[2] || fLocalMax[3]) {
				//	dVarDatatag[j][i] = 2;
				//}
			}
		}
		}

		// Check absolute gradient criteira
		if (param.dMinAbsGrad != 0.0) {
			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dAbsGrad[j][i] < param.dMinAbsGrad) {
					dVarDatatag[j][i] = 0;
				}
			}
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
		if (fVerbose) AnnounceEndBlock("Done");

		// Find the orientation at each spine node
		if (param.fOutputSpineOrientation) {

			if (fVerbose) AnnounceStartBlock("Calculating spine orientation");

			// Create a new orientation array
			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				dOrientation[j][i] = -1.0;
			}
			}

			// Find all spine points
			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dVarDatatag[j][i] == 2) {

					// Build the set of contiguous spine points
					std::set< std::pair<int,int> > setSpinePointsVisited;
					std::vector< std::pair<int,int> > vecSpinePointsToVisit;
					std::vector< std::pair<int,int> > vecSpinePointsToVisitNext;
					vecSpinePointsToVisit.push_back(std::pair<int,int>(j,i));

					for (int p = 0; p < param.nSpineOrientationDist; p++) {
						for (int q = 0; q < vecSpinePointsToVisit.size(); q++) {
							setSpinePointsVisited.insert(vecSpinePointsToVisit[q]);
							for (int m = -1; m <= 1; m++) {
							for (int n = -1; n <= 1; n++) {
								int jx = vecSpinePointsToVisit[q].first + m;
								if ((jx < 0) || (jx >= dimLat->size())) {
									continue;
								}
								int ix = vecSpinePointsToVisit[q].second + n;
								if (ix < 0) {
									ix += dimLon->size();
								}
								if (ix >= dimLon->size()) {
									ix -= dimLon->size();
								}
								if (dVarDatatag[jx][ix] == 2) {
									std::pair<int,int> px(jx,ix);
									if (setSpinePointsVisited.find(px) ==
										setSpinePointsVisited.end()
									) {
										vecSpinePointsToVisitNext.push_back(px);
									}
								}
							}
							}
						}
						vecSpinePointsToVisit = vecSpinePointsToVisitNext;
						vecSpinePointsToVisitNext.clear();
					}

					// Remove isolated spine points
					if (setSpinePointsVisited.size() == 1) {
						dVarDatatag[j][i] = 1;
					}

					// Calculate the orientation of these points
					double dInvN = 1.0 / static_cast<double>(setSpinePointsVisited.size());
					double dSumX = 0.0;
					double dSumY = 0.0;
					double dSumX2 = 0.0;
					double dSumXY = 0.0;

					for (
						auto iter = setSpinePointsVisited.begin();
						iter != setSpinePointsVisited.end(); iter++
					) {
						double dY = dLatRad[iter->first];
						double dX = dLonRad[iter->second];

						dSumX += dX;
						dSumY += dY;
						dSumX2 += dX * dX;
						dSumXY += dX * dY;
					}

					double dCovXY = dSumXY - dInvN * dSumX * dSumY;
					double dVarX = dSumX2 - dInvN * dSumX * dSumX;

					if (dVarX < 1.0e-12) {
						dOrientation[j][i] = 90.0;
					} else {
						dOrientation[j][i] = atan2(dCovXY, dVarX) * 180.0 / M_PI;
					}
				}
			}
			}

			if (fVerbose) AnnounceEndBlock("Done");
		}

		if (fVerbose) AnnounceStartBlock("Writing results");

		// Output tagged cell array
		if (dimTimeOut != NULL) {
			if (varLaplacian != NULL) {
				varLaplacian->set_cur(t, 0, 0);
				varLaplacian->put(&(dLaplacian[0][0]), 1, dimLatOut->size(), dimLonOut->size());
			}
			if (varAbsGrad != NULL) {
				varAbsGrad->set_cur(t, 0, 0);
				varAbsGrad->put(&(dAbsGrad[0][0]), 1, dimLatOut->size(), dimLonOut->size());
			}
			if (varOrientation != NULL) {
				varOrientation->set_cur(t, 0, 0);
				varOrientation->put(&(dOrientation[0][0]), 1, dimLatOut->size(), dimLonOut->size());
			}

			varIWVtag->set_cur(t, 0, 0);
			varIWVtag->put(&(dVarDatatag[0][0]), 1, dimLatOut->size(), dimLonOut->size());

		} else {
			if (varLaplacian != NULL) {
				varLaplacian->set_cur(0, 0);
				varLaplacian->put(&(dLaplacian[0][0]), dimLatOut->size(), dimLonOut->size());
			}
			if (varAbsGrad != NULL) {
				varAbsGrad->set_cur(0, 0);
				varAbsGrad->put(&(dAbsGrad[0][0]), dimLatOut->size(), dimLonOut->size());
			}
			if (varOrientation != NULL) {
				varOrientation->set_cur(t, 0, 0);
				varOrientation->put(&(dOrientation[0][0]), dimLatOut->size(), dimLonOut->size());
			}

			varIWVtag->set_cur(0, 0);
			varIWVtag->put(&(dVarDatatag[0][0]), dimLatOut->size(), dimLonOut->size());
		}

		if (fVerbose) AnnounceEndBlock("Done");

		if (fVerbose) AnnounceEndBlock(NULL);
	}

	if (fVerbose) AnnounceEndBlock("Done");
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

	// Write logs
	std::string strLogDir;

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
		CommandLineInt(arparam.iSpineWidth, "spinewidth", -1);
		CommandLineDouble(arparam.dMinLaplacian, "minlaplacian", 0.5e4);
		CommandLineBool(arparam.fOutputLaplacian, "outputlaplacian");
		CommandLineInt(arparam.nSpineOrientationDist, "spineorientationdist", 3);
		CommandLineBool(arparam.fOutputSpineOrientation, "outputspineorientation");
		CommandLineDouble(arparam.dMinAbsGrad, "minabsgrad", 0.0);
		CommandLineBool(arparam.fOutputAbsGrad, "outputabsgrad");
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
		CommandLineBool(arparam.fRegional, "regional");
		CommandLineString(arparam.strLongitudeName, "lonname", "lon");
		CommandLineString(arparam.strLatitudeName, "latname", "lat");
		CommandLineString(strLogDir, "logdir", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

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

	AnnounceStartBlock("Initializing detector");

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

	Announce("Executing detection with %i threads over %i files",
		nMPISize, vecInputFiles.size());
#endif

	// Create Variable registry
	VariableRegistry varreg;

	// Get the integrated water vapor variable
	Variable varSearchByArg;
	varSearchByArg.ParseFromString(varreg, strSearchByVariable);
	arparam.ixSearchByVar = varreg.FindOrRegister(varSearchByArg);

	AnnounceEndBlock("Done");

	// Loop over all files to be processed
	AnnounceStartBlock("Performing detection");

#if defined(TEMPEST_MPIOMP)
	if (strLogDir == "") {
		Announce("Reporting disabled (if reporting desired, use --logdir)");
	} else {
		Announce("Logs will be written to \"%s\"", strLogDir.c_str());
	}
#endif

#if defined(TEMPEST_MPIOMP)
	// Open log file
	if (strLogDir != "") {
		std::string strLogFile = strLogDir;
		if (strLogFile[strLogFile.length()-1] != '/') {
			strLogFile += "/";
		}
		char szTemp[20];
		sprintf(szTemp, "log%06i.txt", nMPIRank);
		strLogFile += szTemp;

		arparam.fpLog = fopen(strLogFile.c_str(), "w");
	}
#else
	if (strLogDir != "") {
		std::string strLogFile = strLogDir;
		if (strLogFile[strLogFile.length()-1] != '/') {
			strLogFile += "/log000000.txt";
		}

		arparam.fpLog = fopen(strLogFile.c_str(), "w");
	}
#endif

	for (int f = 0; f < vecInputFiles.size(); f++) {
		bool fVerbose = true;

#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
		if (nMPISize > 1) {
			fVerbose = false;
		}
#endif

		// Detect atmospheric rivers
		SpineARs(
			f,
			vecInputFiles[f],
			vecOutputFiles[f],
			varreg,
			arparam,
			fVerbose);
	}

	if (arparam.fpLog != NULL) {
		fclose(arparam.fpLog);
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

