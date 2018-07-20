///////////////////////////////////////////////////////////////////////////////
///
///	\file    SpineARs_v2.cpp
///	\author  Paul Ullrich
///	\version July 17th, 2018
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
#include "SimpleGrid.h"
#include "SparseMatrix.h"
#include "kdtree.h"

#include "DataVector.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <set>

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

void GenerateEqualDistanceSpherePoints(
	double dX0,
	double dY0,
	double dZ0,
	int nPoints,
	double dDist,
	std::vector<double> & dXout,
	std::vector<double> & dYout,
	std::vector<double> & dZout
) { 
	dXout.resize(nPoints);
	dYout.resize(nPoints);
	dZout.resize(nPoints);

	// Pick a quasi-arbitrary reference direction
	double dX1;
	double dY1;
	double dZ1;
	if ((fabs(dX0) >= fabs(dY0)) && (fabs(dX0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0 + 1.0;
		dZ1 = dZ0;
	} else if ((fabs(dY0) >= fabs(dX0)) && (fabs(dY0) >= fabs(dZ0))) {
		dX1 = dX0;
		dY1 = dY0;
		dZ1 = dZ0 + 1.0;
	} else {
		dX1 = dX0 + 1.0;
		dY1 = dY0;
		dZ1 = dZ0;
	}

	// Project perpendicular to detection location
	double dDot = dX1 * dX0 + dY1 * dY0 + dZ1 * dZ0;

	dX1 -= dDot * dX0;
	dY1 -= dDot * dY0;
	dZ1 -= dDot * dZ0;

	// Normalize
	double dMag1 = sqrt(dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1);

	if (dMag1 < 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	double dScale1 = tan(dDist * M_PI / 180.0) / dMag1;

	dX1 *= dScale1;
	dY1 *= dScale1;
	dZ1 *= dScale1;

	// Verify dot product is zero
	dDot = dX0 * dX1 + dY0 * dY1 + dZ0 * dZ1;
	if (fabs(dDot) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	// Cross product (magnitude automatically 
	double dCrossX = dY0 * dZ1 - dZ0 * dY1;
	double dCrossY = dZ0 * dX1 - dX0 * dZ1;
	double dCrossZ = dX0 * dY1 - dY0 * dX1;

	// Generate all points
	for (int j = 0; j < nPoints; j++) {

		// Angle of rotation
		double dAngle = 2.0 * M_PI
			* static_cast<double>(j)
			/ static_cast<double>(nPoints);

		// Calculate new rotated vector
		double dX2 = dX0 + dX1 * cos(dAngle) + dCrossX * sin(dAngle);
		double dY2 = dY0 + dY1 * cos(dAngle) + dCrossY * sin(dAngle);
		double dZ2 = dZ0 + dZ1 * cos(dAngle) + dCrossZ * sin(dAngle);

		double dMag2 = sqrt(dX2 * dX2 + dY2 * dY2 + dZ2 * dZ2);

		dX2 /= dMag2;
		dY2 /= dMag2;
		dZ2 /= dMag2;

		dXout[j] = dX2;
		dYout[j] = dY2;
		dZout[j] = dZ2;
	}
}

///////////////////////////////////////////////////////////////////////////////

void BuildLaplacianOperator(
	const SimpleGrid & grid,
	int nLaplacianPoints,
	double dLaplacianDist,
	double dLaplacianSmoothDistFactor,
	SparseMatrix<double> & opLaplacian
) {
	opLaplacian.Clear();

	int iRef = 0;
/*
	// Convert distance to radians
	dMinDist = dMinDist * M_PI / 180.0;
	dMaxDist = dMaxDist * M_PI / 180.0;

	// Convert distance to chord length
	dMinDist = 2.0 * sin(dMinDist / 2.0);
	dMaxDist = 2.0 * sin(dMaxDist / 2.0);

	// Convert distance to chord length squared
	double dMinDist2 = dMinDist * dMinDist;
	double dMaxDist2 = dMaxDist * dMaxDist;
*/
	// SPH smoothing distance
	double dH = dLaplacianDist * dLaplacianSmoothDistFactor;
	double dHScale = 3.0 / (2.0 * M_PI * dH * dH * dH);

	// Create a kdtree with all nodes in grid
	kdtree * kdGrid = kd_create(3);
	if (kdGrid == NULL) {
		_EXCEPTIONT("Error creating kdtree");
	}

	DataVector<double> dXi(grid.GetSize());
	DataVector<double> dYi(grid.GetSize());
	DataVector<double> dZi(grid.GetSize());
	for (int i = 0; i < grid.GetSize(); i++) {
		double dLat = grid.m_dLat[i] * M_PI / 180.0;
		double dLon = grid.m_dLon[i] * M_PI / 180.0;

		dXi[i] = cos(dLon) * cos(dLat);
		dYi[i] = sin(dLon) * cos(dLat);
		dZi[i] = sin(dLat);

		kd_insert3(kdGrid, dXi[i], dYi[i], dZi[i], (void*)((&iRef)+i));
	}

	// Construct the Laplacian operator using SPH
	for (int i = 0; i < grid.GetSize(); i++) {

		// Generate points for the Laplacian
		std::vector<double> dXout;
		std::vector<double> dYout;
		std::vector<double> dZout;

		GenerateEqualDistanceSpherePoints(
			dXi[i], dYi[i], dZi[i],
			nLaplacianPoints,
			dLaplacianDist,
			dXout, dYout, dZout);

		dXout.push_back(dXi[i]);
		dYout.push_back(dYi[i]);
		dZout.push_back(dZi[i]);
/*
		kdres * kdr = kd_nearest_range3(kdGrid, dXi[i], dYi[i], dZi[i], dMaxDist);
		if (kdr == NULL) {
			_EXCEPTIONT("Error in kd_nearest_range3");
		}
		int nNodes = kd_res_size(kdr);

		std::cout << nNodes << " found at distance " << dMaxDist << std::endl;
*/
		std::vector<int> vecCols;
		std::vector<double> vecWeight;
		double dVolume = 0.0;

		std::set<int> setPoints;

		for (int j = 0; j < nLaplacianPoints; j++) {
			// Find the nearest grid point to the output point
			kdres * kdr = kd_nearest3(kdGrid, dXout[j], dYout[j], dZout[j]);
			if (kdr == NULL) {
				_EXCEPTIONT("NULL return value in call to kd_nearest3");
			}

			void* pData = kd_res_item_data(kdr);
			if (pData == NULL) {
				_EXCEPTIONT("NULL data index");
			}
			int k = ((int*)(pData)) - (&iRef);

			kd_res_free(kdr);

			if (k > dXi.GetRows()) {
				_EXCEPTIONT("Invalid point index");
			}

			// Ensure points are not duplicated
			if (setPoints.find(k) != setPoints.end()) {
				continue;
			} else {
				setPoints.insert(k);
			}
	/*
			// ============== BEGIN DEBUGGING =============================
			double dLon0 = atan2(dYi[i], dXi[i]) * 180.0 / M_PI;
			double dLat0 = asin(dZi[i]) * 180.0 / M_PI;

			double dLon1 = atan2(dYout[j], dXout[j]) * 180.0 / M_PI;
			double dLat1 = asin(dZout[j]) * 180.0 / M_PI;

			double dDist0 =
				sqrt(
					(dXout[j] - dXi[i]) * (dXout[j] - dXi[i])
					+ (dYout[j] - dYi[i]) * (dYout[j] - dYi[i])
					+ (dZout[j] - dZi[i]) * (dZout[j] - dZi[i]));

			std::cout << 2.0 * sin(dDist0 / 2.0) * 180.0 / M_PI << std::endl;

			double dDist1 =
				sqrt(
					(dXi[k] - dXi[i]) * (dXi[k] - dXi[i])
					+ (dYi[k] - dYi[i]) * (dYi[k] - dYi[i])
					+ (dZi[k] - dZi[i]) * (dZi[k] - dZi[i]));

			std::cout << 2.0 * sin(dDist1 / 2.0) * 180.0 / M_PI << std::endl;

			double dLon2 = atan2(dYi[k], dXi[k]) * 180.0 / M_PI;
			double dLat2 = asin(dZi[k]) * 180.0 / M_PI;

			printf("XY: %1.3f %1.3f :: %1.3f %1.3f :: %1.3f %1.3f\n", dXi[i], dYi[i], dXout[j], dYout[j], dXi[k], dYi[k]);
			printf("LL: %1.2f %1.2f :: %1.2f %1.2f\n", dLon0, dLat0, dLon1, dLat1);
			// ============== END DEBUGGING =============================
*/

			// Kernel derivatives (cubic kernel)
			double dX1 = dXi[k] - dXi[i];
			double dY1 = dYi[k] - dYi[i];
			double dZ1 = dZi[k] - dZi[i];

			double dDist2 = dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1;

			double dQ = sqrt(dDist2) / dH;
			double dTwoMinusQ = 2.0 - dQ;

			double dW;
			double dWq;
			double dWqq;
			if (dQ < 1.0) {
				dW = 2.0/3.0 - dQ * dQ * (1.0 - 0.5 * dQ);
				dWq =  -2.0 * dQ + 1.5 * dQ * dQ;
				dWqq = -2.0 + 3.0 * dQ;
			} else if (dQ < 2.0) {
				dW = 1.0/6.0 * dTwoMinusQ * dTwoMinusQ * dTwoMinusQ;
				dWq = -0.5 * dTwoMinusQ * dTwoMinusQ;
				dWqq = dTwoMinusQ;
			} else {
				dW = 0.0;
				dWq = 0.0;
				dWqq = 0.0;
			}

			dW *= dHScale;
			dWq *= dHScale;
			dWqq *= dHScale;

			dVolume += dW;

			vecCols.push_back(k);
			vecWeight.push_back(dWqq / (dH * dH) + dWq * 2.0 / dH);
		}
		if (dVolume <= 0.0) {
			_EXCEPTIONT("Logic error: point with negative volume detected");
		}
		if (setPoints.size() < 5) {
			Announce("WARNING: Fewer than 5 points used for Laplacian in cell %i"
				" -- accuracy may be affected", i);
		}

		for (int j = 0; j < vecWeight.size(); j++) {
			opLaplacian(i, vecCols[j]) = vecWeight[j] / dVolume;
		}
	}

	kd_free(kdGrid);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Grid connectivity file
	std::string strConnectivity;

	// Input integrated water vapor variable name
	std::string strIWVVariable;

	// Output variable name
	std::string strOutputVariable;

	// Number of points in Laplacian calculation
	int nLaplacianPoints;

	// Radius of the Laplacian
	double dLaplacianDist;

	// Laplacian smooth distance
	double dLaplacianSmoothDistFactor;

	// Minimum Laplacian
	double dMinLaplacian;

	// Minimum absolute latitude
	double dMinAbsLat;

	// Maximum absolute latitude of equatorial band
	double dEqBandMaxLat;

	// Minimum pointwise integrated water vapor
	double dMinIWV;

	// Minimum area
	double dMinArea;

	// Maximal areal fraction
	double dMaxArealFraction;

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

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		//CommandLineString(strInputFileList, "inlist", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strIWVVariable, "var", "");
		CommandLineString(strOutputVariable, "outvar", "");
		CommandLineInt(nLaplacianPoints, "laplacianpoints", 8);
		CommandLineDoubleD(dLaplacianDist, "laplaciandist", 10.0, "(degrees)");
		CommandLineDouble(dLaplacianSmoothDistFactor, "laplaciansmoothdist", 1.2);
		//CommandLineDoubleD(dMinLaplacianSize, "minlaplaciansize", 8.0, "(degrees)");
		//CommandLineDoubleD(dMaxLaplacianSize, "maxlaplaciansize", 12.0, "(degrees)");
		CommandLineDouble(dMinLaplacian, "minlaplacian", 0.5e4);
		CommandLineDouble(dMinAbsLat, "minabslat", 15.0);
		CommandLineDouble(dEqBandMaxLat, "eqbandmaxlat", 15.0);
		CommandLineDouble(dMinIWV, "minval", 20.0);
		CommandLineInt(nMinArea, "minarea", 0);
		CommandLineInt(nAddTimeDim, "addtimedim", -1);
		CommandLineString(strAddTimeDimUnits, "addtimedimunits", "");
		CommandLineBool(fOutputLaplacian, "laplacianout");
		CommandLineBool(fRegional, "regional");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if (strInputFile == "") {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check variable
	if (strIWVVariable == "") {
		_EXCEPTIONT("No IWV variable name (--iwvvar) specified");
	}

	// Check output variable
	if (strOutputVariable.length() == 0) {
		strOutputVariable = strIWVVariable + "tag";
	}

	// Open the NetCDF input file
	NcFile ncInput(strInputFile.c_str());

	if (!ncInput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for reading",
			strInputFile.c_str());
	}

	AnnounceStartBlock("Building grid");

	// Define SimpleGrid
	SimpleGrid grid;

	// Dimensions
	int nSize = 0;
	int nLon = 0;
	int nLat = 0;

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

		nSize = grid.GetSize();

	// No connectivity file; check for latitude/longitude dimension
	} else {
		// Get the longitude dimension
		NcDim * dimLon = ncInput.get_dim("lon");
		if (dimLon == NULL) {
			_EXCEPTIONT("Error accessing dimension \"lon\"");
		}

		// Get the longitude variable
		NcVar * varLon = ncInput.get_var("lon");
		if (varLon == NULL) {
			_EXCEPTIONT("Error accessing variable \"lon\"");
		}

		DataVector<double> dLonDeg(dimLon->size());
		varLon->get(&(dLonDeg[0]), dimLon->size());

		// Get the latitude dimension
		NcDim * dimLat = ncInput.get_dim("lat");
		if (dimLat == NULL) {
			_EXCEPTIONT("Error accessing dimension \"lat\"");
		}

		// Get the latitude variable
		NcVar * varLat = ncInput.get_var("lat");
		if (varLat == NULL) {
			_EXCEPTIONT("Error accessing variable \"lat\"");
		}

		DataVector<double> dLatDeg(dimLat->size());
		varLat->get(&(dLatDeg[0]), dimLat->size());

		// Generate the SimpleGrid
		grid.GenerateLatitudeLongitude(dLatDeg, dLonDeg, fRegional);
	}

	AnnounceEndBlock("Done");

	AnnounceStartBlock("Copying metadata to output file");

	// Get the time dimension
	NcDim * dimTime = ncInput.get_dim("time");
	//if (dimTime == NULL) {
	//	_EXCEPTIONT("Error accessing dimension \"time\"");
	//}

	// Get the integrated water vapor variable
	NcVar * varIWV = ncInput.get_var(strIWVVariable.c_str());
	if (varIWV == NULL) {
		_EXCEPTION1("Error accessing variable \"%s\"",
			strIWVVariable.c_str());
	}

	// Open the NetCDF output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Copy over latitude, longitude and time variables to output file
	NcDim * dimTimeOut = NULL;
	if (dimTime != NULL) {
		CopyNcVar(ncInput, ncOutput, "time", true);
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

	// FIX for unstructured grids
	CopyNcVar(ncInput, ncOutput, "lat", true);
	CopyNcVar(ncInput, ncOutput, "lon", true);

	NcDim * dim0 = ncOutput.get_dim("lat");
	if (dim0 == NULL) {
		_EXCEPTIONT("Error copying variable \"lat\" to output file");
	}
	NcDim * dim1 = ncOutput.get_dim("lon");
	if (dim1 == NULL) {
		_EXCEPTIONT("Error copying variable \"lon\" to output file");
	}

	AnnounceEndBlock("Done");

	// Build Laplacian operator
	AnnounceStartBlock("Building Laplacian operator");

	SparseMatrix<double> opLaplacian;
	BuildLaplacianOperator(
		grid,
		nLaplacianPoints,
		dLaplacianDist,
		dLaplacianSmoothDistFactor,
		opLaplacian);

	AnnounceEndBlock("Done");

	// Input buffers and output variables
	DataVector<double> dIWV(grid.GetSize());

	DataVector<int> dIWVtag(grid.GetSize());
	DataVector<double> dLaplacian(grid.GetSize());

	// Create output variables
	NcVar * varIWVtag = NULL;
	NcVar * varLaplacian = NULL;

	if (grid.m_nGridDim.size() == 1) {
		if (dimTimeOut != NULL) {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dimTimeOut,
				dim0);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dimTimeOut,
				dim0);

		} else {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dim0);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dim0);
		}

	} else if (grid.m_nGridDim.size() == 2) {
		if (dimTimeOut != NULL) {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dimTimeOut,
				dim0,
				dim1);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dimTimeOut,
				dim0,
				dim1);

		} else {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dim0,
				dim1);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dim0,
				dim1);
		}

	} else {
		_EXCEPTIONT("Invalid grid dimension -- value must be 1 or 2");
	}

	// Loop through all times
	int nTimes = 1;
	if (dimTime != NULL) {
		nTimes = dimTime->size();
	}

	for (int t = 0; t < nTimes; t++) {

		// Announce
		char szBuffer[20];
		sprintf(szBuffer, "Time %i", t);
		AnnounceStartBlock(szBuffer);

		// Load in data at this time slice
		AnnounceStartBlock("Reading data");
		if (dimTime != NULL) {
			if (grid.m_nGridDim.size() == 1) {
				varIWV->set_cur(t, 0);
				varIWV->get(&(dIWV[0]), 1, grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varIWV->set_cur(t, 0, 0);
				varIWV->get(&(dIWV[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}

		} else {
			if (grid.m_nGridDim.size() == 1) {
				varIWV->set_cur((long)0);
				varIWV->get(&(dIWV[0]), grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varIWV->set_cur(0, 0);
				varIWV->get(&(dIWV[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}
		}
		AnnounceEndBlock("Done");

		// Compute Laplacian
		AnnounceStartBlock("Applying Laplacian");
		opLaplacian.Apply(dIWV, dLaplacian);
		AnnounceEndBlock("Done");

/*

		// Compute zonal threshold
		DataVector<float> dZonalThreshold(dimLat->size());
		for (int j = 0; j < dimLat->size(); j++) {
			float dMaxZonalIWV = dIWV[j][0];
			for (int i = 0; i < dimLon->size(); i++) {
				dZonalThreshold[j] += dIWV[j][i];
				if (dIWV[j][i] > dMaxZonalIWV) {
					dMaxZonalIWV = dIWV[j][i];
				}
			}
			dZonalThreshold[j] /= static_cast<float>(dimLon->size());

			dZonalThreshold[j] =
				dZonalMeanWeight * dZonalThreshold[j]
				+ dZonalMaxWeight * dMaxZonalIWV;
		}

		// Compute meridional threshold
		DataVector<float> dMeridThreshold(dimLon->size());
		for (int i = 0; i < dimLon->size(); i++) {
			float dMaxMeridIWV = dIWV[0][i];
			for (int j = 0; j < dimLat->size(); j++) {
				dMeridThreshold[i] += dIWV[j][i];
				if (dIWV[j][i] > dMaxMeridIWV) {
					dMaxMeridIWV = dIWV[j][i];
				}
			}
			dMeridThreshold[i] /= static_cast<float>(dimLon->size());

			dMeridThreshold[i] =
				dMeridMeanWeight * dMeridThreshold[i]
				+ dMeridMaxWeight * dMaxMeridIWV;
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
			if (dIWV[j][i] < dMinIWV) {
				continue;
			}
			if (dIWV[j][i] < dZonalThreshold[j]) {
				continue;
			}
			if (dIWV[j][i] < dMeridThreshold[i]) {
				continue;
			}

			//dIWVtag[j][i] = 1 + static_cast<int>(dLaplacian[j][i]);
			if (dLaplacian[j][i] > -dMinLaplacian) {
				dIWVtag[j][i] = 0;
			} else {
				dIWVtag[j][i] = 1;
			}
		}
		}

		// Only retain blobs with minimum area
		if (nMinArea > 0) {
			std::set<Point> setBlobs;

			for (int j = 0; j < dimLat->size(); j++) {
			for (int i = 0; i < dimLon->size(); i++) {
				if (dIWVtag[j][i] != 0) {
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
						dIWVtag[iter->first][iter->second] = 0;
					}
				}
			}
		}
*/
/*
		// Remove points connected with equatorial moisture band
		for (int i = 0; i < dimLon->size(); i++) {
			bool fSouthDone = false;
			bool fNorthDone = false;

			for (int j = 0; j < dimLat->size(); j++) {
				if ((!fSouthDone) && (dLatDeg[j] > -dMinAbsLat)) {
					for (int k = j-1; k > 0; k--) {
						if (dIWVtag[k][i] == 0) {
							break;
						}
						if (dLatDeg[k] < -dEqBandMaxLat) {
							break;
						}
						dIWVtag[k][i] = 0;
					}
					fSouthDone = true;
				}
				if ((!fNorthDone) && (dLatDeg[j] > dMinAbsLat)) {
					for (int k = j-1; k < dimLat->size(); k++) {
						if (dIWVtag[k][i] == 0) {
							break;
						}
						if (dLatDeg[k] > dEqBandMaxLat) {
							break;
						}
						dIWVtag[k][i] = 0;
					}
					fNorthDone = true;
				}
			}
		}
*/
		AnnounceStartBlock("Writing results");

		// Output tagged cell array
		if (dimTimeOut != NULL) {
			if (varLaplacian != NULL) {
				if (grid.m_nGridDim.size() == 1) {
					varLaplacian->set_cur(t, 0);
					varLaplacian->put(&(dLaplacian[0]), 1, grid.m_nGridDim[0]);
				} else if (grid.m_nGridDim.size() == 2) {
					varLaplacian->set_cur(t, 0, 0);
					varLaplacian->put(&(dLaplacian[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
				} else {
					_EXCEPTION();
				}
			}

			//varIWVtag->set_cur(t, 0, 0);
			//varIWVtag->put(&(dIWVtag[0][0]), 1, dimLatOut->size(), dimLonOut->size());

		} else {
			if (varLaplacian != NULL) {
				if (grid.m_nGridDim.size() == 1) {
					varLaplacian->set_cur((long)0);
					varLaplacian->put(&(dLaplacian[0]), grid.m_nGridDim[0]);
				} else if (grid.m_nGridDim.size() == 2) {
					varLaplacian->set_cur(0, 0);
					varLaplacian->put(&(dLaplacian[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
				} else {
					_EXCEPTION();
				}
			}

			//varIWVtag->set_cur(0, 0);
			//varIWVtag->put(&(dIWVtag[0][0]), dimLatOut->size(), dimLonOut->size());
		}

		AnnounceEndBlock("Done");

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

