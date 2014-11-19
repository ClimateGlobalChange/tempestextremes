///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectCyclones.cpp
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
#include "DataVector.h"
#include "DataMatrix.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all minima in the given DataMatrix.
///	</summary>
void FindAllLocalMinima(
	const DataMatrix<float> data,
	std::set< std::pair<int, int> > & setMinima
) {
	int nLon = data.GetColumns();
	int nLat = data.GetRows();

	for (int j = 1; j < nLat-1; j++) {
	for (int i = 0; i < nLon; i++) {

		bool fMinimum = true;
		for (int q = -1; q <= 1; q++) {
		for (int p = -1; p <= 1; p++) {
			int ix = (i + nLon + p) % nLon;
			int jx = (j + q);

			if (data[jx][ix] < data[j][i]) {
				fMinimum = false;
				goto DonePressureMinima;
			}
		}
		}

DonePressureMinima:
		if (fMinimum) {
			setMinima.insert(std::pair<int,int>(j,i));
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataMatrix.
///	</summary>
void FindAllLocalMaxima(
	const DataMatrix<float> data,
	std::set< std::pair<int, int> > & setMaxima
) {
	int nLon = data.GetColumns();
	int nLat = data.GetRows();

	for (int j = 1; j < nLat-1; j++) {
	for (int i = 0; i < nLon; i++) {

		bool fMaximum = true;
		for (int q = -1; q <= 1; q++) {
		for (int p = -1; p <= 1; p++) {
			int ix = (i + nLon + p) % nLon;
			int jx = (j + q);

			if (data[jx][ix] > data[j][i]) {
				fMaximum = false;
				goto DonePressureMaxima;
			}
		}
		}

DonePressureMaxima:
		if (fMaximum) {
			setMaxima.insert(std::pair<int,int>(j,i));
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Require temperature maxima at T200 and T500 within this distance
	double dWarmCoreDist;

	// No temperature maxima at T200 and T500 within this distance
	double dNoWarmCoreDist;

	// Minimum Laplacian value
	double dMinLaplacian;

	// Distance to search for maximum wind speed
	double dWindSpDist;

	// Append to output file
	bool fAppend;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineDoubleD(dWarmCoreDist, "warmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dNoWarmCoreDist, "nowarmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLaplacian, "minlaplacian", 0.0, "(Pa / degree^2)");
		CommandLineDoubleD(dWindSpDist, "windspdist", 0.0, "(degrees)");
		CommandLineBool(fAppend, "append");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check input
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input file (--in) specified");
	}

	// Check output
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Check warm core distance
	if ((dWarmCoreDist != 0.0) && (dNoWarmCoreDist != 0.0)) {
		_EXCEPTIONT("Only one of --warmcoredist and --nowarmcoredist"
			   " may be active");
	}

	// Load the netcdf file
	NcFile ncInput(strInputFile.c_str());

	// Get latitude/longitude dimensions
	NcDim * dimLat = ncInput.get_dim("lat");
	NcDim * dimLon = ncInput.get_dim("lon");

	NcVar * varLat = ncInput.get_var("lat");
	NcVar * varLon = ncInput.get_var("lon");

	int nLat = dimLat->size();
	int nLon = dimLon->size();

	DataVector<double> dataLat(nLat);
	varLat->get(dataLat, nLat);

	for (int j = 0; j < nLat; j++) {
		dataLat[j] *= M_PI / 180.0;
	}

	DataVector<double> dataLon(nLon);
	varLon->get(dataLon, nLon);

	for (int i = 0; i < nLon; i++) {
		dataLon[i] *= M_PI / 180.0;
	}

	// Get time dimension
	NcDim * dimTime = ncInput.get_dim("time");
	NcVar * varTime = ncInput.get_var("time");

	int nTime = dimTime->size();

	DataVector<double> dTime;
	dTime.Initialize(nTime);

	varTime->get(dTime, nTime);

	// Get auxiliary variables
	NcVar * varPSL  = ncInput.get_var("PSL");
	NcVar * varU850 = ncInput.get_var("U850");
	NcVar * varV850 = ncInput.get_var("V850");
	NcVar * varT200 = ncInput.get_var("T200");
	NcVar * varT500 = ncInput.get_var("T500");

	// Storage for auxiliary variables
	DataMatrix<float> dataPSL(nLat, nLon);
	DataMatrix<float> dataU850(nLat, nLon);
	DataMatrix<float> dataV850(nLat, nLon);
	DataMatrix<float> dataT200(nLat, nLon);
	DataMatrix<float> dataT500(nLat, nLon);

	DataMatrix<float> dataDel2PSL(nLat, nLon);

	// Loop through all times
	for (int t = 0; t < nTime; t++) {

		char szStartBlock[128];
		sprintf(szStartBlock, "Time %i", t);
		AnnounceStartBlock(szStartBlock);

		// Get the auxiliary variables
		varPSL->set_cur(t,0,0);
		varPSL->get(&(dataPSL[0][0]), 1, nLat, nLon);

		varU850->set_cur(t,0,0);
		varU850->get(&(dataU850[0][0]), 1, nLat, nLon);

		varV850->set_cur(t,0,0);
		varV850->get(&(dataV850[0][0]), 1, nLat, nLon);

		varT200->set_cur(t,0,0);
		varT200->get(&(dataT200[0][0]), 1, nLat, nLon);

		varT500->set_cur(t,0,0);
		varT500->get(&(dataT500[0][0]), 1, nLat, nLon);

		// Tag all pressure minima
		std::set< std::pair<int, int> > setPressureMinima;
		FindAllLocalMinima(dataPSL, setPressureMinima);

		// Total number of pressure minima
		int nTotalPressureMinima = setPressureMinima.size();
		int nRejectedWarmCore = 0;
		int nRejectedNoWarmCore = 0;
		int nRejectedLaplacian = 0;

		// Detect presence of warm core near PSL min
		if ((dWarmCoreDist != 0.0) || (dNoWarmCoreDist != 0.0)) {
			
			std::set< std::pair<int, int> > setT200Maxima;
			FindAllLocalMaxima(dataT200, setT200Maxima);

			std::set< std::pair<int, int> > setT500Maxima;
			FindAllLocalMaxima(dataT500, setT500Maxima);

			// Construct KD tree for T200
			kdtree * kdT200Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterT200
				= setT200Maxima.begin();

			for (; iterT200 != setT200Maxima.end(); iterT200++) {
				double dLat = dataLat[iterT200->first];
				double dLon = dataLon[iterT200->second];

				double dX = sin(dLon) * cos(dLat);
				double dY = cos(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdT200Maxima, dX, dY, dZ, NULL);
			}

			// Construct KD tree for T500
			kdtree * kdT500Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterT500
				= setT500Maxima.begin();

			for (; iterT500 != setT500Maxima.end(); iterT500++) {
				double dLat = dataLat[iterT500->first];
				double dLon = dataLon[iterT500->second];

				double dX = sin(dLon) * cos(dLat);
				double dY = cos(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdT500Maxima, dX, dY, dZ, NULL);
			}

			// Remove pressure minima that are near temperature maxima
			std::set< std::pair<int, int> > setNewPressureMinima;

			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {
				double dLat = dataLat[iterPSL->first];
				double dLon = dataLon[iterPSL->second];

				double dX = sin(dLon) * cos(dLat);
				double dY = cos(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kdres * kdresT200 = kd_nearest3(kdT200Maxima, dX, dY, dZ);
				kdres * kdresT500 = kd_nearest3(kdT500Maxima, dX, dY, dZ);

				double dT200pos[3];
				double dT500pos[3];

				kd_res_item(kdresT200, dT200pos);
				kd_res_item(kdresT500, dT500pos);

				kd_res_free(kdresT200);
				kd_res_free(kdresT500);

				double dT200dist = sqrt(
					  (dX - dT200pos[0]) * (dX - dT200pos[0])
					+ (dY - dT200pos[1]) * (dY - dT200pos[1])
					+ (dZ - dT200pos[2]) * (dZ - dT200pos[2]));

				dT200dist = 2.0 * asin(0.5 * dT200dist) * 180.0 / M_PI;

				double dT500dist = sqrt(
					  (dX - dT500pos[0]) * (dX - dT500pos[0])
					+ (dY - dT500pos[1]) * (dY - dT500pos[1])
					+ (dZ - dT500pos[2]) * (dZ - dT500pos[2]));

				dT500dist = 2.0 * asin(0.5 * dT500dist) * 180.0 / M_PI;

				// Reject storms with warm core
				if (dNoWarmCoreDist != 0.0) {
					if ((dT200dist >= dNoWarmCoreDist) ||
						(dT500dist >= dNoWarmCoreDist)
					) {
						setNewPressureMinima.insert(*iterPSL);
					} else {
						nRejectedWarmCore++;
					}
				}

				// Reject storms with no warm core
				if (dWarmCoreDist != 0.0) {
					if ((dT200dist <= dWarmCoreDist) &&
						(dT500dist <= dWarmCoreDist)
					) {
						setNewPressureMinima.insert(*iterPSL);
					} else {
						nRejectedNoWarmCore++;
					}
				}
			}

			kd_free(kdT200Maxima);
			kd_free(kdT500Maxima);

			setPressureMinima = setNewPressureMinima;
		}

		// Calculate the Laplacian of pressure at the PSL min
		if (dMinLaplacian != 0.0) {
			float dDeltaLat = static_cast<float>(dataLat[1] - dataLat[0]);
			float dDeltaLon = static_cast<float>(dataLon[1] - dataLon[0]);

			std::set< std::pair<int, int> > setNewPressureMinima;

			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {

				int i = iterPSL->second;
				int j = iterPSL->first;

				int inext = (i + 1) % nLon;
				int jnext = (j + 1);

				int ilast = (i + nLon - 1) % nLon;
				int jlast = (j - 1);

				float dDphiPSL =
					(dataPSL[jnext][i] - dataPSL[jlast][i]) / (2.0 * dDeltaLat);
				float dDlambdaPSL =
					(dataPSL[j][inext] - dataPSL[j][ilast]) / (2.0 * dDeltaLon);
				float dD2phiPSL =
					(dataPSL[jnext][i] - 2.0 * dataPSL[j][i] + dataPSL[jlast][i])
					/ (dDeltaLat * dDeltaLat);
				float dD2lambdaPSL =
					(dataPSL[j][inext] - 2.0 * dataPSL[j][i] + dataPSL[j][ilast])
					/ (dDeltaLon * dDeltaLon);

				float dSecLat = 1.0 / cos(dataLat[j]);
			
				float dLaplacian =
					dD2phiPSL - tan(dataLat[j]) * dDphiPSL
					+ dSecLat * dSecLat * dD2lambdaPSL;

				// Convert to Pa / degree^2
				dLaplacian *= (M_PI / 180.0) * (M_PI / 180.0);

				if (dLaplacian >= dMinLaplacian) {
					setNewPressureMinima.insert(*iterPSL);
				} else {
					nRejectedLaplacian++;
				}
			}

			setPressureMinima = setNewPressureMinima;
		}

		Announce("Total candidates: %i", setPressureMinima.size());
		Announce("Rejected (   warm core): %i", nRejectedWarmCore);
		Announce("Rejected (no warm core): %i", nRejectedNoWarmCore);
		Announce("Rejected (   laplacian): %i", nRejectedLaplacian);
		AnnounceEndBlock("Done");
	}

	ncInput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

