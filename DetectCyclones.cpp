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
#include "TimeObj.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all minima in the given DataMatrix.
///	</summary>
void FindAllLocalMinima(
	const DataMatrix<float> & data,
	std::set< std::pair<int, int> > & setMinima
) {
	int nLon = data.GetColumns();
	int nLat = data.GetRows();

	// Check interior nodes
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
/*
	// Check south pole (one node at south pole)
	{
		bool fMinimum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[1][i] < data[0][i]) {
				fMinimum = false;
				break;
			}
		}
		if (fMinimum) {
			_EXCEPTION();
			setMinima.insert(std::pair<int,int>(0,0));
		}
	}

	// Check north pole (one node at north pole)
	{
		bool fMinimum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[nLat-2][i] < data[nLat-1][i]) {
				fMinimum = false;
				break;
			}
		}
		if (fMinimum) {
			_EXCEPTION();
			setMinima.insert(std::pair<int,int>(nLat-1,0));
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataMatrix.
///	</summary>
void FindAllLocalMaxima(
	const DataMatrix<float> & data,
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
/*
	// Check south pole (one node at south pole)
	{
		bool fMaximum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[1][i] > data[0][i]) {
				fMaximum = false;
				break;
			}
		}
		if (fMaximum) {
			_EXCEPTION();
			setMaxima.insert(std::pair<int,int>(0,0));
		}
	}

	// Check north pole (one node at north pole)
	{
		bool fMaximum = true;
		for (int i = 0; i < nLon; i++) {
			if (data[nLat-2][i] > data[nLat-1][i]) {
				fMaximum = false;
				break;
			}
		}
		if (fMaximum) {
			_EXCEPTION();
			setMaxima.insert(std::pair<int,int>(nLat-1,0));
		}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
void FindLocalMaximum(
	const DataMatrix<float> & data,
	const DataVector<double> & dataLon,
	const DataVector<double> & dataLat,
	int iLat,
	int iLon,
	double dMaxDist,
	std::pair<int, int> & prMaximum,
	float & dMaxValue,
	float & dRMax
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Number of latitudes and longitudes
	const int nLat = dataLat.GetRows();
	const int nLon = dataLon.GetRows();

	// Initialize the maximum to the central location
	prMaximum = std::pair<int, int>(iLat, iLon);
	dMaxValue = data[iLat][iLon];
	dRMax = 0.0;

	// Queue of nodes that remain to be visited
	std::queue< std::pair<int, int> > queueNodes;
	queueNodes.push(prMaximum);

	// Set of nodes that have already been visited
	std::set< std::pair<int, int> > setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = dataLat[iLat];
	double dLon0 = dataLon[iLon];

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

		// Check for new local maximum
		if (data[pr.first][pr.second] > dMaxValue) {
			prMaximum = pr;
			dMaxValue = data[pr.first][pr.second];
			dRMax = dR;
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
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a pair of Date values.
///	</summary>
void ParseDate(
	int nDate,
	int nDateSec,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	nDateYear  = nDate / 10000;
	nDateMonth = (nDate % 10000) / 100;
	nDateDay   = (nDate % 100);
	nDateHour  = nDateSec / 3600;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a time value.
///	</summary>
void ParseTimeDouble(
	const std::string & strTimeUnits,
	const std::string & strTimeCalendar,
	double dTime,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	// Get calendar type
	Time::CalendarType cal;
	if (strTimeCalendar == "noleap") {
		cal = Time::CalendarNoLeap;
	} else if (strTimeCalendar == "standard") {
		cal = Time::CalendarStandard;
	} else {
		_EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
	}

	// Time format is "days since ..."
	if ((strTimeUnits.length() >= 11) &&
	    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(11);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dTime);
		time.AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 86400.0);
		time.AddSeconds(nSeconds);

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	// Time format is "hours since ..."
	} else if (
	    (strTimeUnits.length() >= 12) &&
	    (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(12);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 3600.0);
		time.AddSeconds(nSeconds);

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
	//_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {
	// Radius of the Earth
	const float ParamEarthRadius = 6.37122e6;

	// Input file
	std::string strInputFile;

	// Output file
	std::string strOutputFile;

	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// Require temperature maxima at T200 and T500 within this distance
	double dWarmCoreDist;

	// No temperature maxima at T200 and T500 within this distance
	double dNoWarmCoreDist;

	// Require 850hPa vorticity maxima within this distance
	double dVortDist;

	// Minimum Laplacian value
	double dMinLaplacian;

	// Distance to search for maximum wind speed
	double dWindSpDist;

	// Append to output file
	bool fAppend;

	// Output header
	bool fOutputHeader;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineDoubleD(dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineDoubleD(dWarmCoreDist, "warmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dNoWarmCoreDist, "nowarmcoredist", 0.0, "(degrees)");
		CommandLineDoubleD(dVortDist, "vortdist", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLaplacian, "minlaplacian", 0.0, "(Pa / degree^2)");
		CommandLineDoubleD(dWindSpDist, "windspdist", 0.0, "(degrees)");
		//CommandLineBool(fAppend, "append");
		CommandLineBool(fOutputHeader, "out_header");

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

	// Check minimum latitude and maximum latitude
	if (dMaxLatitude < 0.0) {
		_EXCEPTIONT("--maxlat must be nonnegative");
	}
	if (dMinLatitude < 0.0) {
		_EXCEPTIONT("--minlat must be nonnegative");
	}

	dMaxLatitude *= M_PI / 180.0;
	dMinLatitude *= M_PI / 180.0;

	// Load the netcdf file
	NcFile ncInput(strInputFile.c_str());

	// Get latitude/longitude dimensions
	NcDim * dimLat = ncInput.get_dim("lat");
	if (dimLat == NULL) {
		_EXCEPTIONT("No dimension \"lat\" found in input file");
	}

	NcDim * dimLon = ncInput.get_dim("lon");
	if (dimLon == NULL) {
		_EXCEPTIONT("No dimension \"lon\" found in input file");
	}

	NcVar * varLat = ncInput.get_var("lat");
	if (varLat == NULL) {
		_EXCEPTIONT("No variable \"lat\" found in input file");
	}

	NcVar * varLon = ncInput.get_var("lon");
	if (varLon == NULL) {
		_EXCEPTIONT("No variable \"lon\" found in input file");
	}

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
	if (dimTime == NULL) {
		_EXCEPTIONT("No dimension \"time\" found in input file");
	}

	NcVar * varTime = ncInput.get_var("time");
	if (varTime == NULL) {
		_EXCEPTIONT("No variable \"time\" found in input file");
	}

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

	if (varPSL == NULL) {
		_EXCEPTIONT("No variable \"PSL\" found in input file");
	}
	if (varU850 == NULL) {
		_EXCEPTIONT("No variable \"U850\" found in input file");
	}
	if (varV850 == NULL) {
		_EXCEPTIONT("No variable \"V850\" found in input file");
	}
	if (varT200 == NULL) {
		_EXCEPTIONT("No variable \"T200\" found in input file");
	}
	if (varT500 == NULL) {
		_EXCEPTIONT("No variable \"T500\" found in input file");
	}

	// Storage for auxiliary variables
	DataMatrix<float> dataPSL(nLat, nLon);
	DataMatrix<float> dataU850(nLat, nLon);
	DataMatrix<float> dataV850(nLat, nLon);
	DataMatrix<float> dataT200(nLat, nLon);
	DataMatrix<float> dataT500(nLat, nLon);

	DataMatrix<float> dataUMag850(nLat, nLon);

	DataMatrix<float> dataDel2PSL(nLat, nLon);

	// Open output file
	FILE * fpOutput = fopen(strOutputFile.c_str(), "w");
	if (fpOutput == NULL) {
		_EXCEPTION1("Could not open output file \"%s\"",
			strOutputFile.c_str());
	}

	if (fOutputHeader) {
		fprintf(fpOutput, "#day\tmonth\tyear\tcount\thour\n");
		fprintf(fpOutput, "#\t#\ti\tj\tpsl_lon\tpsl_lat\twind_max\tr_wind_max\tpsl_min\n");
	}

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

		// Compute wind magnitude
		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {
			dataUMag850[j][i] = sqrt(
				  dataU850[j][i] * dataU850[j][i]
				+ dataV850[j][i] * dataV850[j][i]);
		}
		}

		// Tag all pressure minima
		std::set< std::pair<int, int> > setPressureMinima;
		FindAllLocalMinima(dataPSL, setPressureMinima);

		// Total number of pressure minima
		int nTotalPressureMinima = setPressureMinima.size();
		int nRejectedLatitude = 0;
		int nRejectedVortMax = 0;
		int nRejectedWarmCore = 0;
		int nRejectedNoWarmCore = 0;
		int nRejectedLaplacian = 0;

		// Eliminate based on maximum latitude
		if (dMaxLatitude != 0.0) {
			std::set< std::pair<int, int> > setNewPressureMinima;

			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {
				double dLat = dataLat[iterPSL->first];

				if (fabs(dLat) <= dMaxLatitude) {
					setNewPressureMinima.insert(*iterPSL);
				} else {
					nRejectedLatitude++;
				}
			}

			setPressureMinima = setNewPressureMinima;
		}

		// Eliminate based on minimum latitude
		if (dMinLatitude != 0.0) {
			std::set< std::pair<int, int> > setNewPressureMinima;

			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {
				double dLat = dataLat[iterPSL->first];

				if (fabs(dLat) >= dMinLatitude) {
					setNewPressureMinima.insert(*iterPSL);
				} else {
					nRejectedLatitude++;
				}
			}

			setPressureMinima = setNewPressureMinima;
		}

		// Eliminate based on presence of vorticity maximum
		if (dVortDist != 0.0) {

			float dDeltaLat = static_cast<float>(dataLat[1] - dataLat[0]);
			float dDeltaLon = static_cast<float>(dataLon[1] - dataLon[0]);

			DataMatrix<float> dataZETA850(nLat, nLon);

			// Compute vorticity
			for (int j = 1; j < nLat-1; j++) {
			for (int i = 0; i < nLon; i++) {

				int inext = (i + 1) % nLon;
				int jnext = (j + 1);

				int ilast = (i + nLon - 1) % nLon;
				int jlast = (j - 1);

				float dDyU850 = (dataU850[jnext][i] - dataU850[jlast][i])
					/ (2.0 * dDeltaLat);
				float dDxV850 = (dataV850[j][inext] - dataV850[j][ilast])
					/ (2.0 * dDeltaLon);

				dataZETA850[j][i] = dDxV850 - dDyU850
					+ dataU850[j][i] / ParamEarthRadius * tan(dataLat[j]);

				if (dataLat[j] < 0.0) {
					dataZETA850[j][i] *= -1.0;
				}
			}
			}

			// Find all vorticity maxima
			std::set< std::pair<int, int> > setZETA850Maxima;
			FindAllLocalMaxima(dataZETA850, setZETA850Maxima);

			// Construct KD tree for ZETA850
			kdtree * kdZETA850Maxima = kd_create(3);
			std::set< std::pair<int, int> >::const_iterator iterZETA850
				= setZETA850Maxima.begin();

			for (; iterZETA850 != setZETA850Maxima.end(); iterZETA850++) {
				double dLat = dataLat[iterZETA850->first];
				double dLon = dataLon[iterZETA850->second];

				double dX = sin(dLon) * cos(dLat);
				double dY = cos(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdZETA850Maxima, dX, dY, dZ, NULL);
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

				kdres * kdresZETA850 = kd_nearest3(kdZETA850Maxima, dX, dY, dZ);

				double dZETA850pos[3];

				kd_res_item(kdresZETA850, dZETA850pos);

				kd_res_free(kdresZETA850);

				double dZETA850dist = sqrt(
					  (dX - dZETA850pos[0]) * (dX - dZETA850pos[0])
					+ (dY - dZETA850pos[1]) * (dY - dZETA850pos[1])
					+ (dZ - dZETA850pos[2]) * (dZ - dZETA850pos[2]));

				dZETA850dist = 2.0 * asin(0.5 * dZETA850dist) * 180.0 / M_PI;

				// Reject storms with no warm core
				if (dZETA850dist <= dVortDist) {
					setNewPressureMinima.insert(*iterPSL);
				} else {
					nRejectedVortMax++;
				}
			}

			kd_free(kdZETA850Maxima);

			setPressureMinima = setNewPressureMinima;
		}

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
		Announce("Rejected (    latitude): %i", nRejectedLatitude);
		Announce("Rejected (no vort max ): %i", nRejectedVortMax);
		Announce("Rejected (   warm core): %i", nRejectedWarmCore);
		Announce("Rejected (no warm core): %i", nRejectedNoWarmCore);
		Announce("Rejected (   laplacian): %i", nRejectedLaplacian);

		// Determine wind maximum at all pressure minima
		AnnounceStartBlock("Searching for maximum winds");
		std::vector<float> vecMaxWindSp;
		std::vector<float> vecRMaxWindSp;
		{
			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {

				std::pair<int, int> prMaximum;

				float dRMaxWind;
				float dMaxWindSp;

				FindLocalMaximum(
					dataUMag850,
					dataLon,
					dataLat,
					iterPSL->first,
					iterPSL->second,
					dWindSpDist,
					prMaximum,
					dMaxWindSp,
					dRMaxWind);

				vecMaxWindSp.push_back(dMaxWindSp);
				vecRMaxWindSp.push_back(dRMaxWind);
			}
		}

		AnnounceEndBlock("Done");

		// Write results to file
		{
			// Parse time information
			//NcVar * varDate = ncInput.get_var("date");
			//NcVar * varDateSec = ncInput.get_var("datesec");

			int nDateYear;
			int nDateMonth;
			int nDateDay;
			int nDateHour;
/*
			if ((varDate != NULL) && (varDateSec != NULL)) {
				int nDate;
				int nDateSec;

				varDate->set_cur(t);
				varDate->get(&nDate, 1);

				varDateSec->set_cur(t);
				varDateSec->get(&nDateSec, 1);

				ParseDate(
					nDate,
					nDateSec,
					nDateYear,
					nDateMonth,
					nDateDay,
					nDateHour);

			} else {
*/
				NcAtt * attTimeUnits = varTime->get_att("units");
				if (attTimeUnits == NULL) {
					_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
				}

				std::string strTimeUnits = attTimeUnits->as_string(0);


				std::string strTimeCalendar = "noleap";
				NcAtt * attTimeCalendar = varTime->get_att("calendar");
				if (attTimeUnits != NULL) {
					strTimeCalendar = attTimeCalendar->as_string(0);
				}


				ParseTimeDouble(
					strTimeUnits,
					strTimeCalendar,
					dTime[t],
					nDateYear,
					nDateMonth,
					nDateDay,
					nDateHour);
//			}

			// Write time information
			fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
				nDateYear,
				nDateMonth,
				nDateDay,
				static_cast<int>(setPressureMinima.size()),
				nDateHour);

			// Write candidate information
			int iCandidateCount = 0;

			std::set< std::pair<int, int> >::const_iterator iterPSL
				= setPressureMinima.begin();
			for (; iterPSL != setPressureMinima.end(); iterPSL++) {
				fprintf(fpOutput, "%i\t%i\t%i\t%3.6f\t%3.6f\t%2.6f\t%3.6f\t%6.6f\n",
					iCandidateCount,
					iterPSL->second,
					iterPSL->first,
					dataLon[iterPSL->second] * 180.0 / M_PI,
					dataLat[iterPSL->first]  * 180.0 / M_PI,
					vecMaxWindSp[iCandidateCount],
					vecRMaxWindSp[iCandidateCount],
					dataPSL[iterPSL->first][iterPSL->second]);

				iCandidateCount++;
			}
		}

		AnnounceEndBlock("Done");
	}

	fclose(fpOutput);

	ncInput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

