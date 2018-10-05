///////////////////////////////////////////////////////////////////////////////
///
///	\file    BlobUtilities.h
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

#ifndef _BLOBUTILITIES_H_
#define _BLOBUTILITIES_H_

#include "Exception.h"
#include "Announce.h"
#include "DataVector.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <vector>
#include <map>
#include <string>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure for storing a latitude / longitude index.
///	</summary>
struct LatLonPair {

	int lat;
	int lon;

	///	<summary>
	///		Constructor.
	///	</summary>
	LatLonPair(int a_lat, int a_lon) :
		lat(a_lat), lon(a_lon)
	{ }

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const LatLonPair & pr) const {
		if (lat < pr.lat) {
			return true;
		} else if (lat > pr.lat) {
			return false;
		}

		if (lon < pr.lon) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const LatLonPair & pr) const {
		if ((lat == pr.lat) && (lon == pr.lon)) {
			return true;
		} else {
			return false;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure for storing a bounding box in latitude / longitude space.
///	</summary>
class LatLonBox {

public:
	///	<summary>
	///		Flag indicating this is a null box.
	///	</summary>
	bool is_null;

	///	<summary>
	///		Bounding latitudes (endpoints are included).
	///	</summar>
	int lat[2];

	///	<summary>
	///		Bounding longitudes (endpoints are included).
	///	</summary>
	int lon[2];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	LatLonBox() :
		is_null(true)
	{ }

	///	<summary>
	///		Width of this LatLonBox.
	///	</summary>
	int Width(
		int nLonCount
	) const {
		if (is_null) {
			return 0;
		}

		if (lon[0] == lon[1]) {
			return (1);

		} else if (lon[0] <= lon[1]) {
			return (lon[1] - lon[0] + 1);
		} else {
			return (lon[1] - lon[0] + nLonCount + 1);
		}
	}

	///	<summary>
	///		Insert a point into this LatLonBox.
	///	</summary>
	void InsertPoint(
		int iLat,
		int iLon,
		int nLatCount,
		int nLonCount
	) {
		if (is_null) {
			is_null = false;
			lat[0] = iLat;
			lat[1] = iLat;
			lon[0] = iLon;
			lon[1] = iLon;
			return;
		}

		// Expand latitudes
		if (iLat > lat[1]) {
			lat[1] = iLat;
		}
		if (iLat < lat[0]) {
			lat[0] = iLat;
		}

		// New longitude lies within existing range
		if (lon[0] <= lon[1]) {
			if ((iLon >= lon[0]) && (iLon <= lon[1])) {
				return;
			}
		} else {
			if ((iLon >= lon[0]) || (iLon <= lon[1])) {
				return;
			}
		}

		// New longitude lies outside of existing range
		LatLonBox boxA(*this);
		boxA.lon[0] = iLon;

		LatLonBox boxB(*this);
		boxB.lon[1] = iLon;

		int iWidthNow = Width(nLonCount);
		int iWidthA = boxA.Width(nLonCount);
		int iWidthB = boxB.Width(nLonCount);

		if ((iWidthA < iWidthNow) || (iWidthB < iWidthNow)) {
			_EXCEPTIONT("Logic error");
		}

		if (iWidthA < iWidthB) {
			(*this) = boxA;
		} else {
			(*this) = boxB;
		}
/*
		if (iLon == nLonCount-1) {
			if (lon[0] <= lon[1]) {
				if (nLonCount - 1 - lon[1] < lon[0] + 1) {
					lon[1] = nLonCount - 1;
				} else {
					lon[0] = nLonCount - 1;
				}
			}

		} else if (iLon == 0) {
			if (lon[0] <= lon[1]) {
				if (nLonCount - lon[1] < lon[0]) {
					lon[1] = 0;
				} else {
					lon[0] = 0;
				}
			}

		} else {
			if (lon[0] <= lon[1]) {
				if (iLon < lon[0]) {
					lon[0] = iLon;
				}
				if (iLon > lon[1]) {
					lon[1] = iLon;
				}

			} else {
				if ((iLon >= lon[0]) ||
					(iLon <= lon[1])
				) {
				} else if (
					iLon - lon[1] < lon[0] - iLon
				) {
					lon[1] = iLon;
				} else {
					lon[0] = iLon;
				}
			}
		}
*/
	}

	///	<summary>
	///		Determine this LatLonBox overlaps with box.
	///	</summary>
	bool Overlaps(const LatLonBox & box) const {

		// Check latitudes
		if (lat[0] > box.lat[1]) {
			return false;
		}
		if (lat[1] < box.lat[0]) {
			return false;
		}

		// Both boxes cross lon 360
		if ((lon[0] > lon[1]) && (box.lon[0] > box.lon[1])) {
			return true;
		}

		// This box crosses lon 360
		if (lon[0] > lon[1]) {
			if (box.lon[1] >= lon[0]) {
				return true;
			}
			if (box.lon[0] <= lon[1]) {
				return true;
			}
			return false;
		}

		// That box crosses lon 360
		if (box.lon[0] > box.lon[1]) {
			if (lon[1] >= box.lon[0]) {
				return true;
			}
			if (lon[0] <= box.lon[1]) {
				return true;
			}
			return false;
		}

		// No boxes cross lon 360
		if (box.lon[1] < lon[0]) {
			return false;
		}
		if (box.lon[0] > lon[1]) {
			return false;
		}
		return true;
	}
};

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

///	<summary>
///		Get a DataVector containing the time variable across a list of
///		input files.
///	</summary>
void GetAllTimes(
	const std::vector<std::string> & vecInputFiles,
	DataVector<double> & dataTimes,
	std::string timename
) {
	std::vector<double> vecTimes;

	for (int f = 0; f < vecInputFiles.size(); f++) {
		NcFile ncFile(vecInputFiles[f].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"",
				vecInputFiles[f].c_str());
		}

		NcDim * dimTime = ncFile.get_dim(timename.c_str());
		if (dimTime == NULL) {
			_EXCEPTIONT("GetAllTimes: No time dimension found in input file");
		}

		int nTime = dimTime->size();

		NcVar * varTime = ncFile.get_var(timename.c_str());
		if (varTime == NULL) {
			for (int t = 0; t < nTime; t++) {
				vecTimes.push_back(static_cast<double>(t));
			}

		} else {
			DataVector<double> dTime;
			dTime.Initialize(nTime);

			varTime->get(dTime, nTime);
			for (int t = 0; t < nTime; t++) {
				vecTimes.push_back(dTime[t]);
			//vecTimes.push_back(dTime[t] + (double)(1440 * f));
			}
		}
	}

	dataTimes.Initialize(vecTimes.size());
	memcpy(&(dataTimes[0]), &(vecTimes[0]), vecTimes.size() * sizeof(double));
}

///////////////////////////////////////////////////////////////////////////////

#endif

