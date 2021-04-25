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
#include "DataArray1D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#include "SimpleGrid.h"

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

void PrepareBlobOutputVar(
	NcFile & ncInput,
	NcFile & ncOutput,
	const std::string & strOutputFile,
	const SimpleGrid & grid,
	const std::string & strTagVar,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName,
	NcType nctype,
	NcDim * dimTime,
	NcDim ** pdim0,
	NcDim ** pdim1,
	NcVar ** pvarTag
) {
	_ASSERT(pdim0 != NULL);
	_ASSERT(pdim1 != NULL);
	_ASSERT(pvarTag != NULL);

	NcDim * dim0 = NULL;
	NcDim * dim1 = NULL;
	NcVar * varTag = NULL;

	// Copy latitude and longitude arrays
	CopyNcVarIfExists(ncInput, ncOutput, strLatitudeName, true);
	CopyNcVarIfExists(ncInput, ncOutput, strLongitudeName, true);

	// Allocate tag variable (unstructured grids)
	if (grid.m_nGridDim.size() == 1) {
		dim0 = ncOutput.get_dim("ncol");
		if (dim0 == NULL) {
			dim0 = ncOutput.add_dim("ncol", grid.GetSize());
		}
		if (dim0 == NULL) {
			_EXCEPTION1("Error creating dim \"ncol\" in file \"%s\"",
				strOutputFile.c_str());
		}
		if (dim0->size() != grid.GetSize()) {
			_EXCEPTION4("Dimension \"%s\" in file \"%s\" has inconsistent length (%i vs %i)",
				dim0->name(),
				strOutputFile.c_str(),
				dim0->size(),
				grid.GetSize());
		}

		if (dimTime != NULL) {
			varTag = ncOutput.add_var(
				strTagVar.c_str(),
				nctype,
				dimTime,
				dim0);

		} else {
			varTag = ncOutput.add_var(
				strTagVar.c_str(),
				nctype,
				dim0);
		}

	// Allocate tag variable (latitude-longitude grids)
	} else if (grid.m_nGridDim.size() == 2) {

		// Copy over latitude and longitude variable
		dim0 = ncOutput.get_dim(strLatitudeName.c_str());
		if (dim0 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				strLatitudeName.c_str());
		}
		dim1 = ncOutput.get_dim(strLongitudeName.c_str());
		if (dim1 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				strLongitudeName.c_str());
		}

		// Create output tag
		if (dimTime != NULL) {
			varTag = ncOutput.add_var(
				strTagVar.c_str(),
				nctype,
				dimTime,
				dim0,
				dim1);

		} else {
			varTag = ncOutput.add_var(
				strTagVar.c_str(),
				nctype,
				dim0,
				dim1);
		}

	} else {
		_EXCEPTIONT("Invalid grid dimension -- value must be 1 or 2");
	}

	_ASSERT(varTag != NULL);

	(*pdim0) = dim0;
	(*pdim1) = dim1;
	(*pvarTag) = varTag;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A structure for storing a bounding box in latitude / longitude space.
///		The box is treated as singularly periodic in the longitudinal direction
///		if lon_periodic is set to true.
///	</summary>
template <typename Type>
class LatLonBox {

public:
	///	<summary>
	///		Flag indicating this is a null box.
	///	</summary>
	bool is_null;

	///	<summary>
	///		Flag indicating this is a regional.
	///	</summary>
	bool lon_periodic;

	///	<summary>
	///		Width of the longitude variable.
	///	</summary>
	Type lon_width;

	///	<summary>
	///		Bounding longitudes (endpoints are included).
	///	</summary>
	Type lon[2];

	///	<summary>
	///		Bounding latitudes (endpoints are included).
	///	</summar>
	Type lat[2];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	LatLonBox(
		bool a_lon_periodic = true,
		Type a_lon_width = static_cast<Type>(360)
	) :
		is_null(true),
		lon_periodic(a_lon_periodic),
		lon_width(a_lon_width)
	{
		lon[0] = static_cast<Type>(0);
		lon[1] = static_cast<Type>(0);
		lat[0] = static_cast<Type>(0);
		lat[1] = static_cast<Type>(0);
	}

	///	<summary>
	///		Width of this LatLonBox.
	///	</summary>
	Type width() const {
		if (is_null) {
			return 0;
		}

		if (!lon_periodic) {
			return (lon[1] - lon[0]);
		}

		if (lon[0] == lon[1]) {
			return static_cast<Type>(0);
		} else if (lon[0] <= lon[1]) {
			return (lon[1] - lon[0]);
		} else {
			return (lon[1] - lon[0] + lon_width);
		}
	}

	///	<summary>
	///		Insert a point into this LatLonBox.
	///	</summary>
	void insert(
		const Type & lat_pt,
		const Type & lon_pt
	) {
		if (is_null) {
			is_null = false;
			lat[0] = lat_pt;
			lat[1] = lat_pt;
			lon[0] = lon_pt;
			lon[1] = lon_pt;
			return;
		}
		if (lon_pt < static_cast<Type>(0)) {
			std::stringstream strError;
			strError << "lon_pt out of range (" << lon_pt << " < 0)" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}
		if (lon_pt > lon_width) {
			std::stringstream strError;
			strError << "lon_pt out of range (" << lon_pt << " > ";
			strError << lon_width << ")" << std::endl;
			_EXCEPTIONT(strError.str().c_str());
		}

		// Expand latitudes
		if (lat_pt > lat[1]) {
			lat[1] = lat_pt;
		}
		if (lat_pt < lat[0]) {
			lat[0] = lat_pt;
		}

		// Expand longitude, if non-periodic
		if (!lon_periodic) {
			if (lon_pt > lon[1]) {
				lon[1] = lon_pt;
			}
			if (lon_pt < lon[0]) {
				lon[0] = lon_pt;
			}
			return;
		}

		// New longitude lies within existing range
		if (lon[0] <= lon[1]) {
			if ((lon_pt >= lon[0]) && (lon_pt <= lon[1])) {
				return;
			}
		} else {
			if ((lon_pt >= lon[0]) || (lon_pt <= lon[1])) {
				return;
			}
		}

		// New longitude lies outside of existing range
		LatLonBox boxA(*this);
		boxA.lon[0] = lon_pt;

		LatLonBox boxB(*this);
		boxB.lon[1] = lon_pt;

		// The updated box is the box of minimum width
		Type dWidthNow = width();
		Type dWidthA = boxA.width();
		Type dWidthB = boxB.width();

		if ((dWidthA - dWidthNow < -1.0e-14) || (dWidthB - dWidthNow < -1.0e-14)) {
			_EXCEPTIONT("Logic error");
		}

		if (dWidthA < dWidthB) {
			(*this) = boxA;
		} else {
			(*this) = boxB;
		}
	}

	///	<summary>
	///		Determine if this LatLonBox contains the given point.
	///	</summary>
	bool contains(
		const Type & lat_pt,
		const Type & lon_pt
	) {
		// Check latitudes
		if (lat[0] > lat_pt) {
			return false;
		}
		if (lat[1] < lat_pt) {
			return false;
		}

		// Check longitudes, if non-periodic
		if (!lon_periodic) {
			if (lon[0] > lon_pt) {
				return false;
			}
			if (lon[1] < lon_pt) {
				return false;
			}
			return true;
		}

		// This box crosses lon 360
		if (lon[0] > lon[1]) {
			if (lon_pt >= lon[0]) {
				return true;
			}
			if (lon_pt <= lon[1]) {
				return true;
			}
			return false;
		}

		// This box does not cross lon 360
		if ((lon_pt >= lon[0]) && (lon_pt <= lon[1])) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Determine if this LatLonBox overlaps with box.
	///	</summary>
	bool overlaps(
		const LatLonBox<Type> & box
	) const {

		// Check flags
		if ((is_null) || (box.is_null)) {
			return false;
		}
		if (lon_periodic != box.lon_periodic) {
			_EXCEPTION2("Inconsistent periodicity flag (%s/%s)",
				(lon_periodic)?("true"):("false"),
				(box.lon_periodic)?("true"):("false"));
		}
		if (lon_width != box.lon_width) {
			_EXCEPTIONT("Inconsistent box lon_width in comparison");
		}

		// Check latitudes
		if (lat[0] > box.lat[1]) {
			return false;
		}
		if (lat[1] < box.lat[0]) {
			return false;
		}

		// Cases when longitude is periodic
		if (lon_periodic) {

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
///		Get a DataArray1D containing the time variable across a list of
///		input files.
///	</summary>
void GetAllTimes(
	const std::vector<std::string> & vecInputFiles,
	DataArray1D<double> & dataTimes,
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
			DataArray1D<double> dTime(nTime);

			varTime->get(dTime, nTime);
			for (int t = 0; t < nTime; t++) {
				vecTimes.push_back(dTime[t]);
			//vecTimes.push_back(dTime[t] + (double)(1440 * f));
			}
		}
	}

	dataTimes.Allocate(vecTimes.size());
	memcpy(&(dataTimes[0]), &(vecTimes[0]), vecTimes.size() * sizeof(double));
}

///////////////////////////////////////////////////////////////////////////////

#endif

