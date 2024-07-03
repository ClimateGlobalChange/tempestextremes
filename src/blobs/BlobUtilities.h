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
#include "LatLonBox.h"

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

		// Get output dimensions
		NcVar * varLat = ncOutput.get_var(strLatitudeName.c_str());
		if ((varLat != NULL) && (varLat->num_dims() == 2)) {
			dim0 = varLat->get_dim(0);
			dim1 = varLat->get_dim(1);

		} else {
			dim0 = ncOutput.get_dim(strLatitudeName.c_str());
			if (dim0 == NULL) {
				_EXCEPTION1("Unable to copy variable \"%s\" to output file",
					strLatitudeName.c_str());
			}
			dim1 = ncOutput.get_dim(strLongitudeName.c_str());
			if (dim1 == NULL) {
				_EXCEPTION1("Unable to copy variable \"%s\" to output file",
					strLongitudeName.c_str());
			}
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
