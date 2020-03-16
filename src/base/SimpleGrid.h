///////////////////////////////////////////////////////////////////////////////
///
///	\file    SimpleGrid.h
///	\author  Paul Ullrich
///	\version August 30, 2019
///
///	<remarks>
///		Copyright 2019 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _SIMPLEGRID_H_
#define _SIMPLEGRID_H_

#include "DataArray1D.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

class kdtree;
class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing the grid, including coordinates of
///		each data point and graph connectivity of elements.
///	</summary>
class SimpleGrid {

public:
	///	<summary>
	///		A identifying the connectivity file format.
	///	</summary>
	static const char * c_szFileIdentifier;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SimpleGrid() :
		m_kdtree(NULL)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~SimpleGrid();

	///	<summary>
	///		Determine if the SimpleGrid is initialized.
	///	</summary>
	bool IsInitialized() const;

public:
	///	<summary>
	///		Generate the unstructured grid information for a
	///		longitude-latitude grid.
	///	</summary>
	void GenerateLatitudeLongitude(
		const DataArray1D<double> & vecLat,
		const DataArray1D<double> & vecLon,
		bool fRegional
	);

	///	<summary>
	///		Try to automatically generate the SimpleGrid from a NetCDF
	///		file with latitude/longitude coordinates.
	///	</summary>
	void GenerateLatitudeLongitude(
		NcFile * ncFile,
		bool fRegional,
		const std::string & strLatitudeName,
		const std::string & strLongitudeName
	);

	///	<summary>
	///		Try to automatically generate the SimpleGrid from a NetCDF
	///		file with latitude/longitude coordinates.
	///	</summary>
	void GenerateLatitudeLongitude(
		NcFile * ncFile,
		bool fRegional
	);

	///	<summary>
	///		Generate the unstructured grid information for a rectilinear
	///		stereographic grid at the given point.
	///	</summary>
	void GenerateRectilinearStereographic(
		double dLonRad0,
		double dLatRad0,
		int nX,
		double dDeltaX,
		bool fCalculateArea = false
	);

	///	<summary>
	///		Generate the unstructured grid information for a radial
	///		stereographic grid at the given point.
	///	</summary>
	void GenerateRadialStereographic(
		double dLonRad0,
		double dLatRad0,
		int nR,
		int nA,
		double dDeltaR,
		bool fCalculateArea = false
	);

	///	<summary>
	///		Initialize the SimpleGrid from a Mesh assuming the mesh is a finite
	///		volume mesh.
	///	</summary>
	void FromMeshFV(
		const Mesh & mesh
	);

	///	<summary>
	///		Initialize the SimpleGrid from a Mesh assuming the mesh is a finite
	///		element mesh.
	///	</summary>
	void FromMeshFE(
		const Mesh & mesh,
		bool fCGLL,
		int nP
	);

	///	<summary>
	///		Read the grid information from a file.
	///	</summary>
	void FromFile(
		const std::string & strConnectivityFile
	);

	///	<summary>
	///		Write the grid information to a file.
	///	</summary>
	void ToFile(
		const std::string & strConnectivityFile
	) const;

	///	<summary>
	///		Get the size of the SimpleGrid (number of points).
	///	</summary>
	size_t GetSize() const {
		return (m_vecConnectivity.size());
	}

	///	<summary>
	///		Convert a coordinate to an index.
	///	</summary>
	int CoordinateVectorToIndex(
		const std::vector<int> & coordvec
	) const;

public:
	///	<summary>
	///		Build a kdtree using this SimpleGrid.
	///	</summary>
	void BuildKDTree();

	///	<summary>
	///		Find the nearest node to the given coordinate.
	///	</summary>
	size_t NearestNode(
		double dLonRad,
		double dLatRad
	) const;

public:
	///	<summary>
	///		Grid dimensions.
	///	</summary>
	std::vector<size_t> m_nGridDim;

	///	<summary>
	///		Longitude of each grid point (in radians).
	///	</summary>
	DataArray1D<double> m_dLon;

	///	<summary>
	///		Latitude of each grid point (in radians).
	///	</summary>
	DataArray1D<double> m_dLat;

	///	<summary>
	///		Area of the grid cell (in m^2).
	///	</summary>
	DataArray1D<double> m_dArea;

	///	<summary>
	///		Connectivity of each grid point.
	///	</summary>
	std::vector< std::vector<int> > m_vecConnectivity;

private:
	///	<summary>
	///		kd tree used for quick lookup of grid points.
	///	</summary>
	kdtree * m_kdtree;
};

///////////////////////////////////////////////////////////////////////////////

#endif

