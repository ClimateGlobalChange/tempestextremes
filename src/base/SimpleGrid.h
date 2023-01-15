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

	///	<summary>
	///		Determine if the SimpleGrid has calculated areas.
	///	</summary>
	bool HasAreas() const {
		if (m_dArea.GetRows() == 0) {
			return false;
		}
		_ASSERT(m_dArea.GetRows() == m_dLon.GetRows());
		return true;
	}

	///	<summary>
	///		Determine if the SimpleGrid has connectivity information.
	///	</summary>
	bool HasConnectivity() const {
		if (m_vecConnectivity.size() == 0) {
			return false;
		}
		_ASSERT(m_vecConnectivity.size() == m_dLon.GetRows());
		return true;
	}

public:
	///	<summary>
	///		Generate connectivity information for a rectilinear grid.
	///	</summary>
	void GenerateRectilinearConnectivity(
		int nLat,
		int nLon,
		bool fRegional,
		bool fDiagonalConnectivity
	);

	///	<summary>
	///		Get the latitude name, variable and dimension from a NcFile.
	///	</summary>
	void GetLatitudeFromNcFile(
		NcFile * ncFile,
		std::string & strLatitudeName,
		NcVar ** pvarLat,
		NcDim ** pdimLat
	);

	///	<summary>
	///		Get the latitude name, variable and dimension from a NcFile.
	///	</summary>
	void GetLongitudeFromNcFile(
		NcFile * ncFile,
		std::string & strLongitudeName,
		NcVar ** pvarLon,
		NcDim ** pdimLon
	);

	///	<summary>
	///		Generate the unstructured grid information for a
	///		longitude-latitude grid.
	///	</summary>
	void GenerateLatitudeLongitude(
		const DataArray1D<double> & vecLat,
		const DataArray1D<double> & vecLon,
		bool fRegional,
		bool fDiagonalConnectivity,
		bool fVerbose
	);

	///	<summary>
	///		Try to automatically generate the SimpleGrid from a NetCDF
	///		file with latitude/longitude coordinates.
	///	</summary>
	void GenerateLatitudeLongitude(
		NcFile * ncFile,
		std::string & strLatitudeName,
		std::string & strLongitudeName,
		bool fRegional,
		bool fDiagonalConnectivity
	);

	///	<summary>
	///		Try to automatically generate the SimpleGrid from a NetCDF
	///		file with latitude/longitude coordinates.
	///	</summary>
	void GenerateLatitudeLongitude(
		NcFile * ncFile,
		bool fRegional,
		bool fDiagonalConnectivity
	);

	///	<summary>
	///		Generate the unstructured grid information for a
	///		longitude-latitude grid.
	///	</summary>
	void GenerateRegionalLatitudeLongitude(
		double dLatRad1,
		double dLatRad2,
		double dLonRad1,
		double dLonRad2,
		int nLat,
		int nLon,
		bool fDiagonalConnectivity
	);

	///	<summary>
	///		Generate the unstructured grid information for a rectilinear
	///		stereographic grid at the given point.
	///	</summary>
	void GenerateRectilinearStereographic(
		double dLonRad0,
		double dLatRad0,
		int nX,
		double dDeltaXRad,
		bool fFlipSouthernHemisphere,
		bool fCalculateArea
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
		double dDeltaRRad,
		bool fCalculateArea
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
	///		Read the coordinate information from a data file (note that connectivity
	///		information won't be available in this case).
	///	</summary>
	void FromUnstructuredDataFile(
		const std::string & strDataFile,
		const std::string & strLatitudeName,
		const std::string & strLongitudeName
	);

	///	<summary>
	///		Get the number of dimensions of the SimpleGrid.
	///	</summary>
	size_t DimCount() const {
		return (m_nGridDim.size());
	}

	///	<summary>
	///		Get the size of the SimpleGrid (number of points).
	///	</summary>
	size_t GetSize() const {
		_ASSERT(m_dLat.GetRows() == m_dLon.GetRows());
		return (m_dLon.GetRows());
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

	///	<summary>
	///		Find the set of nodes within the specified distance (degrees great circle distance)
	///		of the given coordinate.
	///	</summary>
	void NearestNodes(
		double dLonRad,
		double dLatRad,
		double dDistDegGCD,
		std::vector<size_t> & vecNodeIxs
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

public:
	///	<summary>
	///		Area of the grid cell (in sr) (optionally initialized).
	///	</summary>
	DataArray1D<double> m_dArea;

	///	<summary>
	///		Connectivity of each grid point (optionally initialized).
	///	</summary>
	std::vector< std::vector<int> > m_vecConnectivity;

private:
	///	<summary>
	///		kd tree used for quick lookup of grid points (optionally initialized).
	///	</summary>
	kdtree * m_kdtree;
};

///////////////////////////////////////////////////////////////////////////////

#endif

