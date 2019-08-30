///////////////////////////////////////////////////////////////////////////////
///
///	\file    SimpleGrid.cpp
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

#include "SimpleGrid.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

const char * SimpleGrid::c_szFileIdentifier =
	"#TempestGridConnectivityFileV2.0";

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	const DataArray1D<double> & vecLat,
	const DataArray1D<double> & vecLon,
	bool fRegional
) {
	int nLat = vecLat.GetRows();
	int nLon = vecLon.GetRows();

	m_dLat.Allocate(nLon * nLat);
	m_dLon.Allocate(nLon * nLat);
	m_vecConnectivity.resize(nLon * nLat);

	m_nGridDim.resize(2);
	m_nGridDim[0] = nLat;
	m_nGridDim[1] = nLon;

	// Verify units of latitude and longitude
	for (int j = 0; j < nLat; j++) {
		if (fabs(vecLat[j]) > 0.5 * M_PI + 1.0e-12) {
			_EXCEPTIONT("In SimpleGrid, latitude array must be given in radians");
		}
	}

	int ixs = 0;
	for (int j = 0; j < nLat; j++) {
	for (int i = 0; i < nLon; i++) {

		// Vectorize coordinates
		m_dLat[ixs] = vecLat[j];
		m_dLon[ixs] = vecLon[i];

		// Connectivity in each compass direction
		if (j != 0) {
			m_vecConnectivity[ixs].push_back((j-1) * nLon + i);
		}
		if (j != nLat-1) {
			m_vecConnectivity[ixs].push_back((j+1) * nLon + i);
		}

		if ((!fRegional) ||
		    ((i != 0) && (i != nLon-1))
		) {
			m_vecConnectivity[ixs].push_back(
				j * nLon + ((i + 1) % nLon));
			m_vecConnectivity[ixs].push_back(
				j * nLon + ((i + nLon - 1) % nLon));
		}

		ixs++;
	}
	}

}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	NcFile * ncFile,
	bool fRegional,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName
) {
	NcDim * dimLat = ncFile->get_dim(strLatitudeName.c_str());
	if (dimLat == NULL) {
		_EXCEPTION1("No dimension \"%s\" found in input file",
			strLatitudeName.c_str());
	}

	NcDim * dimLon = ncFile->get_dim(strLongitudeName.c_str());
	if (dimLon == NULL) {
		_EXCEPTION1("No dimension \"%s\" found in input file",
			strLongitudeName.c_str());
	}

	NcVar * varLat = ncFile->get_var(strLatitudeName.c_str());
	if (varLat == NULL) {
		_EXCEPTION1("No variable \"%s\" found in input file",
			strLatitudeName.c_str());
	}

	NcVar * varLon = ncFile->get_var(strLongitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION1("No variable \"%s\" found in input file",
			strLongitudeName.c_str());
	}

	int nLat = dimLat->size();
	int nLon = dimLon->size();

	DataArray1D<double> vecLat(nLat);
	varLat->get(vecLat, nLat);

	for (int j = 0; j < nLat; j++) {
		vecLat[j] *= M_PI / 180.0;
	}

	DataArray1D<double> vecLon(nLon);
	varLon->get(vecLon, nLon);

	for (int i = 0; i < nLon; i++) {
		vecLon[i] *= M_PI / 180.0;
	}

	// Generate the SimpleGrid
	GenerateLatitudeLongitude(vecLat, vecLon, fRegional);
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	NcFile * ncFile,
	bool fRegional
) {
	return GenerateLatitudeLongitude(ncFile, fRegional, "lat", "lon");
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::FromFile(
	const std::string & strConnectivityFile
) {
	if ((m_nGridDim.size() != 0) ||
	    (m_dLon.GetRows() != 0) ||
	    (m_dLat.GetRows() != 0) ||
	    (m_dArea.GetRows() != 0) ||
	    (m_vecConnectivity.size() != 0)
	) {
		_EXCEPTIONT("Attempting to read previously initialized SimpleGrid");
	}

	std::ifstream fsGrid(strConnectivityFile.c_str());
	if (!fsGrid.is_open()) {
		_EXCEPTION1("Unable to open file \"%s\"",
			strConnectivityFile.c_str());
	}

	std::string strConnectivityString;
	fsGrid >> strConnectivityString;

	if (strConnectivityString != c_szFileIdentifier) {
		_EXCEPTION1("Invalid connectivity file format \"%s\"",
			strConnectivityFile.c_str());
	}

	size_t sDims;
	fsGrid >> sDims;

	size_t sFaces = 1;

	if ((sDims < 1) || (sDims > 2)) {
		_EXCEPTION1("Invalid connectivity file: %lu dimensions out "
			"of range (expected 1,2)", sDims);
	}

	m_nGridDim.resize(sDims);
	for (size_t s = 0; s < sDims; s++) {
		char cComma;
		fsGrid >> cComma;
		fsGrid >> m_nGridDim[s];

		if (m_nGridDim[s] < 1) {
			_EXCEPTION2("Grid dimension %lu out of range (%lu found)",
				s, m_nGridDim[s]);
		}

		sFaces *= m_nGridDim[s];
	}

	m_dLon.Allocate(sFaces);
	m_dLat.Allocate(sFaces);
	m_dArea.Allocate(sFaces);
	m_vecConnectivity.resize(sFaces);

	for (size_t f = 0; f < sFaces; f++) {
		size_t sNeighbors;
		char cComma;
		fsGrid >> m_dLon[f];
		fsGrid >> cComma;
		fsGrid >> m_dLat[f];
		fsGrid >> cComma;
		fsGrid >> m_dArea[f];
		fsGrid >> cComma;
		fsGrid >> sNeighbors;
		fsGrid >> cComma;

		// Convert to radians
		m_dLon[f] *= M_PI / 180.0;
		m_dLat[f] *= M_PI / 180.0;

		// Load connectivity
		m_vecConnectivity[f].resize(sNeighbors);
		for (size_t n = 0; n < sNeighbors; n++) {
			fsGrid >> m_vecConnectivity[f][n];
			if (n != sNeighbors-1) {
				fsGrid >> cComma;
			}
			m_vecConnectivity[f][n]--;
		}
		if (fsGrid.eof()) {
			if (f != sFaces-1) {
				_EXCEPTIONT("Premature end of file");
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::ToFile(
	const std::string & strConnectivityFile
) {
	std::ofstream fsOutput(strConnectivityFile.c_str());
	if (!fsOutput.is_open()) {
		_EXCEPTION1("Cannot open output file \"%s\"",
			strConnectivityFile.c_str());
	}

	fsOutput << std::setprecision(14);
	fsOutput << std::scientific;

	fsOutput << c_szFileIdentifier << std::endl;
	
	size_t sFaces = 1;
	fsOutput << m_nGridDim.size();
	for (size_t i = 0; i < m_nGridDim.size(); i++) {
		sFaces *= m_nGridDim[i];
		fsOutput << "," << m_nGridDim[i];
	}

	if (m_dLon.GetRows() != sFaces) {
		_EXCEPTIONT("Mangled SimpleGrid structure: m_dLon.size() != size from m_nGridDim");
	}
	if (sFaces != m_dLat.GetRows()) {
		_EXCEPTIONT("Mangled SimpleGrid structure: m_dLon.size() != m_dLat.size()");
	}
	if (sFaces != m_dArea.GetRows()) {
		_EXCEPTIONT("Mangled SimpleGrid structure: m_dLon.size() != m_dArea.size()");
	}
	if (sFaces != m_vecConnectivity.size()) {
		_EXCEPTIONT("Mangled SimpleGrid structure: m_dLon.size() != m_vecConnectivity.size()");
	}

	fsOutput << sFaces << std::endl;
	for (size_t i = 0; i < sFaces; i++) {
		fsOutput << m_dLon[i] << "," << m_dLat[i] << ","
			<< m_dArea[i] << "," << m_vecConnectivity[i].size() << std::endl;
		for (size_t j = 0; j < m_vecConnectivity[i].size(); i++) {
			fsOutput << "," << m_vecConnectivity[i][j];
		}
		fsOutput << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

int SimpleGrid::CoordinateVectorToIndex(
	const std::vector<int> & coordvec
) {
	if (m_nGridDim.size() == 0) {
		_EXCEPTIONT("Invalid SimpleGrid");
	}

	if (coordvec.size() != m_nGridDim.size()) {
		_EXCEPTIONT("Invalid coordinate vector");
	}
	if (coordvec.size() == 1) {
		if (coordvec[0] >= m_nGridDim[0]) {
			_EXCEPTIONT("Coordinate vector out of range");
		}
		return coordvec[0];
	}
	if (coordvec.size() == 2) {
		if (coordvec[0] >= m_nGridDim[0]) {
			_EXCEPTIONT("Coordinate vector out of range");
		}
		if (coordvec[1] >= m_nGridDim[1]) {
			_EXCEPTIONT("Coordinate vector out of range");
		}
	}

	int ix = 0;
	int d = 1;
	for (int i = 0; i < coordvec.size(); i++) {
		if (coordvec[i] >= m_nGridDim[i]) {
			_EXCEPTIONT("Coordinate vector out of range");
		}
		ix = ix + i * d;
		d = d * m_nGridDim[i];
	}

	return ix;
}

///////////////////////////////////////////////////////////////////////////////

