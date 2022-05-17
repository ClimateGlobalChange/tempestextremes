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
#include "GridElements.h"
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "Announce.h"
#include "Constants.h"
#include "CoordTransforms.h"
#include "NetCDFUtilities.h"
#include "kdtree.h"

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

SimpleGrid::~SimpleGrid() {
	if (m_kdtree != NULL) {
		kd_free(m_kdtree);
	}
}

///////////////////////////////////////////////////////////////////////////////

bool SimpleGrid::IsInitialized() const {
	if ((m_nGridDim.size() != 0) ||
	    m_dLon.IsAttached() ||
	    m_dLat.IsAttached() ||
	    m_dArea.IsAttached() ||
	    (m_vecConnectivity.size() != 0) ||
		(m_kdtree != NULL)
	) {
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateRectilinearConnectivity(
	int nLat,
	int nLon,
	bool fRegional,
	bool fDiagonalConnectivity
) {
	m_vecConnectivity.clear();
	m_vecConnectivity.resize(nLon * nLat);

	size_t ixs = 0;
	for (int j = 0; j < nLat; j++) {
	for (int i = 0; i < nLon; i++) {

		// Connectivity in eight directions
		if (fDiagonalConnectivity) {
			if (fRegional) {
				for (int ix = -1; ix <= 1; ix++) {
				for (int jx = -1; jx <= 1; jx++) {
					if ((ix == 0) && (jx == 0)) {
						continue;
					}

					int inew = i + ix;
					int jnew = j + jx;

					if ((inew < 0) || (inew >= nLon)) {
						continue;
					}
					if ((jnew < 0) || (jnew >= nLat)) {
						continue;
					}

					m_vecConnectivity[ixs].push_back(jnew * nLon + inew);
				}
				}

			} else {
				for (int ix = -1; ix <= 1; ix++) {
				for (int jx = -1; jx <= 1; jx++) {
					if ((ix == 0) && (jx == 0)) {
						continue;
					}

					int inew = i + ix;
					int jnew = j + jx;

					if ((jnew < 0) || (jnew >= nLat)) {
						continue;
					}
					if (inew < 0) {
						inew += nLon;
					}
					if (inew >= nLon) {
						inew -= nLon;
					}

					m_vecConnectivity[ixs].push_back(jnew * nLon + inew);
				}
				}
			}

		// Connectivity in the four primary directions
		} else {
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
		}

		ixs++;
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	const DataArray1D<double> & vecLat,
	const DataArray1D<double> & vecLon,
	bool fRegional,
	bool fDiagonalConnectivity,
	bool fVerbose
) {
	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call GenerateLatitudeLongitude() on previously initialized grid");
	}

	int nLat = vecLat.GetRows();
	int nLon = vecLon.GetRows();

	if (nLat < 2) {
		_EXCEPTIONT("At least two latitudes needed to generate grid.");
	}
	if (nLon < 2) {
		_EXCEPTIONT("At least two longitudes needed to generate grid.");
	}

	m_dLat.Allocate(nLon * nLat);
	m_dLon.Allocate(nLon * nLat);
	m_dArea.Allocate(nLon * nLat);

	m_nGridDim.resize(2);
	m_nGridDim[0] = nLat;
	m_nGridDim[1] = nLon;

	// Verify units of latitude and longitude
	bool fCalculateArea = true;

	for (int j = 0; j < nLat; j++) {
		if (fabs(vecLat[j]) > 0.5 * M_PI + 1.0e-12) {
			if (fRegional) {
				Announce("WARNING: Latitude array out of bounds [-90,90] "
					"(%1.5f); defaulting grid areas to 1.",
					RadToDeg(vecLat[j]));
				fCalculateArea = false;
				break;
			} else {
				_EXCEPTION1("Latitude array out of bounds [-90,90] (%1.5f). "
					 "Did you mean to specify \"--regional\"?",
					 RadToDeg(vecLat[j]));
			}
		}
	}
	for (int i = 0; i < nLon; i++) {
		if (fabs(vecLon[i]) > 3.0 * M_PI + 1.0e-12) {
			if (fRegional) {
				Announce("WARNING: Longitude array out of bounds [-540,540] "
					"(%1.5f); defaulting grid areas to 1.",
					RadToDeg(vecLon[i]));
				fCalculateArea = false;
				break;
			} else {
				_EXCEPTION1("Longitude array out of bounds [-540,540] (%1.5f). "
					"Did you mean to specify \"--regional\"?",
					RadToDeg(vecLon[i]));
			}
		}
	}

	// Determine orientation of latitude array
	double dLatOrient = 1.0;
	if (vecLat[1] < vecLat[0]) {
		dLatOrient = -1.0;
	}
	for (int j = 0; j < nLat-1; j++) {
		if (dLatOrient * vecLat[1] < dLatOrient * vecLat[0]) {
			_EXCEPTIONT("Latitude array must be monotone.");
		}
	}

	int ixs = 0;
	for (int j = 0; j < nLat; j++) {
	for (int i = 0; i < nLon; i++) {

		// Vectorize coordinates
		m_dLat[ixs] = vecLat[j];
		m_dLon[ixs] = vecLon[i];

		// Bounds of each volume and associated area
		double dLatRad1;
		double dLatRad2;
		
		double dLonRad1;
		double dLonRad2;

		if (j == 0) {
			if (fRegional) {
				dLatRad1 = vecLat[0] - 0.5 * (vecLat[1] - vecLat[0]);
			} else {
				dLatRad1 = - dLatOrient * 0.5 * M_PI;
			}
		} else {
			dLatRad1 = 0.5 * (vecLat[j-1] + vecLat[j]);
		}
		if (j == nLat-1) {
			if (fRegional) {
				dLatRad2 = vecLat[j] + 0.5 * (vecLat[j] - vecLat[j-1]);
			} else {
				dLatRad2 = dLatOrient * 0.5 * M_PI;
			}
		} else {
			dLatRad2 = 0.5 * (vecLat[j+1] + vecLat[j]);
		}

		if (i == 0) {
			if (fRegional) {
				dLonRad1 = vecLon[0] - 0.5 * (vecLon[1] - vecLon[0]);
			} else {
				dLonRad1 = AverageLongitude_Rad(vecLon[0], vecLon[nLon-1]);
			}
		} else {
			dLonRad1 = AverageLongitude_Rad(vecLon[i-1], vecLon[i]);
		}

		if (i == nLon-1) {
			if (fRegional) {
				dLonRad2 = vecLon[i] + 0.5 * (vecLon[i] - vecLon[i-1]);
			} else {
				dLonRad2 = AverageLongitude_Rad(vecLon[nLon-1], vecLon[0]);
			}
		} else {
			dLonRad2 = AverageLongitude_Rad(vecLon[i], vecLon[i+1]);
		}

		// Because of the way AverageLongitude_Rad works,
		if (dLonRad1 > dLonRad2) {
			dLonRad1 -= 2.0 * M_PI;
		}
		double dDeltaLon = dLonRad2 - dLonRad1;
		if (!fRegional && (dDeltaLon >= M_PI)) {
			_EXCEPTION1("Grid element longitudinal extent too large (%1.7f deg).  Did you mean to specify \"--regional\"?",
				dDeltaLon * 180.0 / M_PI);
		}

		if (fCalculateArea) {
			m_dArea[ixs] = fabs(sin(dLatRad2) - sin(dLatRad1)) * dDeltaLon;
		} else {
			m_dArea[ixs] = 1.0;
		}

		ixs++;
	}
	}

	// Generate connectivity
	GenerateRectilinearConnectivity(nLat, nLon, fRegional, fDiagonalConnectivity);

	// Output total area
	{
		double dTotalArea = 0.0;
		for (size_t i = 0; i < m_dArea.GetRows(); i++) {
			dTotalArea += m_dArea[i];
		}
		if (fVerbose) {
			Announce("Total calculated grid area: %1.15e sr", dTotalArea);
		}
	}

}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	NcFile * ncFile,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName,
	bool fRegional,
	bool fDiagonalConnectivity
) {
	_ASSERT(ncFile != NULL);

	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call GenerateLatitudeLongitude() on previously initialized grid");
	}

	// Load latitude variable
	NcVar * varLat;
	NcDim * dimLat;

	std::vector<std::string> vecLatitudeNames;
	if (strLatitudeName == std::string("[auto]")) {
		vecLatitudeNames.push_back("lat");
		vecLatitudeNames.push_back("latitude");
		vecLatitudeNames.push_back("LAT");
		vecLatitudeNames.push_back("latitude0");
		vecLatitudeNames.push_back("Latitude");
		vecLatitudeNames.push_back("XLAT");
	} else {
		vecLatitudeNames.push_back(strLatitudeName);
	}

	size_t sLatIx = NcGetVarFromList(*ncFile, vecLatitudeNames, &varLat, &dimLat);
	if (sLatIx == vecLatitudeNames.size()) {
		_EXCEPTION1("No variable %s found in input file",
			STLStringHelper::ConcatenateStringVector(vecLatitudeNames, ", ").c_str());
	}
	if ((varLat->num_dims() < 1) || (varLat->num_dims() > 2)) {
		_EXCEPTION1("In NetCDF file latitude variable \"%s\" must have either one or two dimensions",
			varLat->name());

	} else if (varLat->num_dims() == 1) {
		if (dimLat == NULL) {
			_EXCEPTION1("In NetCDF file latitude variable \"%s\" must have dimension with same name",
				varLat->name());
		}
	}

	// Load longitude variable
	NcVar * varLon;
	NcDim * dimLon;

	std::vector<std::string> vecLongitudeNames;
	if (strLatitudeName == std::string("[auto]")) {
		vecLongitudeNames.push_back("lon");
		vecLongitudeNames.push_back("longitude");
		vecLongitudeNames.push_back("LON");
		vecLongitudeNames.push_back("longitude0");
		vecLongitudeNames.push_back("Longitude");
		vecLongitudeNames.push_back("XLONG");
	} else {
		vecLongitudeNames.push_back(strLongitudeName);
	}

	size_t sLonIx = NcGetVarFromList(*ncFile, vecLongitudeNames, &varLon, &dimLon);
	if (sLonIx == vecLatitudeNames.size()) {
		_EXCEPTION1("No variable %s found in input file",
			STLStringHelper::ConcatenateStringVector(vecLatitudeNames, ", ").c_str());
	}
	if ((varLon->num_dims() < 1) || (varLon->num_dims() > 2)) {
		_EXCEPTION1("In NetCDF file longitude variable \"%s\" must have either one or two dimensions",
			varLon->name());

	} else if (varLon->num_dims() == 1) {
		if (dimLon == NULL) {
			_EXCEPTION1("In NetCDF file longitude variable \"%s\" must have dimension with same name",
				varLon->name());
		}
	}

	if (varLat->num_dims() != varLon->num_dims()) {
		_EXCEPTION2("Latitude variable \"%s\" and longitude variable \"%s\" must have same number of dimensions",
			varLat->name(), varLon->name());
	}

	// RLL mesh (lat and lon are vectors)
	if (varLat->num_dims() == 1) {

		// Load lat/lon data
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
		GenerateLatitudeLongitude(
			vecLat,
			vecLon,
			fRegional,
			fDiagonalConnectivity,
			true);                  // Verbosity enabled
	}

	// Other rectilinear projection (lat and lon are arrays)
	if (varLat->num_dims() == 2) {

		long lY = varLat->get_dim(0)->size();
		long lX = varLat->get_dim(1)->size();

		m_nGridDim.resize(2);
		m_nGridDim[0] = lY;
		m_nGridDim[1] = lX;

		if (lY != varLon->get_dim(0)->size()) {
			_EXCEPTION4("Latitude variable \"%s\" and longitude variable \"%s\" have incompatible first dimension size (%li != %li)",
				varLat->name(), varLat->get_dim(0)->size(),
				varLon->name(), varLon->get_dim(0)->size());
		}
		if (lX != varLon->get_dim(1)->size()) {
			_EXCEPTION4("Latitude variable \"%s\" and longitude variable \"%s\" have incompatible second dimension size (%li != %li)",
				varLat->name(), varLat->get_dim(1)->size(),
				varLon->name(), varLon->get_dim(1)->size());
		}

		m_dLon.Allocate(lX * lY);
		m_dLat.Allocate(lX * lY);

		varLat->get(m_dLat, lY, lX);
		varLon->get(m_dLon, lY, lX);

		for (size_t s = 0; s < m_dLon.GetRows(); s++) {
			m_dLon[s] = DegToRad(m_dLon[s]);
			m_dLat[s] = DegToRad(m_dLat[s]);
		}

		// Approximate cell areas
		if ((lX > 1) && (lY > 1)) {
			m_dArea.Allocate(lX * lY);

			double dTotalArea = 0.0;
			size_t s = 0;
			for (size_t j = 0; j < static_cast<size_t>(lY); j++) {
			for (size_t i = 0; i < static_cast<size_t>(lX); i++) {
				double dLonRad0;
				double dLonRad1;
				double dLatRad0;
				double dLatRad1;
				if (i == 0) {
					dLonRad0 = 1.5 * m_dLon[s] - 0.5 * m_dLon[s+1];
				} else {
					dLonRad0 = AverageLongitude_Rad(m_dLon[s-1], m_dLon[s]);
				}
				if (i == static_cast<size_t>(lX)-1) {
					dLonRad1 = 1.5 * m_dLon[s] - 0.5 * m_dLon[s-1];
				} else {
					dLonRad1 = AverageLongitude_Rad(m_dLon[s], m_dLon[s+1]);
				}
				if (j == 0) {
					dLatRad0 = 1.5 * m_dLat[s] - 0.5 * m_dLat[s+static_cast<size_t>(lX)];
				} else {
					dLatRad0 = 0.5 * m_dLat[s-static_cast<size_t>(lX)] + 0.5 * m_dLat[s];
				}
				if (j == static_cast<size_t>(lY)-1) {
					dLatRad1 = 1.5 * m_dLat[s] - 0.5 * m_dLat[s-static_cast<size_t>(lX)];
				} else {
					dLatRad1 = 0.5 * m_dLat[s] + 0.5 * m_dLat[s+static_cast<size_t>(lX)];
				}

				m_dArea[s] =
					fabs(sin(dLatRad1) - sin(dLatRad0))
					* fabs(dLonRad1 - dLonRad0);

				dTotalArea += m_dArea[s];

				s++;
			}
			}

			Announce("Approximate domain area: %1.15e sr", dTotalArea);
 		}

		GenerateRectilinearConnectivity(lY, lX, fRegional, fDiagonalConnectivity);
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateLatitudeLongitude(
	NcFile * ncFile,
	bool fRegional,
	bool fDiagonalConnectivity
) {
	return
		GenerateLatitudeLongitude(
			ncFile,
			"lat",
			"lon",
			fRegional,
			fDiagonalConnectivity);
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateRegionalLatitudeLongitude(
	double dLatRad1,
	double dLatRad2,
	double dLonRad1,
	double dLonRad2,
	int nLat,
	int nLon,
	bool fDiagonalConnectivity
) {
	_ASSERT(nLat >= 1);
	_ASSERT(nLon >= 1);
	if ((dLatRad1 < -0.5 * M_PI) || (dLatRad1 > 0.5 * M_PI)) {
		_EXCEPTION1("Latitude-longitude grid latitude bound out of range (%1.5f)", dLatRad1);
	}
	if ((dLatRad2 < -0.5 * M_PI) || (dLatRad2 > 0.5 * M_PI)) {
		_EXCEPTION1("Latitude-longitude grid latitude bound out of range (%1.5f)", dLatRad2);
	}
	if (fabs(dLatRad1 - dLatRad2) < 1.0e-12) {
		_EXCEPTION1("Latitude-longitude grid latitude bounds must not be equal (%1.5f)", dLatRad1);
	}
	if (fabs(dLonRad1 - dLonRad2) < 1.0e-12) {
		_EXCEPTION1("Latitude-longitude grid longitude bounds must not be equal (%1.5f)", dLonRad1);
	}

	DataArray1D<double> vecLat(nLat);
	DataArray1D<double> vecLon(nLon);

	double dDeltaLatRad = (dLatRad2 - dLatRad1) / static_cast<double>(nLat);
	double dDeltaLonRad = (dLonRad2 - dLonRad1) / static_cast<double>(nLon);

	for (int j = 0; j < nLat; j++) {
		vecLat[j] = dLatRad1 + dDeltaLatRad * (static_cast<double>(j) + 0.5);
	}
	for (int i = 0; i < nLon; i++) {
		vecLon[i] = dLonRad1 + dDeltaLonRad * (static_cast<double>(i) + 0.5);
	}

	GenerateLatitudeLongitude(
		vecLat,
		vecLon,
		true,                   // Regional
		fDiagonalConnectivity,
		false);                 // Verbosity disabled
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateRectilinearStereographic(
	double dLonRad0,
	double dLatRad0,
	int nX,
	double dDeltaXRad,
	bool fFlipSouthernHemisphere,
	bool fCalculateArea
) {
	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call GenerateRectilinearStereographic() on previously initialized grid");
	}

	_ASSERT(nX >= 1);
	_ASSERT(dDeltaXRad > 0.0);
	_ASSERT(fabs(dLatRad0 <= 0.5 * M_PI));

	if (fabs(dLatRad0 - 0.5 * M_PI) < ReferenceTolerance) {
		dLatRad0 = 0.5 * M_PI;
	}
	if (fabs(dLatRad0 + 0.5 * M_PI) < ReferenceTolerance) {
		dLatRad0 = - 0.5 * M_PI;
	}

	m_nGridDim.resize(2);
	m_nGridDim[0] = nX;
	m_nGridDim[1] = nX;

	m_dLon.Allocate(nX * nX);
	m_dLat.Allocate(nX * nX);

	const double dXgcd0 = - 0.5 * dDeltaXRad * static_cast<double>(nX - 1);

	if (dXgcd0 < - M_PI + ReferenceTolerance) {
		_EXCEPTION1("Total angular coverage of rectilinear stereographic "
			"grid too large (%1.5f <= -pi/2)", dXgcd0);
	}

	// Get coordinates in the space of the stereographic projection
	DataArray1D<double> dXs(nX);
	DataArray1D<double> dYs(nX);

	for (int j = 0; j < nX; j++) {
		double dYgcd = dXgcd0 + dDeltaXRad * static_cast<double>(j);
		dYs[j] = sqrt(4.0 * (1.0 - cos(dYgcd)) / (1.0 + cos(dYgcd)));
		if (dYgcd < 0.0) {
			dYs[j] *= -1.0;
		}
	}

	for (int i = 0; i < nX; i++) {
		double dXgcd = dXgcd0 + dDeltaXRad * static_cast<double>(i);
		dXs[i] = sqrt(4.0 * (1.0 - cos(dXgcd)) / (1.0 + cos(dXgcd)));
		if (dXgcd < 0.0) {
			dXs[i] *= -1.0;
		}
	}

	if (fFlipSouthernHemisphere && (dLatRad0 < 0.0)) {
		for (int j = 0; j < nX; j++) {
			dYs[j] *= -1.0;
		}
	}

	// Store longitude and latitude of centerpoints
	int s = 0;
	for (int j = 0; j < nX; j++) {
	for (int i = 0; i < nX; i++) {
		StereographicProjectionInv(
			dLonRad0,
			dLatRad0,
			dXs[i],
			dYs[j],
			m_dLon[s],
			m_dLat[s]);
/*
		// DEBUG: Verify that the great circle distance is correct
		printf("%i %i %1.5f\n", i, j,
			180.0 / M_PI
			* GreatCircleDistance_Rad(
				dLonRad0,
				dLatRad0,
				m_dLon[s],
				m_dLat[s]));
*/
		s++;
	}
	}

	// Calculate the area
	if (fCalculateArea) {
		_EXCEPTIONT("Unable to calculate the area of the RectilinearStereographic grid (not impemented)");
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::GenerateRadialStereographic(
	double dLonRad0,
	double dLatRad0,
	int nR,
	int nA,
	double dDeltaRRad,
	bool fCalculateArea
) {
	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call GenerateRadialStereographic() on previously initialized grid");
	}

	if (nA < 8) {
		_EXCEPTIONT("Minimum of 8 azimuthal slices allowed");
	}

	_ASSERT(nR >= 1);
	_ASSERT(dDeltaRRad > 0.0);
	_ASSERT(fabs(dLatRad0 <= 0.5 * M_PI));

	const double dRgcdmax = (static_cast<double>(nR-1) + 0.5) * dDeltaRRad;

	if (dRgcdmax >= M_PI) {
		_EXCEPTION1("Total angular coverage of radial stereographic "
			"grid too large (%1.5f >= pi)", dRgcdmax);
	}

	m_nGridDim.resize(2);
	m_nGridDim[0] = nA;
	m_nGridDim[1] = nR;

	m_dLon.Allocate(nA * nR);
	m_dLat.Allocate(nA * nR);

	// Get the reference coordinates in the azimuthal and radial directions
	DataArray1D<double> dXs(nA);
	DataArray1D<double> dYs(nA);

	for (int i = 0; i < nA; i++) {
		double dAzimuth = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(nA);
		dXs[i] = cos(dAzimuth);
		dYs[i] = sin(dAzimuth);
	}

	DataArray1D<double> dRs(nR);
	for (int i = 0; i < nR; i++) {
		double dRgcd = (static_cast<double>(i) + 0.5) * dDeltaRRad;
		dRs[i] = 2.0 * sqrt((1.0 - cos(dRgcd)) / (1.0 + cos(dRgcd)));
	}

	// Calculate the lon/lat coordinates
	int s = 0;
	for (int j = 0; j < nR; j++) {
	for (int i = 0; i < nA; i++) {
		StereographicProjectionInv(
			dLonRad0,
			dLatRad0,
			dXs[i] * dRs[j],
			dYs[i] * dRs[j],
			m_dLon[s],
			m_dLat[s]);
/*
		// DEBUG: Verify that the great circle distance is uniform
		// in the radial direction
		printf("%i %i %1.5f %1.5f %1.5f\n", i, j,
			dXs[i] * dRs[j],
			dYs[i] * dRs[j],
			180.0 / M_PI
			* GreatCircleDistance_Rad(
				dLonRad0,
				dLatRad0,
				m_dLon[s],
				m_dLat[s]));
*/
		s++;
	}
	}

	// Calculate the area
	if (fCalculateArea) {
		_EXCEPTIONT("Unable to calculate the area of the RectilinearStereographic grid (not impemented)");
	}

}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::FromMeshFV(
	const Mesh & mesh
) {
	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call FromMeshFV() on previously initialized grid");
	}

	if (mesh.vecFaceArea.GetRows() == 0) {
		_EXCEPTIONT("Mesh::CalculateFaceAreas() must be called prior "
			"to SimpleGrid::FromMeshFV()");
	}
	if (mesh.edgemap.size() == 0) {
		_EXCEPTIONT("Mesh::ConstructEdgeMap() must be called prior "
			"to SimpleGrid::FromMeshFV()");
	}

	// Number of faces
	size_t sFaces = mesh.faces.size();

	// Copy over areas
	m_dArea = mesh.vecFaceArea;

	// Generate connectivity from edge map
	std::vector< std::set<int> > m_vecConnectivitySet;
	m_vecConnectivitySet.resize(sFaces);

	size_t sOneSidedEdges = 0;
	EdgeMapConstIterator iterEdgeMap = mesh.edgemap.begin();
	for (; iterEdgeMap != mesh.edgemap.end(); iterEdgeMap++) {
		const FacePair & facepr = iterEdgeMap->second;
		if ((facepr[0] == InvalidFace) || (facepr[1] == InvalidFace)) {
			sOneSidedEdges++;
			continue;
		}
		if ((facepr[0] < 0) || (facepr[0] >= sFaces)) {
			_EXCEPTION1("EdgeMap FacePair out of range (%i)", facepr[0]);
		}
		if ((facepr[1] < 0) || (facepr[1] >= sFaces)) {
			_EXCEPTION1("EdgeMap FacePair out of range (%i)", facepr[1]);
		}
		m_vecConnectivitySet[facepr[0]].insert(facepr[1]);
		m_vecConnectivitySet[facepr[1]].insert(facepr[0]);
	}

	if (sOneSidedEdges != 0) {
		Announce("One-sided edges: %lu", sOneSidedEdges);
	}

	m_vecConnectivity.resize(sFaces);
	for (int i = 0; i < sFaces; i++) {
		std::set<int>::const_iterator iterConnect =
			m_vecConnectivitySet[i].begin();
		for (; iterConnect != m_vecConnectivitySet[i].end(); iterConnect++) {
			m_vecConnectivity[i].push_back(*iterConnect);
		}
	}

	// Generate centerpoints
	m_dLon.Allocate(sFaces);
	m_dLat.Allocate(sFaces);

	// Mesh is a RLL mesh -- special calculation for centerpoint values
	if (mesh.type == Mesh::MeshType_RLL) {
		_EXCEPTIONT("How to define nGridDim?");

		for (int i = 0; i < sFaces; i++) {

			const Face & face = mesh.faces[i];

			int nNodes = face.edges.size();

			if (nNodes != 4) {
				_EXCEPTIONT("RLL mesh must have exactly 4 nodes per face");
			}

			double dLonc = 0.0;
			double dLatc = 0.0;

			for (int j = 0; j < nNodes; j++) {

				const Node & node = mesh.nodes[face[j]];

				double dLon;
				double dLat;

				XYZtoRLL_Deg(
					node.x,
					node.y,
					node.z,
					dLon,
					dLat);

				// Consider periodicity of longitude
				if ((j != 0) && (fabs(dLonc / static_cast<double>(j) - dLon) > 180.0)) {
					if (dLonc > dLon) {
						dLon += 360.0;
					}
					if (dLonc < dLon) {
						dLon -= 360.0;
					}
				}

				// Sanity check: Longitude still out of range
				if ((j != 0) && (fabs(dLonc / static_cast<double>(j) - dLon) > 180.0)) {
					printf("FATAL:  Mesh face appears to extend more than 180 degrees longitude\n");
					for (int k = 0; k < j; k++) {

						const Node & nodeK = mesh.nodes[face[k]];

						XYZtoRLL_Deg(
							nodeK.x,
							nodeK.y,
							nodeK.z,
							dLon,
							dLat);

						printf("Node %i: %1.15e %1.15e\n", k, dLon, dLat);
					}
					_EXCEPTION();
				}

				dLonc += dLon;
				dLatc += dLat;
			}

			m_dLon[i] = dLonc / static_cast<double>(nNodes);
			m_dLat[i] = dLatc / static_cast<double>(nNodes);
		}

	// Any other kind of mesh, use Carteisan coords and project to surface of sphere
	} else {
		m_nGridDim.resize(1);
		m_nGridDim[0] = sFaces;

		for (size_t i = 0; i < sFaces; i++) {

			const Face & face = mesh.faces[i];

			int nNodes = face.edges.size();

			double dXc = 0.0;
			double dYc = 0.0;
			double dZc = 0.0;

			for (int j = 0; j < nNodes; j++) {
				const Node & node = mesh.nodes[face[j]];

				dXc += node.x;
				dYc += node.y;
				dZc += node.z;
			}

			dXc /= static_cast<double>(nNodes);
			dYc /= static_cast<double>(nNodes);
			dZc /= static_cast<double>(nNodes);

			double dMag = sqrt(dXc * dXc + dYc * dYc + dZc * dZc);

			dXc /= dMag;
			dYc /= dMag;
			dZc /= dMag;

			XYZtoRLL_Deg(dXc, dYc, dZc, m_dLon[i], m_dLat[i]);
		}
	}

	// Output total area
	{
		double dTotalArea = 0.0;
		for (size_t i = 0; i < sFaces; i++) {
			dTotalArea += m_dArea[i];
		}
		Announce("Total calculated area: %1.15e sr", dTotalArea);
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::FromMeshFE(
	const Mesh & mesh,
	bool fCGLL,
	int nP
) {
	if (!fCGLL) {
		_EXCEPTIONT("Sorry, not implemented yet!");
	}

	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call FromMeshFE() on previously initialized grid");
	}

	if (mesh.vecFaceArea.GetRows() == 0) {
		_EXCEPTIONT("Mesh::CalculateFaceAreas() must be called prior "
			"to SimpleGrid::FromMeshFE()");
	}
	if (mesh.edgemap.size() == 0) {
		_EXCEPTIONT("Mesh::ConstructEdgeMap() must be called prior "
			"to SimpleGrid::FromMeshFE()");
	}

	// Number of elements
	size_t nElements = mesh.faces.size();

	// Gauss-Lobatto nodes and weights
	DataArray1D<double> dG(nP);
	DataArray1D<double> dW(nP);

	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Generate coincident node map and Jacobian
	DataArray3D<int> dataGLLnodes(nP, nP, nElements);
	DataArray3D<double> dataGLLJacobian(nP, nP, nElements);

	GenerateMetaData(mesh, nP, true, dataGLLnodes, dataGLLJacobian);

	// Generate areas
	if (fCGLL) {
		GenerateUniqueJacobian(
			dataGLLnodes,
			dataGLLJacobian,
			m_dArea);
	} else {
		GenerateDiscontinuousJacobian(
			dataGLLJacobian,
			m_dArea);
	}

	// Number of faces
	size_t sFaces = m_dArea.GetRows();

	m_nGridDim.resize(1);
	m_nGridDim[0] = sFaces;

	// Generate coordinates
	{
		m_dLon.Allocate(sFaces);
		m_dLat.Allocate(sFaces);

		const NodeVector & nodevec = mesh.nodes;
		for (int k = 0; k < nElements; k++) {
			const Face & face = mesh.faces[k];

			if (face.edges.size() != 4) {
				_EXCEPTIONT("Mesh must only contain quadrilateral elements");
			}

			double dFaceNumericalArea = 0.0;

			for (int j = 0; j < nP; j++) {
			for (int i = 0; i < nP; i++) {
				int ix = dataGLLnodes[j][i][k]-1;

				_ASSERT((ix >= 0) && (ix < m_dLon.GetRows()));

				Node nodeGLL;
				Node dDx1G;
				Node dDx2G;

				ApplyLocalMap(
					face,
					nodevec,
					dG[i],
					dG[j],
					nodeGLL,
					dDx1G,
					dDx2G);

				XYZtoRLL_Deg(
					nodeGLL.x,
					nodeGLL.y,
					nodeGLL.z,
					m_dLon[ix],
					m_dLat[ix]);
			}
			}
		}
	}

	// Generate connectivity
	{
		std::vector< std::set<int> > vecConnectivitySet;
		vecConnectivitySet.resize(sFaces);

		for (size_t f = 0; f < nElements; f++) {

			for (int q = 0; q < nP; q++) {
			for (int p = 0; p < nP; p++) {

				std::set<int> & setLocalConnectivity =
					vecConnectivitySet[dataGLLnodes[q][p][f]-1];

				// Connect in all directions
				// Note that dataGLLnodes is 1-indexed, whereas we need 0-indexing for connectivity
				if (p != 0) {
					setLocalConnectivity.insert(
						dataGLLnodes[q][p-1][f]-1);
				}
				if (p != (nP-1)) {
					setLocalConnectivity.insert(
						dataGLLnodes[q][p+1][f]-1);
				}
				if (q != 0) {
					setLocalConnectivity.insert(
						dataGLLnodes[q-1][p][f]-1);
				}
				if (q != (nP-1)) {
					setLocalConnectivity.insert(
						dataGLLnodes[q+1][p][f]-1);
				}
			}
			}
		}

		m_vecConnectivity.resize(sFaces);
		for (size_t i = 0; i < sFaces; i++) {
			std::set<int>::const_iterator iterConnect =
				vecConnectivitySet[i].begin();
			for (; iterConnect != vecConnectivitySet[i].end(); iterConnect++) {
				m_vecConnectivity[i].push_back(*iterConnect);
			}
		}
	}

	// Output total area
	{
		double dTotalArea = 0.0;
		for (size_t i = 0; i < sFaces; i++) {
			dTotalArea += m_dArea[i];
		}
		Announce("Total calculated area: %1.15e", dTotalArea);
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::FromFile(
	const std::string & strConnectivityFile
) {
	if (IsInitialized()) {
		_EXCEPTIONT("Attempting to call FromFile() on previously initialized grid");
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
			if (m_vecConnectivity[f][n] == 0) {
				_EXCEPTION2("Zero index found in connectivity file \"%s\" for node %lu",
					strConnectivityFile.c_str(), f);
			}
			m_vecConnectivity[f][n]--;
			if (m_vecConnectivity[f][n] == f) {
				_EXCEPTION2("Self-connected node found in connectivity file \"%s\" for node %lu",
					strConnectivityFile.c_str(), f);
			}
			if (m_vecConnectivity[f][n] >= sFaces) {
				_EXCEPTION2("Out-of-range index found in connectivity file \"%s\" for node %lu",
					strConnectivityFile.c_str(), f);
			}
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
) const {
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
	fsOutput << std::endl;

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

	for (size_t i = 0; i < sFaces; i++) {
		fsOutput << m_dLon[i] << "," << m_dLat[i] << ","
			<< m_dArea[i] << "," << m_vecConnectivity[i].size();
		for (size_t j = 0; j < m_vecConnectivity[i].size(); j++) {
			fsOutput << "," << (m_vecConnectivity[i][j]+1);
		}
		fsOutput << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::FromUnstructuredDataFile(
	const std::string & strDataFile,
	const std::string & strLatitudeName,
	const std::string & strLongitudeName
) {
	NcFile ncfile(strDataFile.c_str());
	if (!ncfile.is_valid()) {
		_EXCEPTION1("Unable to open data file \"%s\"", strDataFile.c_str());
	}

	NcVar * varLon = ncfile.get_var(strLongitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION2("Unable to read longitude variable \"%s\" from data file \"%s\"",
			strLongitudeName.c_str(), strDataFile.c_str());
	}
	if (varLon->num_dims() != 1) {
		_EXCEPTION2("Longitude variable \"%s\" in data file \"%s\" must have dimension 1",
			strLongitudeName.c_str(), strDataFile.c_str());
	}

	NcVar * varLat = ncfile.get_var(strLatitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION2("Unable to read latitude variable \"%s\" from data file \"%s\"",
			strLatitudeName.c_str(), strDataFile.c_str());
	}
	if (varLat->num_dims() != 1) {
		_EXCEPTION2("Latitude variable \"%s\" in data file \"%s\" must have dimension 1",
			strLatitudeName.c_str(), strDataFile.c_str());
	}

	if (varLon->get_dim(0)->size() != varLat->get_dim(0)->size()) {
		_EXCEPTION4("Longitude variable \"%s\" (size %li) and latitude variable \"%s\" (size %li) must have same size",
			varLon->name(), varLon->get_dim(0)->size(),
			varLat->name(), varLat->get_dim(0)->size());
	}

	m_dLon.Allocate(varLon->get_dim(0)->size());
	m_dLat.Allocate(varLat->get_dim(0)->size());

	varLon->get(&(m_dLon[0]), varLon->get_dim(0)->size());
	varLat->get(&(m_dLat[0]), varLat->get_dim(0)->size());

	m_nGridDim.resize(1);
	m_nGridDim[0] = varLon->get_dim(0)->size();
}

///////////////////////////////////////////////////////////////////////////////

int SimpleGrid::CoordinateVectorToIndex(
	const std::vector<int> & coordvec
) const {
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

void SimpleGrid::BuildKDTree() {
	if (m_kdtree != NULL) {
		_EXCEPTIONT("kdtree already exists");
	}
	if (m_dLon.GetRows() == 0) {
		_EXCEPTIONT("At least one grid cell needed in SimpleGrid");
	}

	_ASSERT(m_dLon.GetRows() == m_dLat.GetRows());

	// Create the kd tree
	m_kdtree = kd_create(3);
	if (m_kdtree == NULL) {
		_EXCEPTIONT("kd_create(3) failed");
	}

	// Insert all nodes from this SimpleGrid into the kdtree
	for (size_t i = 0; i < m_dLon.GetRows(); i++) {
		double dX = cos(m_dLon[i]) * cos(m_dLat[i]);
		double dY = sin(m_dLon[i]) * cos(m_dLat[i]);
		double dZ = sin(m_dLat[i]);

		kd_insert3(m_kdtree, dX, dY, dZ, (void*)(i));
	}
}

///////////////////////////////////////////////////////////////////////////////

size_t SimpleGrid::NearestNode(
	double dLonRad,
	double dLatRad
) const {
	if (m_kdtree == NULL) {
		_EXCEPTIONT("BuildKDTree() must be called before NearestNode()");
	}

	// Find the nearest node from a given point
	double dX, dY, dZ;
	RLLtoXYZ_Rad(dLonRad, dLatRad, dX, dY, dZ);

	kdres * kdresNearest = kd_nearest3(m_kdtree, dX, dY, dZ);
	if (kdresNearest == NULL) {
		_EXCEPTIONT("kd_nearest3() failed");
	}
	int nResSize = kd_res_size(kdresNearest);
	if (nResSize != 1) {
		_EXCEPTION1("kd_nearest3() returned incorrect result size (%i)", nResSize);
	}

	size_t i = (size_t)(kd_res_item_data(kdresNearest));

	kd_res_free(kdresNearest);

	return i;
}

///////////////////////////////////////////////////////////////////////////////

void SimpleGrid::NearestNodes(
	double dLonRad,
	double dLatRad,
	double dDistDegGCD,
	std::vector<size_t> & vecNodeIxs
) const {
	if (m_kdtree == NULL) {
		_EXCEPTIONT("BuildKDTree() must be called before NearestNode()");
	}

	// Clear the index array
	vecNodeIxs.clear();

	// Find the nearest node from a given point
	double dX, dY, dZ;
	RLLtoXYZ_Rad(dLonRad, dLatRad, dX, dY, dZ);

	double dDistXYZ = 2.0 * sin(DegToRad(dDistDegGCD) / 2.0) + ReferenceTolerance;

	kdres * kdresNearestRange = kd_nearest_range3(m_kdtree, dX, dY, dZ, dDistXYZ);
	if (kdresNearestRange == NULL) {
		_EXCEPTIONT("kd_nearest_range3() failed");
	}
	int nResSize = kd_res_size(kdresNearestRange);
	if (nResSize != 0) {
		for (;;) {
			size_t i = (size_t)(kd_res_item_data(kdresNearestRange));
			vecNodeIxs.push_back(i);

			if (kd_res_next(kdresNearestRange) == 0) {
				break;
			}
		}
	}

	kd_res_free(kdresNearestRange);

}

///////////////////////////////////////////////////////////////////////////////

