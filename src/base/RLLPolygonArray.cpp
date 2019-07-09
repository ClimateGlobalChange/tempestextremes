///////////////////////////////////////////////////////////////////////////////
///
///	\file    RLLPolygonArray.h
///	\author  Paul Ullrich
///	\version July 2, 2019
///
///	<remarks>
///		Copyright 2000-2019 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "RLLPolygonArray.h"
#include "Exception.h"
#include "STLStringHelper.h"

#include <iostream>
#include <fstream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void RLLPolygonArray::FromFile(
	const std::string & strFilename
) {
	if ((m_vecNames.size() != 0) ||
	    (m_vecNodes.size() != 0)
	) {
		_EXCEPTIONT("RLLPolygonArray already initialized");
	}

	// String buffer
	std::string strBuffer;

	// Open the file as an input stream
	std::ifstream ifInput(strFilename);
	if (!ifInput.is_open()) {
		_EXCEPTION1("Unable to open input file \"%s\"", strFilename.c_str());
	}

	// Read mode
	enum {
		RLLPolygonArray_ReadModeCount,
		RLLPolygonArray_ReadModePolygons,
		RLLPolygonArray_ReadModeDone
	} eReadMode = RLLPolygonArray_ReadModeCount;

	// Count
	int nPolygonCount = (-1);

	// Number of polygons read
	int nPolygonRead = 0;

	// Loop through all lines
	int iLine = 0;
	for (;;) {
		iLine++;
		getline(ifInput, strBuffer);
		if (ifInput.eof()) {
			break;
		}

		if (strBuffer.length() == 0) {
			continue;
		}
		if (strBuffer[0] == '#') {
			continue;
		}

		// Check for end of file
		if (eReadMode == RLLPolygonArray_ReadModeDone) {
			_EXCEPTION1("Too many entries found in file \"%s\"",
				strFilename.c_str());
		}

		// Get the number of regions
		if (eReadMode == RLLPolygonArray_ReadModeCount) {
			nPolygonCount = atoi(strBuffer.c_str());
			if ((nPolygonCount <= 0) || (nPolygonCount > 100000)) {
				_EXCEPTION2("Sanity check failed: RLLPolygonArray file \"%s\" "
					"indicates %i polygons present (max 100000)",
					strFilename.c_str(), nPolygonCount);
			}
			m_vecNodes.resize(nPolygonCount);
			eReadMode = RLLPolygonArray_ReadModePolygons;
			continue;
		}

		// Read in a polygon
		if (eReadMode == RLLPolygonArray_ReadModePolygons) {

			enum {
				RLLPolygonArray_ReadModePolyInfo,
				RLLPolygonArray_ReadModePolyCount,
				RLLPolygonArray_ReadModePolyLon,
				RLLPolygonArray_ReadModePolyLat,
				RLLPolygonArray_ReadModePolyDone
			} eReadModePoly = RLLPolygonArray_ReadModePolyInfo;

			std::string strPolygonName;
			std::vector<double> vecPolygonLon;
			std::vector<double> vecPolygonLat;
			int nVertexCount = 0;
			int nLonRead = 0;
			int nLatRead = 0;
			int iLast = 0;
			for (int i = 0; i <= strBuffer.length(); i++) {

				if (eReadModePoly == RLLPolygonArray_ReadModePolyDone) {
					_EXCEPTION2("Spurious characters found on line %i of \"%s\"",
						iLine, strFilename.c_str());
				}

				// Check comma
				if ((i == strBuffer.length()) || (strBuffer[i] == ',')) {
					std::string strText = strBuffer.substr(iLast, i-iLast);
					iLast = i + 1;

					// Information about polygon
					if (eReadModePoly == RLLPolygonArray_ReadModePolyInfo) {
						strPolygonName = strText;
						eReadModePoly = RLLPolygonArray_ReadModePolyCount;
						continue;
					}

					// Number of vertices
					if (eReadModePoly == RLLPolygonArray_ReadModePolyCount) {
						nVertexCount = atoi(strText.c_str());
						if ((nVertexCount <= 0) || (nVertexCount > 1000)) {
							_EXCEPTION3("Sanity check failed: RLLPolygonArray file \"%s\" "
								"line %i indicates %i vertices present (max 1000)",
								strFilename.c_str(), iLine, nVertexCount);
						}
						eReadModePoly = RLLPolygonArray_ReadModePolyLon;
						continue;
					}

					// Longitude
					if (eReadModePoly == RLLPolygonArray_ReadModePolyLon) {
						if (!STLStringHelper::IsFloat(strText)) {
							_EXCEPTION3("Syntax error on line %i of \"%s\":"
								"Invalid float \"%s\"",
								iLine, strFilename.c_str(), strText.c_str());
						}
						vecPolygonLon.push_back(atof(strText.c_str()));

						if (vecPolygonLon.size() == nVertexCount) {
							eReadModePoly = RLLPolygonArray_ReadModePolyLat;
						}
						continue;
					}

					// Latitude
					if (eReadModePoly == RLLPolygonArray_ReadModePolyLat) {
						if (!STLStringHelper::IsFloat(strText)) {
							_EXCEPTION3("Syntax error on line %i of \"%s\":"
								"Invalid float \"%s\"",
								iLine, strFilename.c_str(), strText.c_str());
						}
						double dLat = atof(strText.c_str());
						vecPolygonLat.push_back(atof(strText.c_str()));

						if (fabs(dLat) > 90.0) {
							_EXCEPTION3("Error on line %i of \"%s\": "
								"Latitude of vertex out of range [-90,90] (%1.15e)",
								iLine, strFilename.c_str(), dLat);
						}

						if (vecPolygonLat.size() == nVertexCount) {
							eReadModePoly = RLLPolygonArray_ReadModePolyDone;
						}
						continue;
					}
				}

				// Check quotation
				if (strBuffer[i] == '\"') {
					if (eReadModePoly != RLLPolygonArray_ReadModePolyInfo) {
						_EXCEPTION2("Quotation marks only allowed in polygon name on line %i of \"%s\"",
							iLine, strFilename.c_str());
					}
					for (i++; i < strBuffer.length(); i++) {
						if (strBuffer[i] == '\"') {
							break;
						}
					}
					if (i == strBuffer.length()) {
						_EXCEPTION2("Unterminated quotation mark on line %i of \"%s\"",
							iLine, strFilename.c_str());
					}
				}
			}

			// Check polygon vertices
			if (vecPolygonLon.size() != nVertexCount) {
				_EXCEPTION2("Insufficient data on line %i of \"%s\"",
					iLine, strFilename.c_str());
			}
			if (vecPolygonLat.size() != nVertexCount) {
				_EXCEPTION2("Insufficient data on line %i of \"%s\"",
					iLine, strFilename.c_str());
			}

			// Check for line segments of zero length
			for (int v = 0; v < nVertexCount; v++) {
				int v1 = (v+1)%nVertexCount;

				double dLon0 = vecPolygonLon[v];
				double dLon1 = vecPolygonLon[v1];

				double dLat0 = vecPolygonLat[v];
				double dLat1 = vecPolygonLat[v1];

				if ((fabs(dLon0 - dLon1) < 1.0e-13) && (fabs(dLat0 - dLat1) < 1.0e-13)) {
					_EXCEPTION2("Line segment of zero length found on line %i of \"%s\"",
						iLine, strFilename.c_str());
				}
			}

			// Verify polygons don't exceed 360 degrees of longitude
			double dMinLon = vecPolygonLon[0];
			double dMaxLon = vecPolygonLon[0];
			for (int v = 1; v < nVertexCount; v++) {
				if (vecPolygonLon[v] < dMinLon) {
					dMinLon = vecPolygonLon[v];
				}
				if (vecPolygonLon[v] > dMaxLon) {
					dMaxLon = vecPolygonLon[v];
				}
			}
			if (dMaxLon - dMinLon > 360.0) {
				_EXCEPTION2("Polygons must not span more than 360 degrees of "
					"longitude on line %i of \"%s\"", iLine, strFilename.c_str());
			}

			// Insert polygon vertices into array
			int ixPoly = m_vecNames.size();

			m_vecNames.push_back(strPolygonName);
			m_vecNodes.push_back(RLLPointVector());
			for (int v = 0; v < nVertexCount; v++) {
				RLLPoint pt;
				pt.lon = vecPolygonLon[v];
				pt.lat = vecPolygonLat[v];
				m_vecNodes[ixPoly].push_back(pt);
			}

			if (m_vecNames.size() == nPolygonCount) {
				eReadMode = RLLPolygonArray_ReadModeDone;
			}
		}
	}

	if (m_vecNames.size() != nPolygonCount) {
		_EXCEPTION3("Number of polygons (%i) does not match count (%i) in \"%s\"",
			strFilename.c_str(), m_vecNames.size(), nPolygonCount);
	}
}

///////////////////////////////////////////////////////////////////////////////

const std::string & RLLPolygonArray::NameOfRegionContainingPoint(
	const RLLPoint & pt_in_degrees
) {
	static const double Threshold = 1.0e-13;

	const RLLPoint & pt = pt_in_degrees;

	_ASSERT(m_vecNames.size() != m_vecNodes.size());

	// Loop through all regions starting at the top until we find one
	// that matches.
	for (int r = 0; r < m_vecNames.size(); r++) {

		_ASSERT(m_vecNodes[r].size() > 0);

		// Minimum longitude of polygon
		double dMinLon = m_vecNodes[r][0].lon;
		double dMaxLon = m_vecNodes[r][0].lon;
		for (int i = 0; i < m_vecNodes[r].size(); i++) {
			if (m_vecNodes[r][i].lon < dMinLon) {
				dMinLon = m_vecNodes[r][i].lon;
			}
			if (m_vecNodes[r][i].lon > dMaxLon) {
				dMaxLon = m_vecNodes[r][i].lon;
			}
		}
		_ASSERT(dMaxLon - dMinLon <= 360.0);

		// If a point is inside a closed region, then a ray cast from the point
		// along a line of constant longitude will pass through an odd number
		// of edges.
		int nIntersections = 0;

		for (int i = 0; i < m_vecNodes[r].size(); i++) {

			int i1 = (i+1)%(m_vecNodes[r].size());

			// Move search to [dMinLon, dMinLon + 360.0]
			double dLon0 = pt.lon;
			double dLat0 = pt.lat;

			double dLon0Ix = floor((dLon0 - dMinLon) / 360.0);
			dLon0 = dLon0 - dLon0Ix * 360.0;
			_ASSERT((dLon0 >= dMinLon) && (dLon0 <= dMinLon + 360.0));
			dLon0 -= dMinLon;

			double dLon1 = m_vecNodes[r][i].lon - dMinLon;
			double dLat1 = m_vecNodes[r][i].lat;

			double dLon2 = m_vecNodes[r][i1].lon - dMinLon;
			double dLat2 = m_vecNodes[r][i1].lat;

			// Check for line segment of zero length
			if ((fabs(dLon1 - dLon2) < Threshold) &&
			    (fabs(dLat1 - dLat2) < Threshold)
			) {
				_EXCEPTIONT("Line segment of zero length detected");
			}

			// Check latitudes
			if (fabs(dLat0) > 90.0 + Threshold) {
				_EXCEPTION1("Latitude of point %1.5e out of range [-90,90]", dLat0);
			}
			if (fabs(dLat1) > 90.0 + Threshold) {
				_EXCEPTION1("Latitude of edge %1.5e out of range [-90,90]", dLat1);
			}
			if (fabs(dLat2) > 90.0 + Threshold) {
				_EXCEPTION1("Latitude of edge %1.5e out of range [-90,90]", dLat2);
			}

			// Edge of constant longitude
			if (fabs(dLon1 - dLon2) < Threshold) {

				// Check if the point lies along the edge
				if (fabs(dLon1 - dLon0) < Threshold) {
					double dT = (dLat0 - dLat1) / (dLat2 - dLat1);
					if ((dT > -Threshold) && (dT < 1.0 + Threshold)) {
						return m_vecNames[r];
					}
				}
				continue;

			// Not an edge of constant longitude
			} else {

				// Calculate longitude of intersection
				double dT = (dLon0 - dLon1) / (dLon2 - dLon1);

				// Latitude of intersection
				double dLatI = dLat1 + dT * (dLat2 - dLat1);

				// Check if point lies along edge
				if (fabs(dLatI - dLat0) < Threshold) {
					return m_vecNames[r];
				}

				// Check if ray intersects edge north of the point
				if ((dLatI > dLat0) && (dT > Threshold) && (dT < 1.0 - Threshold)) {
					nIntersections++;
				}
			}
		}

		// If the number of intersections is odd then the point is inside the polygon
		if ((nIntersections % 2) == 1) {
			return m_vecNames[r];
		}
	}

	return m_strDefault;
}

///////////////////////////////////////////////////////////////////////////////

