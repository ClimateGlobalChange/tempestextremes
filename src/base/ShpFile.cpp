///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShpFile.cpp
///	\author  Paul Ullrich
///	\version May 18, 2022
///
///	<remarks>
///		Copyright 2022 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"
#include "order32.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////

static const int32_t SHPFileCodeRef = 0x0000270a;

static const int32_t SHPVersionRef = 1000;

static const int32_t SHPPolygonType = 5;

struct SHPHeader {
	int32_t iFileCode;
	int32_t iUnused[5];
	int32_t iFileLength;
	int32_t iVersion;
	int32_t iShapeType;
};

struct SHPBounds {
	double dXmin;
	double dYmin;
	double dXmax;
	double dYmax;

	double dZmin;
	double dZmax;

	double dMmin;
	double dMmax;
};

struct SHPRecordHeader {
	int32_t iNumber;
	int32_t nLength;
};

struct SHPPolygonHeader {
	double dXmin;
	double dYmin;
	double dXmax;
	double dYmax;
	int32_t nNumParts;
	int32_t nNumPoints;
};

///////////////////////////////////////////////////////////////////////////////

int32_t SwapEndianInt32(const int32_t num) {

	int32_t res;
	const char * pnum = (const char *)(&num);
	char * pres = (char *)(&res);

	pres[0] = pnum[3];
	pres[1] = pnum[2];
	pres[2] = pnum[1];
	pres[3] = pnum[0];

	return res;
}

///////////////////////////////////////////////////////////////////////////////

double SwapEndianDouble(const double num) {

	double res;
  	const char * pnum = (const char *)(&num);
	char * pres = (char *)(&res);

	pres[0] = pnum[7];
	pres[1] = pnum[6];
	pres[2] = pnum[5];
	pres[3] = pnum[4];
	pres[4] = pnum[3];
	pres[5] = pnum[2];
	pres[6] = pnum[1];
	pres[7] = pnum[0];

	return res;
}

///////////////////////////////////////////////////////////////////////////////

void ReadShpFileAsMesh(
	const std::string & strInputFile,
	Mesh & mesh,
	bool fVerbose
) {
	// Units for each dimension
	const std::string strXYUnits = "lonlat";

	// Convexify the mesh
	bool fConvexify = false;

	// First polygon
	int iPolygonFirst = (-1);

	// Last polygon
	int iPolygonLast = (-1);

	// Coarsen
	int nCoarsen = 1;

	// Check file names
	if (strInputFile == "") {
		_EXCEPTIONT("No input file specified");
	}
	if (nCoarsen < 1) {
		_EXCEPTIONT("--coarsen must be greater than or equal to 1");
	}

	// Clear the mesh
	mesh.Clear();

	// Load shapefile
	if (fVerbose) AnnounceStartBlock("Loading SHP file \"%s\"", strInputFile.c_str());
	std::ifstream shpfile(
		strInputFile.c_str(), std::ios::in | std::ios::binary);
	if (!shpfile.is_open()) {
		_EXCEPTION1("Unable to open SHP file \"%s\"", strInputFile.c_str());
	}

	SHPHeader shphead;
	shpfile.read((char*)(&shphead), sizeof(SHPHeader));

	if (O32_HOST_ORDER == O32_LITTLE_ENDIAN) {
		shphead.iFileCode = SwapEndianInt32(shphead.iFileCode);
		shphead.iFileLength = SwapEndianInt32(shphead.iFileLength);
	} else if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
		shphead.iVersion = SwapEndianInt32(shphead.iVersion);
		shphead.iShapeType = SwapEndianInt32(shphead.iShapeType);
	} else {
		_EXCEPTIONT("Invalid system Endian");
	}

	if (shphead.iFileCode != SHPFileCodeRef) {
		_EXCEPTIONT("Input file does not appear to be a ESRI Shapefile: "
			"File code mismatch");
	}
	if (shphead.iVersion != SHPVersionRef) {
		_EXCEPTIONT("Input file error: Version mismatch");
	}
	if (shphead.iShapeType != SHPPolygonType) {
		_EXCEPTIONT("Input file error: Polygon type expected");
	}

	SHPBounds shpbounds;
	shpfile.read((char*)(&shpbounds), sizeof(SHPBounds));

	if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
		shpbounds.dXmin = SwapEndianDouble(shpbounds.dXmin);
		shpbounds.dYmin = SwapEndianDouble(shpbounds.dYmin);
		shpbounds.dXmax = SwapEndianDouble(shpbounds.dXmax);
		shpbounds.dYmax = SwapEndianDouble(shpbounds.dXmax);

		shpbounds.dZmin = SwapEndianDouble(shpbounds.dZmin);
		shpbounds.dZmax = SwapEndianDouble(shpbounds.dZmax);

		shpbounds.dMmin = SwapEndianDouble(shpbounds.dMmin);
		shpbounds.dMmax = SwapEndianDouble(shpbounds.dMmax);
	}

	// Current position (in 16-bit words)
	int32_t iCurrentPosition = 50;

	int32_t iPolygonIx = 1;

	// Load records
	while (iCurrentPosition < shphead.iFileLength) {

		// Read the record header
		SHPRecordHeader shprechead;
		shpfile.read((char*)(&shprechead), sizeof(SHPRecordHeader));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_LITTLE_ENDIAN) {
			shprechead.iNumber = SwapEndianInt32(shprechead.iNumber);
			shprechead.nLength = SwapEndianInt32(shprechead.nLength);
		}

		iPolygonIx++;

		if (fVerbose) AnnounceStartBlock("Polygon %i", shprechead.iNumber);

		iCurrentPosition += shprechead.nLength;

		// Read the shape type
		int32_t iShapeType;
		shpfile.read((char*)(&iShapeType), sizeof(int32_t));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			iShapeType = SwapEndianInt32(iShapeType);
		}
		if (iShapeType != SHPPolygonType) {
			_EXCEPTIONT("Input file error: Record Polygon type expected");
		}

		// Read the polygon header
		SHPPolygonHeader shppolyhead;
		shpfile.read((char*)(&shppolyhead), sizeof(SHPPolygonHeader));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			shppolyhead.dXmin = SwapEndianDouble(shppolyhead.dXmin);
			shppolyhead.dYmin = SwapEndianDouble(shppolyhead.dYmin);
			shppolyhead.dXmax = SwapEndianDouble(shppolyhead.dXmax);
			shppolyhead.dYmax = SwapEndianDouble(shppolyhead.dYmax);
			shppolyhead.nNumParts = SwapEndianInt32(shppolyhead.nNumParts);
			shppolyhead.nNumPoints = SwapEndianInt32(shppolyhead.nNumPoints);
		}

		// Sanity check
		if (shppolyhead.nNumParts > 0x1000000) {
			_EXCEPTION1("Polygon NumParts exceeds sanity bound (%i)",
				shppolyhead.nNumParts);
		}
		if (shppolyhead.nNumPoints > 0x1000000) {
			_EXCEPTION1("Polygon NumPoints exceeds sanity bound (%i)",
				shppolyhead.nNumPoints);
		}
		if (fVerbose) {
			Announce("containing %i part(s) with %i points",
				shppolyhead.nNumParts,
				shppolyhead.nNumPoints);
			Announce("Xmin: %3.5f", shppolyhead.dXmin);
			Announce("Ymin: %3.5f", shppolyhead.dYmin);
			Announce("Xmax: %3.5f", shppolyhead.dXmax);
			Announce("Ymax: %3.5f", shppolyhead.dYmax);
		}

		if (shppolyhead.nNumParts != 1) {
			Announce("WARNING: Only polygons with 1 part currently supported"
				" in Exodus format; ignoring remaining parts");
			//_EXCEPTIONT("Only polygons with 1 part currently supported"
			//	" in Exodus format");
		}

		DataArray1D<int32_t> iPartBeginIx(shppolyhead.nNumParts+1);
		shpfile.read((char*)&(iPartBeginIx[0]),
			shppolyhead.nNumParts * sizeof(int32_t));
		if (shpfile.eof()) {
			break;
		}
		iPartBeginIx[shppolyhead.nNumParts] = shppolyhead.nNumPoints;

		DataArray1D<double> dPoints(shppolyhead.nNumPoints * 2);
		shpfile.read((char*)&(dPoints[0]),
			shppolyhead.nNumPoints * 2 * sizeof(double));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			for (int i = 0; i < shppolyhead.nNumParts; i++) {
				iPartBeginIx[i] = SwapEndianInt32(iPartBeginIx[i]);
			}
			for (int i = 0; i < shppolyhead.nNumPoints * 2; i++) {
				dPoints[i] = SwapEndianDouble(dPoints[i]);
			}
		}

		if ((iPolygonFirst != (-1)) && (shprechead.iNumber < iPolygonFirst)) {
			AnnounceEndBlock("Done");
			continue;
		}

		if ((iPolygonLast != (-1)) && (shprechead.iNumber > iPolygonLast)) {
			AnnounceEndBlock("Done");
			continue;
		}

		// Find the largest part
		int iLargestPartBeginIx = 0;
		int iLargestPartEndIx = shppolyhead.nNumPoints;
		if (shppolyhead.nNumParts > 1) {
			int iMaxPartIx = 0;
			int nMaxPartSize = 0;
			for (int i = 0; i < iPartBeginIx.GetRows()-1; i++) {
				int nPartSize = iPartBeginIx[i+1] - iPartBeginIx[i];
				if (nPartSize > nMaxPartSize) {
					iMaxPartIx = i;
					nMaxPartSize = nPartSize;
					iLargestPartBeginIx = iPartBeginIx[i];
					iLargestPartEndIx = iPartBeginIx[i+1];
				}
			}
			if (fVerbose) {
				Announce("Using part %i with %i points [%i, %i]",
					iMaxPartIx, nMaxPartSize,
					iLargestPartBeginIx, iLargestPartEndIx-1);
			}
		}
/*
		FILE * fp = fopen("log.txt", "w");
		for (int i = 0; i < shppolyhead.nNumPoints; i++) {
			fprintf(fp, "%1.15e %1.15e\n", dPoints[2*i], dPoints[2*i+1]);
		}
		fclose(fp);
*/
		// Convert to Mesh
		int nFaces = mesh.faces.size();
		int nNodes = mesh.nodes.size();

		int nShpNodes = (iLargestPartEndIx - iLargestPartBeginIx) / nCoarsen;

		mesh.faces.resize(nFaces+1);
		mesh.nodes.resize(nNodes + nShpNodes);

		mesh.faces[nFaces] = Face(nShpNodes);

		// Convert from longitude/latitude to XYZ
		if (strXYUnits == "lonlat") {
			for (int i = 0; i < nShpNodes; i++) {
				int ix = iLargestPartBeginIx + i * nCoarsen;

				double dLonRad = dPoints[2*ix] / 180.0 * M_PI;
				double dLatRad = dPoints[2*ix+1] / 180.0 * M_PI;

				mesh.nodes[nNodes+i].x = cos(dLatRad) * cos(dLonRad);
				mesh.nodes[nNodes+i].y = cos(dLatRad) * sin(dLonRad);
				mesh.nodes[nNodes+i].z = sin(dLatRad);

				// Note that shapefile polygons are specified in clockwise order,
				// whereas Exodus files request polygons to be specified in
				// counter-clockwise order.  Hence we need to reorient these Faces.
				mesh.faces[nFaces].SetNode(
					nShpNodes - i - 1, nNodes + i);
			}

		} else {
			_EXCEPTION1("Invalid units \"%s\"", strXYUnits.c_str());
		}

		if (fVerbose) AnnounceEndBlock("Done");
	}

	if (fVerbose) AnnounceEndBlock("Done");
}

///////////////////////////////////////////////////////////////////////////////



