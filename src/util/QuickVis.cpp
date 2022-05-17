///////////////////////////////////////////////////////////////////////////////
///
///	\file    QuickVis.cpp
///	\author  Paul Ullrich
///	\version May 13, 2022
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
#include "Exception.h"
#include "Announce.h"
#include "CoordTransforms.h"
#include "NetCDFUtilities.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "GridElements.h"
#include "Variable.h"
#include "NcFileVector.h"
#include "STLStringHelper.h"

#include "netcdfcpp.h"
#include "kdtree.h"
#include "lodepng.h"


///////////////////////////////////////////////////////////////////////////////

void GenerateMap(
	const std::string & strInputGrid,
	const std::string & strInputConnect,
	const std::string & strInputData,
	const std::string & strLonName,
	const std::string & strLatName,
	const DataArray1D<double> & dSampleLon,
	const DataArray1D<double> & dSampleLat,
	DataArray2D<int> & dImageMap
) {

	int nResX = dSampleLon.GetRows();
	int nResY = dSampleLat.GetRows();

	// Allocate the image map
	dImageMap.Allocate(nResY, nResX);

	// Reference integer for memory indexing
	int iRef = 0;

	// Load the data file, if provided
	NcFile * pncdata = NULL;
	if (strInputData.length() != 0) {
		pncdata = new NcFile(strInputData.c_str());
		if (!pncdata->is_valid()) {
			_EXCEPTION1("Unable to open data file \"%s\"", strInputData.c_str());
		}
	}

	// Build the kdtree
	kdtree * kd = kd_create(3);
	if (kd == NULL) {
		_EXCEPTIONT("kd_create(3) failed");
	}

	// Generate kdtree from grid
	if (strInputGrid.length() != 0) {
		AnnounceStartBlock("Generating kdtree from gridfile");

		Mesh meshInput(strInputGrid);

		Announce("Gridfile contains %lu faces", meshInput.faces.size());

		for (int i = 0; i < meshInput.faces.size(); i++) {
			const Face & face = meshInput.faces[i];

			Node nodeCentroid(0.0, 0.0, 0.0);

			for (int j = 0; j < face.edges.size(); j++) {
				const Node & nodeVertex = meshInput.nodes[face[j]];
				nodeCentroid.x += nodeVertex.x;
				nodeCentroid.y += nodeVertex.y;
				nodeCentroid.z += nodeVertex.z;
			}

			nodeCentroid /= nodeCentroid.Magnitude();

			kd_insert3(kd, nodeCentroid.x, nodeCentroid.y, nodeCentroid.z, (void*)(&iRef + i));
		}

		AnnounceEndBlock("Done");

	// Generate kdtree from datafile
	} else if (strInputData.length() != 0) {
		_ASSERT(pncdata != NULL);

		AnnounceStartBlock("Generating kdtree from datafile");

		NcVar * varLon = pncdata->get_var(strLonName.c_str());
		if (varLon == NULL) {
			_EXCEPTION2("Unable to read longitude variable \"%s\" from data file \"%s\"",
				strLonName.c_str(), strInputData.c_str());
		}
		if (varLon->num_dims() != 1) {
			_EXCEPTION2("Longitude variable \"%s\" in data file \"%s\" must have dimension 1",
				strLonName.c_str(), strInputData.c_str());
		}

		NcVar * varLat = pncdata->get_var(strLatName.c_str());
		if (varLon == NULL) {
			_EXCEPTION2("Unable to read latitude variable \"%s\" from data file \"%s\"",
				strLatName.c_str(), strInputData.c_str());
		}
		if (varLat->num_dims() != 1) {
			_EXCEPTION2("Latitude variable \"%s\" in data file \"%s\" must have dimension 1",
				strLatName.c_str(), strInputData.c_str());
		}

		if (varLon->get_dim(0)->size() != varLat->get_dim(0)->size()) {
			_EXCEPTION4("Longitude variable \"%s\" (size %li) and latitude variable \"%s\" (size %li) must have same size",
				varLon->name(), varLon->get_dim(0)->size(),
				varLat->name(), varLat->get_dim(0)->size());
		}

		DataArray1D<double> dLon(varLon->get_dim(0)->size());
		DataArray1D<double> dLat(varLat->get_dim(0)->size());

		varLon->get(&(dLon[0]), varLon->get_dim(0)->size());
		varLat->get(&(dLat[0]), varLat->get_dim(0)->size());

		long iReportSize = varLon->get_dim(0)->size() / 100;
		for (long i = 0; i < varLon->get_dim(0)->size(); i++) {
			double dX;
			double dY;
			double dZ;

			RLLtoXYZ_Deg(dLon[i], dLat[i], dX, dY, dZ);

			kd_insert3(kd, dX, dY, dZ, (void*)(&iRef + i));

			if (i % iReportSize == 0) {
				Announce("%li%% complete", i / iReportSize);
			}
		}

		AnnounceEndBlock("Done");

	} else {
		_EXCEPTIONT("Not implemented");
	}

	// Query the kdtree to build image map
	{
		long lTotal = 0;
		long lReportSize = static_cast<long>(dSampleLat.GetRows()) * static_cast<long>(dSampleLon.GetRows()) / 100;
		AnnounceStartBlock("Querying data points within the kdtree");
		for (int j = 0; j < dSampleLat.GetRows(); j++) {
		for (int i = 0; i < dSampleLon.GetRows(); i++) {
			double dX;
			double dY;
			double dZ;

			RLLtoXYZ_Deg(dSampleLon[i], dSampleLat[j], dX, dY, dZ);

			kdres * kdres = kd_nearest3(kd, dX, dY, dZ);
			if (kdres == NULL) {
				_EXCEPTIONT("Unexpected query failure in kdtree");
			}

			dImageMap[j][i] = static_cast<int>(reinterpret_cast<int*>(kd_res_item_data(kdres)) - &iRef);

			lTotal++;
			if (lTotal % lReportSize == 0) {
				Announce("%li%% complete", lTotal / lReportSize);
			}
		}
		}
		AnnounceEndBlock("Done");
	}

	// Clean up
	kd_free(kd);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input grid
	std::string strInputGrid;

	// Input connectivity file
	std::string strInputConnect;

	// Input map
	std::string strInputMap;

	// Input data
	std::string strInputData;

	// Output PNG
	std::string strOutputPNG;

	// Output map
	std::string strOutputMap;

	// Variable to plot
	std::string strVarName;

	// Time
	int iTime;

	// X resolution of output image
	int nResX;

	// Y resolution of output image
	int nResY;

	// Minimum longitude of plot
	double dMinLon;

	// Maximum longitude of plot
	double dMaxLon;

	// Minimum latitude of plot
	double dMinLat;

	// Maximum latitude of plot
	double dMaxLat;

	// Minimum value in color range
	std::string strMinRange;

	// Maximum value in color range
	std::string strMaxRange;

	// Name of longitude variable
	std::string strLonName;

	// Name of latitude variable
	std::string strLatName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputGrid, "in_grid", "");
		CommandLineString(strInputConnect, "in_connect", "");
		CommandLineString(strInputMap, "in_map", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strOutputPNG, "out_png", "");
		CommandLineString(strOutputMap, "out_map", "");
		CommandLineString(strVarName, "var", "");
		CommandLineInt(iTime, "time", 0);
		CommandLineInt(nResX, "nx", 720);
		CommandLineInt(nResY, "ny", 360);
		CommandLineDouble(dMinLon, "minlon", 0.0);
		CommandLineDouble(dMaxLon, "maxlon", 360.0);
		CommandLineDouble(dMinLat, "minlat", -90.0);
		CommandLineDouble(dMaxLat, "maxlat", 90.0);
		CommandLineString(strMinRange, "minrange", "");
		CommandLineString(strMaxRange, "maxrange", "");
		CommandLineString(strLonName, "lonname", "lon");
		CommandLineString(strLatName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Validate arguments
	int nInputArguments =
		  ((strInputGrid.length() != 0)?(1):(0))
		+ ((strInputConnect.length() != 0)?(1):(0))
		+ ((strInputMap.length() != 0)?(1):(0));
	if (nInputArguments > 1) {
		_EXCEPTIONT("Only one of --in_grid, --in_connect, or --in_map can be specified");
	}

	if ((strInputMap.length() != 0) && (strOutputMap.length() != 0)) {
		_EXCEPTIONT("Only one of --in_map or --out_map can be specified");
	}

	// Variable registry
	VariableRegistry varreg;
	int ixVariable = varreg.FindOrRegister(strVarName);

	// Load in the benchmark file
	NcFileVector vecFiles;

	// Define the SimpleGrid for the input data
	SimpleGrid grid;

	if (strInputConnect != "") {
		grid.FromFile(strInputConnect);
	} else if (strInputData != "") {
		grid.FromUnstructuredDataFile(strInputData, strLatName, strLonName);
	}

	// Image map
	DataArray2D<int> dImageMap;

	DataArray1D<double> dSampleLon;
	DataArray1D<double> dSampleLat;

	// Build image map
	if (strInputMap == "") {
		dSampleLon.Allocate(nResX);
		for (int i = 0; i < nResX; i++) {
			dSampleLon[i] = dMinLon + (dMaxLon - dMinLon) * ((static_cast<double>(i) + 0.5) / static_cast<double>(nResX));
		}

		dSampleLat.Allocate(nResY);
		for (int i = 0; i < nResY; i++) {
			dSampleLat[i] = dMinLat + (dMaxLat - dMinLat) * ((static_cast<double>(i) + 0.5) / static_cast<double>(nResY));
		}

		GenerateMap(
			strInputGrid,
			strInputConnect,
			strInputData,
			strLonName,
			strLatName,
			dSampleLon,
			dSampleLat,
			dImageMap);

	// Read image map
	} else {
		AnnounceStartBlock("Reading image map");

		NcFile ncMapFile(strInputMap.c_str());
		if (!ncMapFile.is_valid()) {
			_EXCEPTION1("Unable to open output map file \"%s\"",
				strOutputMap.c_str());
		}

		NcVar * varLon = ncMapFile.get_var("lon");
		if (varLon == NULL) {
			_EXCEPTION1("Unable to read variable \"lon\" in input map file \"%s\"",
				strOutputMap.c_str());
		}
		if (varLon->num_dims() != 1) {
			_EXCEPTION1("Variable \"lon\" in input map file \"%s\" must have dimension 1",
				strOutputMap.c_str());
		}
		nResX = varLon->get_dim(0)->size();
		dSampleLon.Allocate(nResX);

		NcVar * varLat = ncMapFile.get_var("lat");
		if (varLat == NULL) {
			_EXCEPTION1("Unable to read variable \"lat\" in input map file \"%s\"",
				strOutputMap.c_str());
		}
		if (varLat->num_dims() != 1) {
			_EXCEPTION1("Variable \"lat\" in input map file \"%s\" must have dimension 1",
				strOutputMap.c_str());
		}
		nResY = varLat->get_dim(0)->size();
		dSampleLat.Allocate(nResY);

		NcVar * varMap = ncMapFile.get_var("image_map");
		if (varMap == NULL) {
			_EXCEPTION1("Unable to read variable \"image_map\" in input map file \"%s\"",
				strOutputMap.c_str());
		}
		if (varMap->num_dims() != 2) {
			_EXCEPTION1("Variable \"image_map\" in input map file \"%s\" must have dimension 2",
				strOutputMap.c_str());
		}
		if (varMap->get_dim(0)->size() != varLat->get_dim(0)->size()) {
			_EXCEPTION3("Variable \"image_map\" in input map file \"%s\" length of first dimension (%li) "
				"must match length of variable \"lat\" (%li)",
				strOutputMap.c_str(), varMap->get_dim(0)->size(), varLat->get_dim(0)->size());
		}
		if (varMap->get_dim(1)->size() != varLon->get_dim(0)->size()) {
			_EXCEPTION3("Variable \"image_map\" in input map file \"%s\" length of second dimension (%li) "
				"must match length of variable \"lon\" (%li)",
				strOutputMap.c_str(), varMap->get_dim(1)->size(), varLon->get_dim(0)->size());
		}
		dImageMap.Allocate(varLat->get_dim(0)->size(), varLon->get_dim(0)->size());
		varMap->get(&(dImageMap[0][0]), nResY, nResX);

		AnnounceEndBlock("Done");
	}

	// Write image map
	if (strOutputMap != "") {
		AnnounceStartBlock("Writing image map");

		NcFile ncMapFile(strOutputMap.c_str(), NcFile::Replace);
		if (!ncMapFile.is_valid()) {
			_EXCEPTION1("Unable to open output map file \"%s\"",
				strOutputMap.c_str());
		}

		NcDim * dimLon = ncMapFile.add_dim("lon", nResX);
		if (dimLon == NULL) {
			_EXCEPTION1("Unable to add dimension \"lon\" to output map file \"%s\"",
				strOutputMap.c_str());
		}

		NcDim * dimLat = ncMapFile.add_dim("lat", nResY);
		if (dimLat == NULL) {
			_EXCEPTION1("Unable to add dimension \"lat\" to output map file \"%s\"",
				strOutputMap.c_str());
		}

		NcVar * varLon = ncMapFile.add_var("lon", ncDouble, dimLon);
		if (varLon == NULL) {
			_EXCEPTION1("Unable to add variable \"lon\" to output map file \"%s\"",
				strOutputMap.c_str());
		}
		varLon->put(&(dSampleLon[0]), dimLon->size());

		NcVar * varLat = ncMapFile.add_var("lat", ncDouble, dimLat);
		if (varLon == NULL) {
			_EXCEPTION1("Unable to add variable \"lat\" to output map file \"%s\"",
				strOutputMap.c_str());
		}
		varLat->put(&(dSampleLat[0]), dimLat->size());

		NcVar * varMap = ncMapFile.add_var("image_map", ncInt, dimLat, dimLon);
		if (varMap == NULL) {
			_EXCEPTION1("Unable to add variable \"image_map\" to output map file \"%s\"",
				strOutputMap.c_str());
		}
		varMap->put(&(dImageMap[0][0]), dimLat->size(), dimLon->size());

		AnnounceEndBlock("Done");
	}

	// Compose image
	if (strOutputPNG != "") {
		AnnounceStartBlock("Producing image");

		// Archive input data filename
		vecFiles.ParseFromString(strInputData);

		// Get time information
		const NcTimeDimension & vecTimes = vecFiles.GetNcTimeDimension(0);
		if (vecTimes.size() == 0) {
			_EXCEPTION1("Input file \"%s\" contain no time information (%s)",
				 strInputData.c_str());
		}

		if (iTime >= vecTimes.size()) {
			_EXCEPTION2("Time index (%i) out of range (< %lu)", iTime, vecTimes.size());
		}

		// Load grid data
		Variable & var = varreg.Get(ixVariable);
		vecFiles.SetTime(vecTimes[iTime]);
		var.LoadGridData(varreg, vecFiles, grid);
		const DataArray1D<float> & data = var.GetData();

		_ASSERT(data.GetRows() > 0);

		// Get minimum range
		double dMinRange = 0.0;
		double dMaxRange = 0.0;

		bool fIsMinRangeGiven = STLStringHelper::IsFloat(strMinRange);
		bool fIsMaxRangeGiven = STLStringHelper::IsFloat(strMaxRange);

		if (fIsMinRangeGiven && fIsMaxRangeGiven) {
			dMinRange = stof(strMinRange);
			dMaxRange = stof(strMaxRange);

		} else if (fIsMinRangeGiven) {
			dMinRange = stof(strMinRange);
			dMaxRange = data[0];
			for (int i = 1; i < data.GetRows(); i++) {
				if (data[i] > dMaxRange) {
					dMaxRange = data[i];
				}
			}

		} else if (fIsMaxRangeGiven) {
			dMaxRange = stof(strMaxRange);
			dMinRange = data[0];
			for (int i = 1; i < data.GetRows(); i++) {
				if (data[i] < dMinRange) {
					dMinRange = data[i];
				}
			}

		} else {
			dMinRange = data[0];
			dMaxRange = data[0];
			for (int i = 1; i < data.GetRows(); i++) {
				if (data[i] < dMinRange) {
					dMinRange = data[i];
				}
				if (data[i] > dMaxRange) {
					dMaxRange = data[i];
				}
			}
		}
		
		Announce("Colorbar spans range [%f, %f]", dMinRange, dMaxRange);

		// Build image
		std::vector<unsigned char> image;
		image.resize(nResX * nResY * 4);
		for (int j = 0; j < nResY; j++) {
		for (int i = 0; i < nResX; i++) {
			int jx = nResY - j - 1;
			double dFrac = 255.0 * (data[dImageMap[j][i]] - dMinRange) / (dMaxRange - dMinRange);
			image[4 * nResX * jx + 4 * i + 0] = static_cast<unsigned char>(dFrac);
			image[4 * nResX * jx + 4 * i + 1] = static_cast<unsigned char>(dFrac);
			image[4 * nResX * jx + 4 * i + 2] = static_cast<unsigned char>(dFrac);
			image[4 * nResX * jx + 4 * i + 3] = 255;
		}
		}

		AnnounceEndBlock(NULL);

		AnnounceStartBlock("Encoding PNG");

		unsigned error = lodepng::encode(strOutputPNG, image, nResX, nResY);
		if (error) {
			_EXCEPTION2("Encoder error (%u): %s", error, lodepng_error_text(error));
		}

		AnnounceEndBlock("Done");
	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

