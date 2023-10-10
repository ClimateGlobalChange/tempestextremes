///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShapefileMask.cpp
///	\author  Paul Ullrich
///	\version August 29, 2023
///
///	<remarks>
///		Copyright 2023 Paul Ullrich
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
#include "NetCDFUtilities.h"
#include "LatLonBox.h"
#include "GridElements.h"
#include "ShpFile.h"
#include "STLStringHelper.h"

#include "netcdfcpp.h"

#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

bool FaceContainsNode(
	const Face & face,
	const NodeVector & nodevec,
	double dLatDeg,
	double dLonDeg
) {
	Node n0;
	RLLtoXYZ_Deg(dLonDeg, dLatDeg, n0.x, n0.y, n0.z);

	return face.Contains(n0, nodevec);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Source data
	std::string strSourceData;

	// Shapefile
	std::string strShapefile;

	// Output data
	std::string strOutputData;

	// Variable name
	std::string strVariable;

	// Longitude name
	std::string strLongitudeName;

	// Latitude name
	std::string strLatitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strSourceData, "in_data", "");
		CommandLineString(strShapefile, "shp", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check arguments
	if (strSourceData.length() == 0) {
		_EXCEPTIONT("No input file (--in_data) specified");
	}
	if (strShapefile.length() == 0) {
		_EXCEPTIONT("No shapefile (--shp) specified");
	}
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if (strVariable.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}

	// Parse variable list (--var)
	std::vector<std::string> vecVariableStrings;
	STLStringHelper::ParseVariableList(strVariable, vecVariableStrings);

	std::vector<std::string> vecVariableNames;
	std::vector< std::vector<std::string> > vecVariableSpecifiedDims;

	STLStringHelper::SplitVariableStrings(
		vecVariableStrings,
		vecVariableNames,
		vecVariableSpecifiedDims
	);

	// Load the shapefile
	AnnounceStartBlock("Loading shapefile");
	Mesh mesh;
	ReadShpFileAsMesh(strShapefile, mesh);
	AnnounceEndBlock("Done");

	// Number of regions
	const size_t sShpRegionCount = mesh.faces.size();

	// Build latlon boxes for each shapefile
	AnnounceStartBlock("Generating latlon boxes for each region");
	std::vector< LatLonBox<double> > vecLatLonBox(sShpRegionCount);
	for (int f = 0; f < sShpRegionCount; f++) {
		Face & face = mesh.faces[f];

		for (int i = 0; i < face.edges.size(); i++) {

			Node & node = mesh.nodes[face[i]];

			double dFaceLonDeg, dFaceLatDeg;
			XYZtoRLL_Deg(node.x, node.y, node.z, dFaceLonDeg, dFaceLatDeg);

			dFaceLonDeg = LonDegToStandardRange(dFaceLonDeg);

			vecLatLonBox[f].insert(dFaceLatDeg, dFaceLonDeg);
		}
	}
	AnnounceEndBlock("Done");

	// Initialize input and output files
	AnnounceStartBlock("Initializing files");

	NcFile ncInput(strSourceData.c_str());
	if (!ncInput.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\"", strSourceData.c_str());
	}

	std::vector<double> dLonDeg;
	std::vector<double> dLatDeg;

	// Get longitude data
	NcVar * varLon = ncInput.get_var(strLongitudeName.c_str());
	if (varLon == NULL) {
		_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
			strSourceData.c_str(), strLongitudeName.c_str());
	}
	if (varLon->num_dims() != 1) {
		_EXCEPTION2("File \"%s\" variable \"%s\" must contain one dimension",
			strSourceData.c_str(), strLongitudeName.c_str());
	}
	dLonDeg.resize(varLon->get_dim(0)->size());
	varLon->get(&(dLonDeg[0]), dLonDeg.size());

	// Get latitude data
	NcVar * varLat = ncInput.get_var(strLatitudeName.c_str());
	if (varLat == NULL) {
		_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
			strSourceData.c_str(), strLatitudeName.c_str());
	}
	if (varLat->num_dims() != 1) {
		_EXCEPTION2("File \"%s\" variable \"%s\" must contain one dimension",
			strSourceData.c_str(), strLatitudeName.c_str());
	}
	dLatDeg.resize(varLat->get_dim(0)->size());
	varLat->get(&(dLatDeg[0]), dLatDeg.size());

	// Open output file
	NcFile ncOutput(strOutputData.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputData.c_str());
	}

	// Create a new variable
	NcDim * dimShp = ncOutput.add_dim("shpix", sShpRegionCount);
	if (dimShp == NULL) {
		_EXCEPTIONT("Error creating dimension \"shpix\" in output file");
	}

	AnnounceEndBlock("Done");

	// Map from input data to shapefiles
	std::vector<size_t> sShpMap;
	std::vector<size_t> sDataCount(sShpRegionCount, 0);

	// Loop through all variables
	AnnounceStartBlock("Processing data");

	for (size_t v = 0; v < vecVariableStrings.size(); v++) {
		AnnounceStartBlock("Variable \"%s\"", vecVariableStrings[v].c_str());

		// Get data values
		NcVar * varData = ncInput.get_var(vecVariableStrings[v].c_str());
		if (varData == NULL) {
			_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
				strSourceData.c_str(), vecVariableStrings[v].c_str());
		}
		if ((varData->num_dims() < 2) || (varData->num_dims() > 3)) {
			_EXCEPTION2("File \"%s\" variable \"%s\" must contain two (time, ncol) or three dimensions (time, lat, lon)",
				strSourceData.c_str(), vecVariableStrings[v].c_str());
		}
		if (varData->num_dims() == 3) {
			if (varData->get_dim(1)->size() != dLatDeg.size()) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have latitude as its second dimension",
					strSourceData.c_str(), vecVariableStrings[v].c_str());
			}
			if (varData->get_dim(2)->size() != dLonDeg.size()) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have longitude as its third dimension",
					strSourceData.c_str(), vecVariableStrings[v].c_str());
			}
		}
		if (varData->num_dims() == 2) {
			if (varData->get_dim(1)->size() != dLatDeg.size()) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have a second dimension with the same length as \"lat\"",
					strSourceData.c_str(), vecVariableStrings[v].c_str());
			}
			if (dLonDeg.size() != dLatDeg.size()) {
				_EXCEPTION1("File \"%s\" dimension \"lat\" must be the same length as dimension \"lon\"",
					strSourceData.c_str());
			}
		}

		// Copy first dimension
		NcDim * dimVarFirst = varData->get_dim(0);
		_ASSERT(dimVarFirst != NULL);

		NcDim * dimVarFirstOut = ncOutput.get_dim(dimVarFirst->name());
		if (dimVarFirstOut == NULL) {
			dimVarFirstOut = ncOutput.add_dim(dimVarFirst->name(), dimVarFirst->size());
			if (dimVarFirstOut == NULL) {
				_EXCEPTION2("Error creating dimension \"%s\" in file \"%s\"",
					dimVarFirst->name(), strOutputData.c_str());
			}
			CopyNcVarIfExists(ncInput, ncOutput, dimVarFirst->name());
		}

		// Allocate data and count
		std::vector<float> dData;

		if (varData->num_dims() == 3) {
			dData.resize(dLatDeg.size() * dLonDeg.size());
		}
		if (varData->num_dims() == 2) {
			dData.resize(dLatDeg.size());
		}

		// For the first variable, create the map
		if (v == 0) {
			AnnounceStartBlock("Creating map");

			if (varData->num_dims() == 3) {
				sShpMap.resize(dLatDeg.size() * dLonDeg.size(), static_cast<size_t>(-1));
				for (size_t j = 0; j < dLatDeg.size(); j++) {
				for (size_t i = 0; i < dLonDeg.size(); i++) {
					size_t k = j * dLonDeg.size() + i;
		
					for (size_t s = 0; s < sShpRegionCount; s++) {
						double dStandardLonDeg = LonDegToStandardRange(dLonDeg[i]);
						if (vecLatLonBox[s].contains(dLatDeg[j], dStandardLonDeg)) {
							if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg[j], dStandardLonDeg)) {
								sShpMap[k] = s;
								break;
							}
						}
					}
				}
				}
			}

			if (varData->num_dims() == 2) {
				sShpMap.resize(dLatDeg.size(), static_cast<size_t>(-1));
				for (size_t k = 0; k < dLatDeg.size(); k++) {
					for (size_t s = 0; s < sShpRegionCount; s++) {
						double dStandardLonDeg = LonDegToStandardRange(dLonDeg[k]);
						if (vecLatLonBox[s].contains(dLatDeg[k], dStandardLonDeg)) {
							if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg[k], dStandardLonDeg)) {
								sShpMap[k] = s;
								break;
							}
						}
					}
				}
			}

			AnnounceEndBlock("Done");
		}

		// Create output variable
		NcVar * varDataOut = ncOutput.add_var(vecVariableStrings[v].c_str(), ncFloat, dimVarFirstOut, dimShp);
		if (varDataOut == NULL) {
			_EXCEPTION1("Error creating variable \"%s\" in output file", vecVariableStrings[v].c_str());
		}
		CopyNcVarAttributes(varData, varDataOut);

		// Fillvalue
		float dFillValue = FLT_MAX;
		NcAtt * attFillValue = varData->get_att("_FillValue");
		if (attFillValue != NULL) {
			dFillValue = attFillValue->as_float(0);
		}

		// Loop through all times
		for (size_t t = 0; t < dimVarFirst->size(); t++) {
	
			Announce("Time %lu/%lu", t, dimVarFirst->size());

			// Reset the counts
			std::fill(sDataCount.begin(), sDataCount.end(), 0);
	
			// Allocate output data
			std::vector<float> dDataOut(sShpRegionCount, 0.0f);
	
			// Load data
			if (varData->num_dims() == 3) {
				varData->set_cur(t,0,0);
				varData->get(&(dData[0]), 1, dLatDeg.size(), dLonDeg.size());
			} else {
				varData->set_cur(t,0);
				varData->get(&(dData[0]), 1, dLatDeg.size());
			}
	
			// Remap data
			for (size_t k = 0; k < sShpMap.size(); k++) {
				if ((!std::isnan(dData[k])) && (dData[k] != dFillValue)) {
					//printf("%1.5f %lu\n", dData[k], sShpMap[k]);
					if (sShpMap[k] != static_cast<size_t>(-1)) {
						dDataOut[sShpMap[k]] += static_cast<double>(dData[k]);
						sDataCount[sShpMap[k]]++;
					}
				}
			}
	
			// Normalize data
			for (size_t s = 0; s < sShpRegionCount; s++) {
				if (sDataCount[s] != 0) {
					dDataOut[s] /= static_cast<float>(sDataCount[s]);
				}
				//if (dDataOut[s] != 0.0) {
				//	printf("%1.5f\n", dDataOut[s]);
				//}
			}
	
			// Write to disk
			varDataOut->set_cur(t,0);
			varDataOut->put(&(dDataOut[0]), 1, sShpRegionCount);
		}

		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}

///////////////////////////////////////////////////////////////////////////////

