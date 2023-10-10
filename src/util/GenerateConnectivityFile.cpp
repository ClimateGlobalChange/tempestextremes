///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateConnectivityFile.cpp
///	\author  Paul Ullrich
///	\version January 28, 2019
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
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
#include "STLStringHelper.h"
#include "GridElements.h"
#include "SimpleGrid.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input mesh file
	std::string strMeshFile;

	// Input data file
	std::string strDataFile;

	// True if the input mesh contains concave elements
	bool fConcave;

	// Polynomial degree
	int nP;

	// Connectivity file type
	std::string strConnectType;

	// Output data file
	std::string strConnectFile;

	// Latitude name
	std::string strLatitudeName;

	// Longitude name
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshFile, "in_mesh", "");
		CommandLineString(strDataFile, "in_data", "");
		CommandLineBool(fConcave, "in_concave");
		CommandLineStringD(strConnectType, "out_type", "FV", "[FV|CGLL|DGLL]");
		CommandLineInt(nP, "out_np", 4);
		CommandLineString(strConnectFile, "out_connect", "");
		CommandLineString(strLatitudeName, "latname", "lat");
		CommandLineString(strLongitudeName, "lonname", "lon");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Convert connectivity type to lowercase
	STLStringHelper::ToLower(strConnectType);

	// Check arguments
	if ((strConnectType != "fv") &&
	    (strConnectType != "cgll") &&
	    (strConnectType != "dgll")
	) {
		_EXCEPTIONT("Invalid --out_type.  Expected \"FV|CGLL|DGLL\"");
	}

	if ((nP < 2) && (strConnectType == "cgll")) {
		_EXCEPTIONT("Invalid --out_np for --out_type CGLL.  Expected np > 1.");
	}

	if ((nP < 1) && (strConnectType == "dgll")) {
		_EXCEPTIONT("Invalid --out_np for --out_type DGLL.  Expected np > 0.");
	}

	// Input mesh file
	if (strMeshFile != "") {

		// Load in the mesh
		AnnounceStartBlock("Loading mesh");
		Mesh mesh(strMeshFile);
		AnnounceEndBlock("Done");

		// Calculate connectivity information
		AnnounceStartBlock("Calculating connectivity information");
		mesh.RemoveZeroEdges();
		mesh.CalculateFaceAreas(fConcave);
		mesh.ConstructEdgeMap();
		AnnounceEndBlock("Done");

		// Generate SimpleGrid
		AnnounceStartBlock("Converting mesh to connectivity format");
		SimpleGrid grid;
		if (strConnectType == "fv") {
			grid.FromMeshFV(mesh);
		} else {
			grid.FromMeshFE(mesh, (strConnectType == "cgll"), nP);
		}
		AnnounceEndBlock("Done");

		// Writing data to file
		AnnounceStartBlock("Writing connectivity file");
		grid.ToFile(strConnectFile);
		AnnounceEndBlock("Done");
	}

	// Input data file
	if (strDataFile != "") {

		NcFile ncfilein(strDataFile.c_str());
		if (!ncfilein.is_valid()) {
			_EXCEPTION1("Unable to open data file \"%s\"", strDataFile.c_str());
		}

		NcVar * varLat = ncfilein.get_var(strLatitudeName.c_str());
		if (varLat == NULL) {
			_EXCEPTION2("Data file \"%s\" does not contain variable \"%s\"", strDataFile.c_str(), strLatitudeName.c_str());
		}

		NcVar * varLon = ncfilein.get_var(strLongitudeName.c_str());
		if (varLon == NULL) {
			_EXCEPTION2("Data file \"%s\" does not contain variable \"%s\"", strDataFile.c_str(), strLongitudeName.c_str());
		}

		if (varLat->get_dim(0)->size() != varLon->get_dim(0)->size()) {
			_EXCEPTIONT("Data file latitude and longitude must have same length");
		}

		// Load data
		std::vector<double> dLat(varLat->get_dim(0)->size());
		varLat->get(&(dLat[0]), dLat.size());

		std::vector<double> dLon(varLon->get_dim(0)->size());
		varLon->get(&(dLon[0]), dLon.size());

		bool fRadians = false;

		NcAtt * attUnits = varLat->get_att("units");
		if (attUnits != NULL) {
			std::string strUnits = attUnits->as_string(0);
			if ((strUnits == "radians") || (strUnits == "radian")) {
				fRadians = true;
			}
		}

		if (fRadians) {
			for (size_t s = 0; s < dLon.size(); s++) {
				dLon[s] = RadToDeg(dLon[s]);
				dLat[s] = RadToDeg(dLat[s]);
			}
		}

		// Generate the grid
		SimpleGrid grid;
		grid.m_nGridDim.resize(1);
		grid.m_nGridDim[0] = dLon.size();

		grid.m_dLon.Allocate(dLon.size());
		memcpy(&(grid.m_dLon[0]), &(dLon[0]), dLon.size() * sizeof(double));

		grid.m_dLat.Allocate(dLon.size());
		memcpy(&(grid.m_dLat[0]), &(dLat[0]), dLat.size() * sizeof(double));

		grid.m_dArea.Allocate(dLon.size());
		for (size_t s = 0; s < dLon.size(); s++) {
			grid.m_dArea[s] = 1.0;
		}

		grid.m_vecConnectivity.resize(dLon.size());

		// Writing data to file
		AnnounceStartBlock("Writing connectivity file");
		grid.ToFile(strConnectFile);
		AnnounceEndBlock("Done");

	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

