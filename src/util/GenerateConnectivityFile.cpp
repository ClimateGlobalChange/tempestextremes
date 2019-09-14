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

	// True if the input mesh contains concave elements
	bool fConcave;

	// Polynomial degree
	int nP;

	// Connectivity file type
	std::string strConnectType;

	// Output data file
	std::string strConnectFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshFile, "in_mesh", "");
		CommandLineBool(fConcave, "in_concave");
		CommandLineStringD(strConnectType, "out_type", "FV", "[FV|CGLL|DGLL]");
		CommandLineInt(nP, "out_np", 4);
		CommandLineString(strConnectFile, "out_connect", "");

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

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

