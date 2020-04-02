///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileCompose.cpp
///	\author  Paul Ullrich
///	\version March 13, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
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
#include "Variable.h"
#include "AutoCurator.h"
#include "DataArray2D.h"
#include "ArgumentTree.h"
#include "STLStringHelper.h"
#include "NodeFileUtilities.h"
#include "NetCDFUtilities.h"
#include "ClosedContourOp.h"
#include "SimpleGridUtilities.h"
#include "GridElements.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>
#include <cmath>
/*
#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif
*/

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
/*
#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Not yet implemented
	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
	if (nMPISize > 1) {
		std::cout << "Sorry!  Parallel processing with MPI is not yet implemented" << std::endl;
		return (-1);
	}
#endif
*/
	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);
/*
	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();
*/
try {

	// Input nodefile
	std::string strInputNodeFile;

	// List of input nodefiles
	std::string strInputNodeFileList;

	// Input nodefile type
	std::string strInputNodeFileType;

	// Input format (columns of in_file)
	std::string strInputFormat;

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Connectivity file
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Output grid type
	std::string strOutputGrid;

	// Output data file
	std::string strOutputData;

	// Start of time range to use
	std::string strTimeBegin;

	// End of time range to use
	std::string strTimeEnd;

	// List of variables to output
	std::string strVariables;

	// Grid spacing of output (Cartesian great-circle or radial distance)
	double dDeltaX;

	// Resolution of output (Cartesian grid size or number of radial points)
	int nResolutionX;

	// Resolution of the output (azimuthal points)
	int nResolutionA;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_file", "");
		//CommandLineString(strInputNodeFileList, "in_file_list", "");
		CommandLineStringD(strInputNodeFileType, "in_file_type", "SN", "[DCU|SN]");
		CommandLineString(strInputFormat, "in_fmt", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");

		CommandLineStringD(strOutputGrid, "out_grid", "XY", "[XY|RAD]");
		CommandLineString(strOutputData, "out_data", "");

		//CommandLineString(strTimeBegin, "time_begin", "");
		//CommandLineString(strTimeEnd, "time_end", "");

		CommandLineString(strVariables, "var", "");

		CommandLineDouble(dDeltaX, "dx", 0.5);
		CommandLineInt(nResolutionX, "resx", 11);
		CommandLineInt(nResolutionA, "resa", 16);

		CommandLineString(strLatitudeName, "latname", "lat");
		CommandLineString(strLongitudeName, "lonname", "lon");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registries
	VariableRegistry varregIn;
	VariableRegistry varregOut;

	// Create autocurator
	AutoCurator autocurator;

	// Check arguments
	if (strInputNodeFile.length() == 0) {
		_EXCEPTIONT("No input file (--in_nodefile) specified");
	}
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if (strVariables.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}

	// Output grid options
	STLStringHelper::ToLower(strOutputGrid);
	if ((strOutputGrid != "xy") && (strOutputGrid != "rad")) {
		_EXCEPTIONT("Grid type of output (--out_grid) must be \"xy\" or \"rad\"");
	}
	if (dDeltaX <= 0.0) {
		_EXCEPTIONT("Grid spacing of output (--dx) must be nonnegative");
	}
	if (nResolutionX < 1) {
		_EXCEPTIONT("Resolution of output (--resx) must be nonnegative");
	}
	if ((strOutputGrid == "rad") && (nResolutionA < 8)) {
		_EXCEPTIONT("Resolution of output (--resa) must be >= 8");
	}

	// Input file type
	NodeFile::PathType iftype;
	if (strInputNodeFileType == "DCU") {
		iftype = NodeFile::PathTypeDCU;
	} else if (strInputNodeFileType == "SN") {
		iftype = NodeFile::PathTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_nodefile_type, expected \"SN\" or \"DCU\"");
	}

	// NodeFile
	NodeFile nodefile;

	// Parse --in_fmt string
	ColumnDataHeader cdhInput;
	cdhInput.Parse(strInputFormat);

	// Parse --var argument
	VariableIndexVector vecVarIxIn;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varixIn;
			int iLast = varregIn.FindOrRegisterSubStr(strVariablesTemp, &varixIn) + 1;

			vecVarIxIn.push_back(varixIn);

			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
 
	// Generate a list of dependent base variables for each variable
	std::vector< std::vector<std::string> > vecvecDependentVarNames;
	vecvecDependentVarNames.resize(vecVarIxIn.size());
	for (int v = 0; v < vecVarIxIn.size(); v++) {
		varregIn.GetDependentVariableNames(
			vecVarIxIn[v],
			vecvecDependentVarNames[v]);

		if (vecvecDependentVarNames[v].size() == 0) {
			_EXCEPTION1("Variable \"%s\" has no dependent base variables",
				varregIn.GetVariableString(vecVarIxIn[v]).c_str());
		}
	}

	// Curate input data
	if (strInputData.length() != 0) {
		AnnounceStartBlock("Autocurating in_data");
		autocurator.IndexFiles(strInputData);

	} else {
		AnnounceStartBlock("Autocurating in_data_list");
		std::ifstream ifInputDataList(strInputDataList.c_str());
		if (!ifInputDataList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputDataList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			Announce(strFileLine.c_str());
			autocurator.IndexFiles(strFileLine);
		}
	}
	AnnounceEndBlock("Done");

	// Define the SimpleGrid for the input
	SimpleGrid grid;

	// Check for connectivity file
	if (strConnectivity != "") {
		AnnounceStartBlock("Generating grid information from connectivity file");
		grid.FromFile(strConnectivity);
		AnnounceEndBlock("Done");

	// No connectivity file; check for latitude/longitude dimension
	} else {
		AnnounceStartBlock("No connectivity file specified");
		Announce("Attempting to generate latitude-longitude grid from data file");
		const std::vector<std::string> & vecFiles = autocurator.GetFilenames();

		if (vecFiles.size() < 1) {
			_EXCEPTIONT("No data files specified; unable to generate grid");
		}

		NcFile ncFile(vecFiles[0].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", vecFiles[0].c_str());
		}

		grid.GenerateLatitudeLongitude(
			&ncFile,
			fRegional,
			strLatitudeName,
			strLongitudeName);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Build the KD tree for the grid
	AnnounceStartBlock("Generating KD tree on grid");
	grid.BuildKDTree();
	AnnounceEndBlock("Done");

	// Load input file list
	std::vector<std::string> vecInputNodeFiles;

	if (strInputNodeFile.length() != 0) {
		vecInputNodeFiles.push_back(strInputNodeFile);

	} else {
		std::ifstream ifInputFileList(strInputNodeFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputNodeFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputNodeFiles.push_back(strFileLine);
		}
	}

	// Open output file
	AnnounceStartBlock("Preparing output file");
	NcFile ncoutfile(strOutputData.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputData.c_str());
	}

	// Write dimensions
	NcDim * dimX;
	NcDim * dimY;

	if (strOutputGrid == "xy") {
		dimX = ncoutfile.add_dim("x", nResolutionX);
		dimY = ncoutfile.add_dim("y", nResolutionX);

		DataArray1D<double> dX(nResolutionX);
		for (int i = 0; i < nResolutionX; i++) {
			dX[i] = dDeltaX * (
				static_cast<double>(i)
				- 0.5 * static_cast<double>(nResolutionX-1));
		}

		NcVar * varX = ncoutfile.add_var("x", ncDouble, dimX);
		varX->put(&(dX[0]), nResolutionX);

		NcVar * varY = ncoutfile.add_var("y", ncDouble, dimY);
		varY->put(&(dX[0]), nResolutionX);

	} else if (strOutputGrid == "rad") {
		dimX = ncoutfile.add_dim("az", nResolutionA);
		dimY = ncoutfile.add_dim("r", nResolutionX);

		DataArray1D<double> dAz(nResolutionA);
		for (int i = 0; i < nResolutionA; i++) {
			dAz[i] = 360.0 * static_cast<double>(i)
				/ static_cast<double>(nResolutionA);
		}

		DataArray1D<double> dR(nResolutionX);
		for (int i = 0; i < nResolutionA; i++) {
			dR[i] = dDeltaX * (static_cast<double>(i) + 0.5);
		}

		NcVar * varAz = ncoutfile.add_var("az", ncDouble, dimX);
		varAz->put(&(dAz[0]), nResolutionA);

		NcVar * varR = ncoutfile.add_var("r", ncDouble, dimY);
		varR->put(&(dR[0]), nResolutionX);

	} else {
		_EXCEPTIONT("Invalid grid");
	}
	AnnounceEndBlock("Done");

	// Vector of output data
	int nDataInstances = 0;

	std::vector< DataArray1D<float> > vecOutputData;
	vecOutputData.resize(vecVarIxIn.size());

	// Vector of output variables
	std::vector<NcVar*> vecOutputVar;

	// A map from Times to file lines, used for StitchNodes formatted output
	TimeToPathNodeMap & mapTimeToPathNode = nodefile.GetTimeToPathNodeMap();

	// Vector of Path information
	PathVector & pathvec = nodefile.GetPathVector();

	// Flag indicating output is initialized
	bool fOutputInitialized = false;

	// Loop over all nodefiles
	for (int f = 0; f < vecInputNodeFiles.size(); f++) {

		AnnounceStartBlock("Processing input (%s)", vecInputNodeFiles[f].c_str());

		// Read contents of NodeFile into PathVector
		AnnounceStartBlock("Reading file");
		nodefile.Read(
			vecInputNodeFiles[f],
			iftype,
			cdhInput,
			grid,
			autocurator.GetCalendarType());
		AnnounceEndBlock("Done");

		// Loop through all Times in the NodeFile
		for (auto iter = mapTimeToPathNode.begin(); iter != mapTimeToPathNode.end(); iter++) {

			// Generate a NcFileVector at this Time
			AnnounceStartBlock("Time %s", iter->first.ToString().c_str());

			NcFileVector vecncDataFiles;
			int iTime;

			autocurator.FindFilesAtTime(
				iter->first,
				vecncDataFiles,
				iTime);

			if (vecncDataFiles.size() == 0) {
				_EXCEPTION1("Time (%s) does not exist in input data fileset",
					iter->first.ToString().c_str());
			}

			// If this is the first time through the loop load the aux dimension info
			// and generate the output data structures and file.
			if (!fOutputInitialized) {

				AnnounceStartBlock("Initializing output variables");

				// Loop through all variables
				for (int v = 0; v < vecVarIxIn.size(); v++) {

					// Variable name
					std::string strVarName =
						varregIn.GetVariableString(vecVarIxIn[v]);

					// Get auxiliary dimension info and verify consistency
					DimInfoVector vecAuxDimInfo;

					varregIn.GetAuxiliaryDimInfoAndVerifyConsistency(
						vecncDataFiles,
						grid,
						vecvecDependentVarNames[v],
						vecAuxDimInfo);

					// Generate output variables
					if (strOutputGrid == "xy") {
						vecAuxDimInfo.push_back(DimInfo("y", nResolutionX));
						vecAuxDimInfo.push_back(DimInfo("x", nResolutionX));

						vecOutputData[v].Allocate(nResolutionX * nResolutionX);

					} else if (strOutputGrid == "rad") {
						vecAuxDimInfo.push_back(DimInfo("r", nResolutionX));
						vecAuxDimInfo.push_back(DimInfo("az", nResolutionA));

						vecOutputData[v].Allocate(nResolutionA * nResolutionX);
					}

					std::vector<NcDim *> vecNcDim;
					vecNcDim.resize(vecAuxDimInfo.size());
					for (int d = 0; d < vecAuxDimInfo.size(); d++) {
						vecNcDim[d] = ncoutfile.get_dim(vecAuxDimInfo[d].name.c_str());
						if (vecNcDim[d] == NULL) {
							vecNcDim[d] = ncoutfile.add_dim(
								vecAuxDimInfo[d].name.c_str(),
								vecAuxDimInfo[d].size);
						} else {
							if (vecNcDim[d]->size() != vecAuxDimInfo[d].size) {
								_EXCEPTION4("Dimension size mismatch when initializing variable \"%s\": Expected dimension \"%s\" to have size \"%li\" (found \"%li\")",
									strVarName.c_str(),
									vecAuxDimInfo[d].name.c_str(),
									vecAuxDimInfo[d].size,
									vecNcDim[d]->size());
							}
						}
					}

					// Generate variable
					NcVar * pvar =
						ncoutfile.add_var(
							strVarName.c_str(),
							ncFloat,
							vecNcDim.size(),
							const_cast<const NcDim**>(&(vecNcDim[0])));

					vecOutputVar.push_back(pvar);
				}

				// Done
				fOutputInitialized = true;

				AnnounceEndBlock("Done");
			}

			// Generate the SimpleGrid for each node
			AnnounceStartBlock("Building composites");
			const PathNodeIndexVector & vecPathNodes = iter->second;

			// Loop through all Variables
			for (int v = 0; v < vecVarIxIn.size(); v++) {

				// Load the data for the search variable
				Variable & var = varregIn.Get(vecVarIxIn[v]);
				var.LoadGridData(varregIn, vecncDataFiles, grid, iTime);
				const DataArray1D<float> & dataState = var.GetData();
				_ASSERT(dataState.GetRows() == grid.GetSize());

				// Loop through all PathNodes
				for (int p = 0; p < vecPathNodes.size(); p++) {
					nDataInstances++;

					const Path & path = pathvec[vecPathNodes[p].first];
					const PathNode & pathnode = path[vecPathNodes[p].second];
					int ixOrigin = static_cast<int>(pathnode.m_gridix);

					_ASSERT((ixOrigin >= 0) && (ixOrigin < grid.GetSize()));

					// Generate the SimpleGrid for this pathnode
					SimpleGrid gridNode;
					if (strOutputGrid == "xy") {
						gridNode.GenerateRectilinearStereographic(
							grid.m_dLon[ixOrigin],
							grid.m_dLat[ixOrigin],
							nResolutionX,
							dDeltaX);

					} else if (strOutputGrid == "rad") {
						gridNode.GenerateRadialStereographic(
							grid.m_dLon[ixOrigin],
							grid.m_dLat[ixOrigin],
							nResolutionX,
							nResolutionA,
							dDeltaX);
					}

					// Sample all points on the stereographic grid
					for (int i = 0; i < gridNode.GetSize(); i++) {

						//printf("%1.5f %1.5f %1.5f %1.5f\n",
						//	grid.m_dLon[ixOrigin],
						//	grid.m_dLat[ixOrigin],
						//	gridNode.m_dLon[i],
						//	gridNode.m_dLat[i]);

						int ixGridIn =
							grid.NearestNode(
								gridNode.m_dLon[i],
								gridNode.m_dLat[i]);

						vecOutputData[v][i] +=
							dataState[ixGridIn];
					}
				}
			}

			AnnounceEndBlock("Done");
			AnnounceEndBlock("Done");
		}

		// Average all Variables
		if (nDataInstances != 0) {
			for (int v = 0; v < vecVarIxIn.size(); v++) {
				for (int i = 0; i < vecOutputData[v].GetRows(); i++) {
					vecOutputData[v][i] /= static_cast<float>(nDataInstances);
				}
			}
		}

		// Write output
		AnnounceStartBlock("Writing output");
		_ASSERT(vecVarIxIn.size() == vecOutputData.size());
		if (strOutputGrid == "xy") {
			for (int v = 0; v < vecVarIxIn.size(); v++) {
				vecOutputVar[v]->put(
					&(vecOutputData[v][0]),
					nResolutionX,
					nResolutionX);
			}

		} else if (strOutputGrid == "rad") {
			for (int v = 0; v < vecVarIxIn.size(); v++) {
				vecOutputVar[v]->put(
					&(vecOutputData[v][0]),
					nResolutionX,
					nResolutionA);
			}
		}
		AnnounceEndBlock("Done");

		// Done processing this nodefile
		AnnounceEndBlock("Done");
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
/*
#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
*/
}

///////////////////////////////////////////////////////////////////////////////


