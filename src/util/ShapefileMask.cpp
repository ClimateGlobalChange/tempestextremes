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
#include "Variable.h"
#include "SimpleGrid.h"

#include "netcdfcpp.h"

#include <cfloat>
#include <algorithm>

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
	std::string strInputData;

	// Connectivity file
	std::string strConnectivity;

	// Grid file
	std::string strGridFile;

	// Data is regional
	bool fRegional;

	// Shapefile
	std::string strShapefile;

	// Lonlat bounds
	std::string strLonLatBounds;

	// Output data
	std::string strOutputData;

	// Output CSV
	std::string strOutputCSV;

	// Variable name
	std::string strVariables;

	// Output variable name
	std::string strOutputVariables;

	// Write the shape ids
	bool fWriteShapeIds;

	// Operation
	std::string strOperation;

	// Longitude name
	std::string strLongitudeName;

	// Latitude name
	std::string strLatitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strGridFile, "in_grid", "");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strShapefile, "shp", "");
		CommandLineStringD(strLonLatBounds, "lonlatbounds", "", "\"lon0,lat0,lon1,lat1\" (in degrees)");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputCSV, "out_csv", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputVariables, "varout", "");
		CommandLineBool(fWriteShapeIds, "writeshapeids");
		CommandLineStringD(strOperation, "op", "mask", "(mask|mean|sum|q25|median|q75)");
		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Check arguments
	if (strInputData.length() == 0) {
		_EXCEPTIONT("No input file (--in_data) specified");
	}
	if (strShapefile.length() == 0) {
		_EXCEPTIONT("No shapefile (--shp) specified");
	}
	if ((strOutputData.length() == 0) && (strOutputCSV.length() == 0)) {
		_EXCEPTIONT("No output file (--out_data) or (--out_csv) specified");
	}
	if (strVariables.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}
	if ((strOutputData.length() == 0) && (fWriteShapeIds)) {
		_EXCEPTIONT("Output file (--out_data) required when using (--writeshapeids)");
	}

	// Check operator
	STLStringHelper::ToLower(strOperation);
	if ((strOperation != "mask") &&
		(strOperation != "mean") &&
		(strOperation != "sum") &&
	    (strOperation != "q25") &&
	    (strOperation != "median") &&
	    (strOperation != "q75")
	) {
		_EXCEPTIONT("Invalid --op, expected \"mask\", \"mean\", \"q25\", \"median\" or \"q75\")");
	}

	// Open output text file
	FILE * fpCSV = NULL;
	if (strOutputCSV != "") {
		if (strOperation == "mask") {
			_EXCEPTIONT("--out_csv cannot be combined with --op \"mask\"");
		}
		fpCSV = fopen(strOutputCSV.c_str(), "w");
		if (fpCSV == NULL) {
			_EXCEPTION1("Unable to open output file \"%s\" for writing", strOutputCSV.c_str());
		}
	}

	// Load the shapefile
	AnnounceStartBlock("Loading shapefile");
	Mesh mesh;
	ReadShpFileAsMesh(strShapefile, mesh);
	Announce("Shapefile contains %lu polygons", mesh.faces.size());
	AnnounceEndBlock("Done");

	// Parse --lonlatbounds
	bool fHasLonLatBounds = false;
	LatLonBox<double> lonlatbounds;
	if (strLonLatBounds != "") {
		AnnounceStartBlock("Parsing lonlat bounds");

		std::vector<std::string> vecLonLatBoundsArg;
		STLStringHelper::ParseVariableList(strLonLatBounds, vecLonLatBoundsArg,",;");
		if (vecLonLatBoundsArg.size() != 4) {
			_EXCEPTIONT("Malformed --lonlatbounds: expected form \"lon0,lat0,lon1,lat1\"");
		}
		if (!STLStringHelper::IsFloat(vecLonLatBoundsArg[0]) ||
		    !STLStringHelper::IsFloat(vecLonLatBoundsArg[1]) ||
		    !STLStringHelper::IsFloat(vecLonLatBoundsArg[2]) ||
		    !STLStringHelper::IsFloat(vecLonLatBoundsArg[3])
		) {
			_EXCEPTIONT("Malformed --lonlatbounds: expected form \"lon0,lat0,lon1,lat1\"");
		}

		double dLonDeg0 = LonDegToStandardRange(std::stof(vecLonLatBoundsArg[0]));
		double dLatDeg0 = std::stof(vecLonLatBoundsArg[1]);
		double dLonDeg1 = LonDegToStandardRange(std::stof(vecLonLatBoundsArg[2]));
		double dLatDeg1 = std::stof(vecLonLatBoundsArg[3]);

		if ((fabs(dLatDeg0) > 90.0) || (fabs(dLatDeg1) > 90.0) || (dLatDeg1 < dLatDeg0)) {
			_EXCEPTION2("Malformed --lonlatbounds: Latitude must be in range [-90,90] with lat0 < lat1 (given: %1.5f %1.5f)",
				dLatDeg0, dLatDeg1);
		}

		lonlatbounds.set(dLonDeg0, dLonDeg1, dLatDeg0, dLatDeg1);
		fHasLonLatBounds = true;

		Announce("Limiting shapefiles to longitude [%g, %g] and latitude [%g, %g]", dLonDeg0, dLonDeg1, dLatDeg0, dLatDeg1);
		AnnounceEndBlock(NULL);
	}

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
	AnnounceStartBlock("Initializing");

	// Initialize input files
	NcFileVector vecFiles;
	vecFiles.ParseFromString(strInputData);
	if (vecFiles.size() != 1) {
		_EXCEPTIONT("Exactly one file can be processed by ShapefileMask at a time");
	}
	NcFile & ncInput = *(vecFiles[0]);

	// Define the SimpleGrid
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

		if (strGridFile == "") {
			grid.GenerateLatitudeLongitude(
				vecFiles[0],
				strLatitudeName,
				strLongitudeName,
				fRegional,
				false);
		} else {
			NcFile ncgridfile(strGridFile.c_str(), NcFile::ReadOnly);
			if (!ncgridfile.is_valid()) {
				_EXCEPTION1("Unable to open --in_grid file \"%s\"",
					strGridFile.c_str());
			}

			grid.GenerateLatitudeLongitude(
				&ncgridfile,
				strLatitudeName,
				strLongitudeName,
				fRegional,
				false);
		}

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Parse --var argument and build VariableRegistry
	std::vector<std::string> vecVariableNamesIn;
	std::vector<VariableIndex> vecVarIxIn;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varix;
			int iLast = varreg.FindOrRegisterSubStr(strVariablesTemp, &varix) + 1;

			vecVarIxIn.push_back(varix);
			vecVariableNamesIn.push_back( strVariablesTemp.substr(0,iLast-1) );
			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}
	}
	_ASSERT(vecVariableNamesIn.size() == vecVarIxIn.size());
 
	// Parse --varout argument
	std::vector<std::string> vecVariableNamesOut;
	if (strOutputVariables != "") {
		STLStringHelper::ParseVariableList(strOutputVariables, vecVariableNamesOut);
	} else {
		vecVariableNamesOut = vecVariableNamesIn;
	}
	if (vecVariableNamesOut.size() != vecVariableNamesIn.size()) {
		_EXCEPTION2("Inconsistent number of variables in --var (%lu) and --varout (%lu)",
			vecVarIxIn.size(), vecVariableNamesOut.size());
	}

	// Open output file
	NcFile * pncOutput = NULL;
	NcDim * dimShp = NULL;
	NcDim * dim0 = NULL;
	NcDim * dim1 = NULL;

	if (strOutputData != "") {
		Announce("Preparing output file");

		pncOutput = new NcFile(strOutputData.c_str(), NcFile::Replace);
		if (!pncOutput->is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
				strOutputData.c_str());
		}

		// Create shape index dimension in output file
		dimShp = pncOutput->add_dim("shpix", sShpRegionCount);
		if (dimShp == NULL) {
			_EXCEPTIONT("Error creating dimension \"shpix\" in output file");
		}

		// Create shape ids variable in output file
		if ((fWriteShapeIds) || (strOperation == "mask")) {
			if (grid.m_nGridDim.size() == 1) {
				dim0 = pncOutput->add_dim("ncol", grid.m_nGridDim[0]);

			} else if (grid.m_nGridDim.size() == 2) {
				dim0 = pncOutput->add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
				dim1 = pncOutput->add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);

				CopyNcVarIfExists(
					ncInput,
					*pncOutput,
					strLatitudeName.c_str(),
					true,
					true);

				CopyNcVarIfExists(
					ncInput,
					*pncOutput,
					strLongitudeName.c_str(),
					true,
					true);

			} else {
				_EXCEPTIONT("Only 1D or 2D spatial data supported");
			}
		}
	}

	AnnounceEndBlock("Done");

	// Get time information
	bool fFileHasTime = (NcGetTimeVariable(ncInput) != NULL);
	NcTimeDimension vecTimes;
	NcDim * dimTime = NULL;
	if (fFileHasTime) {
		ReadCFTimeDataFromNcFile(
			vecFiles[0],
			strInputData,
			vecTimes,
			false);

		if (pncOutput != NULL) {
			WriteCFTimeDataToNcFile(
				pncOutput,
				strOutputData,
				vecTimes);

			dimTime = NcGetTimeDimension(*pncOutput);
			_ASSERT(dimTime != NULL);
		}

	} else {
		vecTimes.push_back(Time(Time::CalendarNone));
	}

	// Build the map from the input data to the shapefiles
	AnnounceStartBlock("Building map");

	std::vector<int> nShpMap;

	nShpMap.resize(grid.m_dLat.GetRows(), static_cast<int>(-1));
	for (size_t k = 0; k < grid.m_dLat.GetRows(); k++) {
		double dLatDeg = RadToDeg(grid.m_dLat[k]);
		double dStandardLonDeg = LonDegToStandardRange(RadToDeg(grid.m_dLon[k]));
		if (fHasLonLatBounds) {
			if (!lonlatbounds.contains(dLatDeg, dStandardLonDeg)) {
				continue;
			}
		}
		for (size_t s = 0; s < sShpRegionCount; s++) {
			if (vecLatLonBox[s].contains(dLatDeg, dStandardLonDeg)) {
				if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg, dStandardLonDeg)) {
					nShpMap[k] = static_cast<int>(s);
					break;
				}
			}
		}
	}

	AnnounceEndBlock("Done");

	// Write the shape ids to a file
	if (fWriteShapeIds) {
		_ASSERT(pncOutput != NULL);

		AnnounceStartBlock("Writing shape ids to file");

		NcVar * varShapeIds = NULL;
		if (grid.m_nGridDim.size() == 1) {
			varShapeIds = pncOutput->add_var("shapeid", ncInt, dim0);
			varShapeIds->add_att("_FillValue", (int)(-1));
			if (varShapeIds == NULL) {
				_EXCEPTION1("Unable to create variable \"shapeid\" in NetCDF file \"%s\"", strOutputData.c_str());
			}
			varShapeIds->set_cur((long)0);
			varShapeIds->put(&(nShpMap[0]), dim0->size());

		} else {
			varShapeIds = pncOutput->add_var("shapeid", ncInt, dim0, dim1);
			varShapeIds->add_att("_FillValue", (int)(-1));
			if (varShapeIds == NULL) {
				_EXCEPTION1("Unable to create variable \"shapeid\" in NetCDF file \"%s\"", strOutputData.c_str());
			}
			varShapeIds->set_cur(0, 0);
			varShapeIds->put(&(nShpMap[0]), dim0->size(), dim1->size());
		}

		AnnounceEndBlock("Done");
	}

	// Build up the processing queue and define output variables
	varreg.ClearProcessingQueue();

	std::vector<NcVar *> vecNcVarOut(vecVarIxIn.size());
	for (int v = 0; v < vecVarIxIn.size(); v++) {
		NcVar * varOut = NULL;

		DimInfoVector vecAuxDimInfo;
		varreg.AppendVariableToProcessingQueue(
			vecFiles,
			grid,
			vecVarIxIn[v],
			&vecAuxDimInfo);

		// Create output variable
		if (pncOutput != NULL) {
			std::vector<NcDim *> vecDimOut;

			if (fFileHasTime) {
				vecDimOut.push_back(dimTime);
			}

			for (long d = 0; d < vecAuxDimInfo.size(); d++) {
				NcDim * dim =
					AddNcDimOrUseExisting(
						*pncOutput,
						vecAuxDimInfo[d].name,
						vecAuxDimInfo[d].size);

				vecDimOut.push_back(dim);
			}

			vecDimOut.push_back(dimShp);
			if (strOperation == "mask") {
				vecDimOut.push_back(dim0);
				if (grid.m_nGridDim.size() == 2) {
					vecDimOut.push_back(dim1);
				}
			}

			vecNcVarOut[v] = pncOutput->add_var(
				vecVariableNamesOut[v].c_str(),
				ncFloat,
				vecDimOut.size(),
				const_cast<const NcDim**>(&(vecDimOut[0])));

			NcVar * ncvarIn = ncInput.get_var(vecVariableNamesIn[v].c_str());
			if (ncvarIn != NULL) {
				CopyNcVarAttributes(ncvarIn, vecNcVarOut[v]);
			}

			if (vecNcVarOut[v] == NULL) {
				_EXCEPTION2("Error creating variable \"%s\" in file \"%s\"",
					vecVariableNamesOut[v].c_str(),
					strOutputData.c_str());
			}
		}
	}

	// Loop through all times in this file
	for (int t = 0; t < vecTimes.size(); t++) {
		if (fFileHasTime) {
			AnnounceStartBlock("Time %s (%lu/%lu)",
				vecTimes[t].ToString().c_str(),
				t+1, vecTimes.size());
		}

		vecFiles.SetTime(vecTimes[t]);

		// Loop through all variables
		varreg.ResetProcessingQueue();
		while(varreg.AdvanceProcessingQueue()) {
			size_t v = varreg.GetProcessingQueueVarPos();

			_ASSERT(v < vecVariableNamesIn.size());
			_ASSERT(v < vecVariableNamesOut.size());

			// Load the data for the search variable
			Variable & var = varreg.GetProcessingQueueVariable();
			Announce("Processing \"%s\" -> \"%s\"",
				var.ToString(varreg).c_str(),
				vecVariableNamesOut[v].c_str());

			var.LoadGridData(varreg, vecFiles, grid);
			const DataArray1D<float> & dataState = var.GetData();

			// Check _FillValue
			float dFillValue = var.GetFillValueFloat();

			// Get the output position and size
			VariableAuxIndex lPos = varreg.GetProcessingQueueAuxIx();
			if (fFileHasTime) {
				lPos.insert(lPos.begin(), t);
			}
			VariableAuxIndex lSize(lPos.size(), 1);

			// Mask the data
			if (strOperation == "mask") {

				_ASSERT(pncOutput != NULL);

				size_t sShpIxPos = lPos.size();
				lPos.push_back(0);
				lSize.push_back(1);

				lPos.push_back(0);
				lSize.push_back(dim0->size());
				if (grid.m_nGridDim.size() == 2) {
					lPos.push_back(0);
					lSize.push_back(dim1->size());
				}

				for (int s = 0; s < static_cast<int>(sShpRegionCount); s++) {
					std::vector<float> dataOut(dataState.GetRows(), dFillValue);
					for (size_t k = 0; k < nShpMap.size(); k++) {
						if (nShpMap[k] == s) {
							dataOut[k] = dataState[k];
						}
					}

					lPos[sShpIxPos] = static_cast<long>(s);
					vecNcVarOut[v]->set_cur(&(lPos[0]));
					vecNcVarOut[v]->put(&(dataOut[0]), &(lSize[0]));
				}

			// Perform reduce operation on the data
			} else {

				lPos.push_back(0);
				lSize.push_back(sShpRegionCount);
		
				// Allocate output data
				std::vector<float> dataOut(sShpRegionCount, dFillValue);

				if (fpCSV != NULL) {
					if (fFileHasTime) {
						fprintf(fpCSV, "%s, \"%s\"", vecVariableNamesOut[v].c_str(), vecTimes[t].ToString().c_str());
					} else {
						fprintf(fpCSV, "%s", vecVariableNamesOut[v].c_str());
					}
				}

				// Calculate the mean within all shapes
				if ((strOperation == "mean") || (strOperation == "sum")) {
					_ASSERT(grid.m_dArea.GetRows() == nShpMap.size());
					//std::vector<size_t> sDataCount(sShpRegionCount, 0);
					std::vector<double> dAccumulatedArea(sShpRegionCount, 0.0);
					std::fill(dataOut.begin(), dataOut.end(), 0.0f);
					for (size_t k = 0; k < nShpMap.size(); k++) {
						if ((!std::isnan(dataState[k])) && (dataState[k] != dFillValue)) {
							if (nShpMap[k] != static_cast<size_t>(-1)) {
								dataOut[nShpMap[k]] += static_cast<double>(dataState[k]) * grid.m_dArea[k];
								dAccumulatedArea[nShpMap[k]] += grid.m_dArea[k];
								//sDataCount[nShpMap[k]]++;
							}
						}
					}
					bool fNormalizeByArea = (strOperation == "mean");
					for (size_t s = 0; s < sShpRegionCount; s++) {
						if (fNormalizeByArea && (dAccumulatedArea[s] != 0)) {
							dataOut[s] /= static_cast<float>(dAccumulatedArea[s]);
						}
						if (fpCSV != NULL) {
							fprintf(fpCSV, ", %g", dataOut[s]);
						}
					}
	
				// Generate arrays of data in each shapefile
				} else {
					std::vector< std::vector<float> > vecShpData;
					vecShpData.resize(sShpRegionCount);
	
					for (size_t k = 0; k < nShpMap.size(); k++) {
						if ((!std::isnan(dataState[k])) && (dataState[k] != dFillValue)) {
							if (nShpMap[k] != (-1)) {
								vecShpData[nShpMap[k]].push_back(dataState[k]);
							}
						}
					}
	
					for (size_t s = 0; s < vecShpData.size(); s++) {
						if (vecShpData[s].size() == 0) {
							dataOut[s] = 0.0;
							continue;
						}
						if (vecShpData[s].size() == 1) {
							dataOut[s] = vecShpData[s][0];
							continue;
						}
						if (vecShpData[s].size() == 2) {
							dataOut[s] = 0.5 * (vecShpData[s][0] + vecShpData[s][1]);
							continue;
						}
	
						std::sort(vecShpData[s].begin(), vecShpData[s].end());

						size_t ix0 = 0;
						size_t ix1 = 0;
						if ((strOperation == "q25") || (fpCSV != NULL)) {
							ix1 = (vecShpData[s].size()+1) / 4;
	
							if ((vecShpData[s].size()+1) % 4 == 0) {
								ix0 = ix1;
							} else {
								ix0 = ix1 - 1;
							}
	
							dataOut[s] = 0.5 * (vecShpData[s][ix0] + vecShpData[s][ix1]);
	
							if (fpCSV != NULL) {
								fprintf(fpCSV, ", %g", dataOut[s]);
							}
						}
	
						if ((strOperation == "median") || (fpCSV != NULL)) {
							ix1 = (vecShpData[s].size()+1) / 2;
							if ((vecShpData[s].size()+1) % 2 == 0) {
								ix0 = ix1;
							} else {
								ix0 = ix1 - 1;
							}
	
							dataOut[s] = 0.5 * (vecShpData[s][ix0] + vecShpData[s][ix1]);
	
							if (fpCSV != NULL) {
								fprintf(fpCSV, ", %g", dataOut[s]);
							}
						}
	
						if ((strOperation == "q75") || (fpCSV != NULL)) {
							ix1 = (vecShpData[s].size()+1) * 3 / 4;
							if ((vecShpData[s].size()+1) % 4 == 0) {
								ix0 = ix1;
							} else {
								ix0 = ix1 - 1;
							}
	
							dataOut[s] = 0.5 * (vecShpData[s][ix0] + vecShpData[s][ix1]);
	
							if (fpCSV != NULL) {
								fprintf(fpCSV, ", %g", dataOut[s]);
							}
						}
					}
				}

				if (pncOutput != NULL) {
					vecNcVarOut[v]->set_cur(&(lPos[0]));
					vecNcVarOut[v]->put(&(dataOut[0]), &(lSize[0]));
				}

				if (fpCSV != NULL) {
					fprintf(fpCSV, "\n");
				}
			}

			NcError err;
			if (err.get_err() != NC_NOERR) {
				_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
			}
		}

		if (fFileHasTime) {
			AnnounceEndBlock(NULL);
		}
	}

	AnnounceBanner();

	// Cleanup
	if (pncOutput != NULL) {
		delete pncOutput;
	}

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}

///////////////////////////////////////////////////////////////////////////////

