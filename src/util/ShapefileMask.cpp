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

	// Data is regional
	bool fRegional;

	// Shapefile
	std::string strShapefile;

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
		CommandLineBool(fRegional, "regional");

		CommandLineString(strShapefile, "shp", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputCSV, "out_csv", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputVariables, "varout", "");
		CommandLineBool(fWriteShapeIds, "writeshapeids");
		CommandLineStringD(strOperation, "op", "mask", "(mask|mean|q25|median|q75)");
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
	if (strOutputData.length() == 0) {
		_EXCEPTIONT("No output file (--out_data) specified");
	}
	if (strVariables.length() == 0) {
		_EXCEPTIONT("No variables (--var) specified");
	}

	// Check operator
	STLStringHelper::ToLower(strOperation);
	if ((strOperation != "mask") &&
		(strOperation != "mean") &&
	    (strOperation != "q25") &&
	    (strOperation != "median") &&
	    (strOperation != "q75")
	) {
		_EXCEPTIONT("Invalid --op, expected \"mask\", \"mean\", \"q25\", \"median\" or \"q75\")");
	}
/*
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
*/
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

		grid.GenerateLatitudeLongitude(
			vecFiles[0],
			strLatitudeName,
			strLongitudeName,
			fRegional,
			false);

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
	NcFile ncOutput(strOutputData.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputData.c_str());
	}

	// Create shape index dimension in output file
	NcDim * dimShp = ncOutput.add_dim("shpix", sShpRegionCount);
	if (dimShp == NULL) {
		_EXCEPTIONT("Error creating dimension \"shpix\" in output file");
	}

	// Create shape ids variable in output file
	NcDim * dim0 = NULL;
	NcDim * dim1 = NULL;
	if ((fWriteShapeIds) || (strOperation == "mask")) {
		if (grid.m_nGridDim.size() == 1) {
			dim0 = ncOutput.add_dim("ncol", grid.m_nGridDim[0]);

		} else if (grid.m_nGridDim.size() == 2) {
			dim0 = ncOutput.add_dim(strLatitudeName.c_str(), grid.m_nGridDim[0]);
			dim1 = ncOutput.add_dim(strLongitudeName.c_str(), grid.m_nGridDim[1]);

			CopyNcVarIfExists(
				ncInput,
				ncOutput,
				strLatitudeName.c_str(),
				true,
				true);

			CopyNcVarIfExists(
				ncInput,
				ncOutput,
				strLongitudeName.c_str(),
				true,
				true);

		} else {
			_EXCEPTIONT("Only 1D or 2D spatial data supported");
		}
	}

	AnnounceEndBlock("Done");

	// Get time information
	bool fFileHasTime = (ncInput.get_var("time") != NULL);
	NcTimeDimension vecTimes;
	NcDim * dimTime = NULL;
	if (fFileHasTime) {
		ReadCFTimeDataFromNcFile(
			vecFiles[0],
			strInputData,
			vecTimes,
			false);

		WriteCFTimeDataToNcFile(
			&ncOutput,
			strOutputData,
			vecTimes);

		dimTime = ncOutput.get_dim("time");
		_ASSERT(dimTime != NULL);

	} else {
		vecTimes.push_back(Time(Time::CalendarNone));
	}

	// Build up the processing queue and define output variables
	varreg.ClearProcessingQueue();

	std::vector<NcVar *> vecNcVarOut(vecVarIxIn.size());
	std::vector<bool> vecNcVarOutFillValueCheck(vecVarIxIn.size(), false);
	for (int v = 0; v < vecVarIxIn.size(); v++) {
		NcVar * varOut = NULL;

		DimInfoVector vecAuxDimInfo;
		varreg.AppendVariableToProcessingQueue(
			vecFiles,
			grid,
			vecVarIxIn[v],
			&vecAuxDimInfo);

		std::vector<NcDim *> vecDimOut;

		if (fFileHasTime) {
			vecDimOut.push_back(dimTime);
		}

		for (long d = 0; d < vecAuxDimInfo.size(); d++) {
			NcDim * dim =
				AddNcDimOrUseExisting(
					ncOutput,
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

		vecNcVarOut[v] = ncOutput.add_var(
			vecVariableNamesOut[v].c_str(),
			ncFloat,
			vecDimOut.size(),
			const_cast<const NcDim**>(&(vecDimOut[0])));

		if (vecNcVarOut[v] == NULL) {
			_EXCEPTION2("Error creating variable \"%s\" in file \"%s\"",
				vecVariableNamesOut[v].c_str(),
				strOutputData.c_str());
		}
	}

	// Add FillValue attributes
	std::vector<float> vecVarFillValues(vecVarIxIn.size(), 0.0f);
	varreg.PopulateFillValues(vecFiles);
	for (int v = 0; v < vecNcVarOut.size(); v++) {
		const Variable & var = varreg.Get(vecVarIxIn[v]);
		if (var.HasExplicitFillValue()) {
			vecNcVarOut[v]->add_att("_FillValue", var.GetFillValueFloat());
			vecVarFillValues[v] = var.GetFillValueFloat();
		}
	}

	// Build the map from the input data to the shapefiles
	AnnounceStartBlock("Building map");

	std::vector<int> nShpMap;
	std::vector<size_t> sDataCount(sShpRegionCount, 0);

	nShpMap.resize(grid.m_dLat.GetRows(), static_cast<int>(-1));
	for (size_t k = 0; k < grid.m_dLat.GetRows(); k++) {
		double dLatDeg = RadToDeg(grid.m_dLat[k]);
		double dStandardLonDeg = LonDegToStandardRange(RadToDeg(grid.m_dLon[k]));
		for (size_t s = 0; s < sShpRegionCount; s++) {
			if (vecLatLonBox[s].contains(dLatDeg, dStandardLonDeg)) {
				if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg, dStandardLonDeg)) {
					nShpMap[k] = static_cast<int>(s);
					sDataCount[s]++;
					break;
				}
			}
		}
	}

	// Write the shape ids to a file
	if (fWriteShapeIds) {
		AnnounceStartBlock("Writing to file");

		NcVar * varShapeIds = NULL;
		if (grid.m_nGridDim.size() == 1) {
			varShapeIds = ncOutput.add_var("shapeid", ncInt, dim0);
			varShapeIds->add_att("_FillValue", (int)(-1));
			if (varShapeIds == NULL) {
				_EXCEPTION1("Unable to create variable \"shapeid\" in NetCDF file \"%s\"", strOutputData.c_str());
			}
			varShapeIds->set_cur((long)0);
			varShapeIds->put(&(nShpMap[0]), dim0->size());

		} else {
			varShapeIds = ncOutput.add_var("shapeid", ncInt, dim0, dim1);
			varShapeIds->add_att("_FillValue", (int)(-1));
			if (varShapeIds == NULL) {
				_EXCEPTION1("Unable to create variable \"shapeid\" in NetCDF file \"%s\"", strOutputData.c_str());
			}
			varShapeIds->set_cur(0, 0);
			varShapeIds->put(&(nShpMap[0]), dim0->size(), dim1->size());
		}

		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

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

			// Get the output position and size
			VariableAuxIndex lPos = varreg.GetProcessingQueueAuxIx();
			if (fFileHasTime) {
				lPos.insert(lPos.begin(), t);
			}
			VariableAuxIndex lSize(lPos.size(), 1);

			// Mask the data
			if (strOperation == "mask") {

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
					std::vector<float> dataOut(dataState.GetRows(), vecVarFillValues[v]);
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
				std::vector<float> dataOut(sShpRegionCount, vecVarFillValues[v]);

				if (fpCSV != NULL) {
					if (fFileHasTime) {
						fprintf(fpCSV, "%s, \"%s\"", vecVariableNamesOut[v].c_str(), vecTimes[t].ToString().c_str());
					} else {
						fprintf(fpCSV, "%s", vecVariableNamesOut[v].c_str());
					}
				}

				// Calculate the mean within all shapes
				if (strOperation == "mean") {
					for (size_t k = 0; k < nShpMap.size(); k++) {
						if ((!std::isnan(dataState[k])) && (dataState[k] != vecVarFillValues[v])) {
							if (nShpMap[k] != static_cast<size_t>(-1)) {
								dataOut[nShpMap[k]] += static_cast<double>(dataState[k]);
							}
						}
					}
					for (size_t s = 0; s < sShpRegionCount; s++) {
						if (sDataCount[s] != 0) {
							dataOut[s] /= static_cast<float>(sDataCount[s]);
						}
						if (fpCSV != NULL) {
							fprintf(fpCSV, ", %1.5f", dataOut[s]);
						}
					}
	
				// Generate arrays of data in each shapefile
				} else {
					std::vector< std::vector<float> > vecShpData;
					vecShpData.resize(sShpRegionCount);
	
					for (size_t k = 0; k < nShpMap.size(); k++) {
						if ((!std::isnan(dataState[k])) && (dataState[k] != vecVarFillValues[v])) {
							if (nShpMap[k] != static_cast<size_t>(-1)) {
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
/*
						if (s == 0) {
							FILE * fpx = fopen("test.txt", "w");
							for (size_t i = 0; i < vecShpData[s].size(); i++) {
								fprintf(fpx, "%1.5f\n", vecShpData[s][i]);
							}
							fclose(fpx);
						}
*/
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
								fprintf(fpCSV, ", %1.5f", dataOut[s]);
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
								fprintf(fpCSV, ", %1.5f", dataOut[s]);
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
								fprintf(fpCSV, ", %1.5f", dataOut[s]);
							}
						}
					}
				}

				vecNcVarOut[v]->set_cur(&(lPos[0]));
				vecNcVarOut[v]->put(&(dataOut[0]), &(lSize[0]));

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
/*
// Loop through all variables
	AnnounceStartBlock("Processing data");

	for (size_t v = 0; v < vecVariableNamesIn.size(); v++) {
		AnnounceStartBlock("Variable \"%s\"", vecVariableNamesIn[v].c_str());

		bool fFirstRecordDim = true;

		// Get data values
		NcVar * varData = ncInput.get_var(vecVariableNamesIn[v].c_str());
		if (varData == NULL) {
			_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
				strInputData.c_str(), vecVariableNamesIn[v].c_str());
		}
		if ((varData->num_dims() < 2) || (varData->num_dims() > 3)) {
			_EXCEPTION2("File \"%s\" variable \"%s\" must contain two (time, ncol) or three dimensions (time, lat, lon)",
				strInputData.c_str(), vecVariableNamesIn[v].c_str());
		}
		if (varData->num_dims() == 3) {
			if (varData->get_dim(1)->size() != dLatDeg.size()) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have latitude as its second dimension",
					strInputData.c_str(), vecVariableNamesIn[v].c_str());
			}
			if (varData->get_dim(2)->size() != dLonDeg.size()) {
				_EXCEPTION2("File \"%s\" variable \"%s\" must have longitude as its third dimension",
					strInputData.c_str(), vecVariableNamesIn[v].c_str());
			}
		}
		if (varData->num_dims() == 2) {
			if ((varData->get_dim(0)->size() == dLatDeg.size()) &&
			    (varData->get_dim(1)->size() == dLonDeg.size())
			) {
				fFirstRecordDim = false;

			} else {
				if (varData->get_dim(1)->size() != dLatDeg.size()) {
					_EXCEPTION2("File \"%s\" variable \"%s\" must have a second dimension with the same length as \"lat\"",
						strInputData.c_str(), vecVariableNamesIn[v].c_str());
				}
				if (dLonDeg.size() != dLatDeg.size()) {
					_EXCEPTION1("File \"%s\" dimension \"lat\" must be the same length as dimension \"lon\"",
						strInputData.c_str());
				}
			}
		}

		// Number of time slices on this variable
		size_t sTimeCount = 1;

		// Copy first dimension if it is record dimension
		NcDim * dimVarFirst = varData->get_dim(0);
		_ASSERT(dimVarFirst != NULL);

		NcDim * dimVarFirstOut = NULL;

		if (fFirstRecordDim) {
			sTimeCount = dimVarFirst->size();

			dimVarFirstOut = ncOutput.get_dim(dimVarFirst->name());
			if (dimVarFirstOut == NULL) {
				dimVarFirstOut = ncOutput.add_dim(dimVarFirst->name(), dimVarFirst->size());
				if (dimVarFirstOut == NULL) {
					_EXCEPTION2("Error creating dimension \"%s\" in file \"%s\"",
						dimVarFirst->name(), strOutputData.c_str());
				}
				CopyNcVarIfExists(ncInput, ncOutput, dimVarFirst->name());
			}
		}

		// Allocate data and count
		std::vector<float> dData;

		if (varData->num_dims() == 3) {
			dData.resize(dLatDeg.size() * dLonDeg.size());
		}
		if (varData->num_dims() == 2) {
			if (fFirstRecordDim) {
				dData.resize(dLatDeg.size());
			} else {
				dData.resize(dLatDeg.size() * dLonDeg.size());
			}
		}

		// For the first variable, create the map
		if (v == 0) {
			AnnounceStartBlock("Creating map");

			// Rectilinear data
			if ((varData->num_dims() == 3) || (!fFirstRecordDim)) {
				nShpMap.resize(dLatDeg.size() * dLonDeg.size(), static_cast<size_t>(-1));
				for (size_t j = 0; j < dLatDeg.size(); j++) {
				for (size_t i = 0; i < dLonDeg.size(); i++) {
					size_t k = j * dLonDeg.size() + i;
		
					for (size_t s = 0; s < sShpRegionCount; s++) {
						double dStandardLonDeg = LonDegToStandardRange(dLonDeg[i]);
						if (vecLatLonBox[s].contains(dLatDeg[j], dStandardLonDeg)) {
							if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg[j], dStandardLonDeg)) {
								nShpMap[k] = s;
								break;
							}
						}
					}
				}
				}

			// Unstructured data
			} else {
				nShpMap.resize(dLatDeg.size(), static_cast<size_t>(-1));
				for (size_t k = 0; k < dLatDeg.size(); k++) {
					for (size_t s = 0; s < sShpRegionCount; s++) {
						double dStandardLonDeg = LonDegToStandardRange(dLonDeg[k]);
						if (vecLatLonBox[s].contains(dLatDeg[k], dStandardLonDeg)) {
							if (FaceContainsNode(mesh.faces[s], mesh.nodes, dLatDeg[k], dStandardLonDeg)) {
								nShpMap[k] = s;
								break;
							}
						}
					}
				}
			}

			AnnounceEndBlock("Done");
		}

		// Create output variable
		NcVar * varDataOut = NULL;
		if (fFirstRecordDim) {
			varDataOut = ncOutput.add_var(vecVariableNames[v].c_str(), ncFloat, dimVarFirstOut, dimShp);
		} else {
			varDataOut = ncOutput.add_var(vecVariableNames[v].c_str(), ncFloat, dimShp);
		}

		if (varDataOut == NULL) {
			_EXCEPTION1("Error creating variable \"%s\" in output file", vecVariableNames[v].c_str());
		}
		CopyNcVarAttributes(varData, varDataOut);

		// Fillvalue
		float dFillValue = FLT_MAX;
		NcAtt * attFillValue = varData->get_att("_FillValue");
		if (attFillValue != NULL) {
			dFillValue = attFillValue->as_float(0);
		}

		// Loop through all times
		for (size_t t = 0; t < sTimeCount; t++) {
	
			Announce("Time %lu/%lu", t, sTimeCount);

			// Output to CSV
			if (fpCSV != NULL) {
				if (sTimeCount == 1) {
					fprintf(fpCSV, "%s", vecVariableNames[v].c_str());
				} else {
					fprintf(fpCSV, "%s(%lu)", vecVariableNames[v].c_str(), t);
				}
			}

			// Reset the counts
			std::fill(sDataCount.begin(), sDataCount.end(), 0);
	
			// Allocate output data
			std::vector<float> dataOut(sShpRegionCount, 0.0f);
	
			// Load data
			if (varData->num_dims() == 3) {
				varData->set_cur(t,0,0);
				varData->get(&(dData[0]), 1, dLatDeg.size(), dLonDeg.size());
			} else {
				if (fFirstRecordDim) {
					varData->set_cur(t,0);
					varData->get(&(dData[0]), 1, dLatDeg.size());
				} else {
					varData->set_cur(0,0);
					varData->get(&(dData[0]), dLatDeg.size(), dLonDeg.size());
				}
			}
	
			// Calculate the mean within all shapes
			if (strOperation == "mean") {
				for (size_t k = 0; k < nShpMap.size(); k++) {
					if ((!std::isnan(dData[k])) && (dData[k] != dFillValue)) {
						if (nShpMap[k] != static_cast<size_t>(-1)) {
							dataOut[nShpMap[k]] += static_cast<double>(dData[k]);
							sDataCount[nShpMap[k]]++;
						}
					}
				}
				for (size_t s = 0; s < sShpRegionCount; s++) {
					if (sDataCount[s] != 0) {
						dataOut[s] /= static_cast<float>(sDataCount[s]);
					}
				}

			// Generate arrays of data in each shapefile
			} else {
				std::vector< std::vector<float> > vecShpData;
				vecShpData.resize(sShpRegionCount);

				for (size_t k = 0; k < nShpMap.size(); k++) {
					if ((!std::isnan(dData[k])) && (dData[k] != dFillValue)) {
						if (nShpMap[k] != static_cast<size_t>(-1)) {
							vecShpData[nShpMap[k]].push_back(dData[k]);
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

					if (s == 0) {
						FILE * fpx = fopen("test.txt", "w");
						for (size_t i = 0; i < vecShpData[s].size(); i++) {
							fprintf(fpx, "%1.5f\n", vecShpData[s][i]);
						}
						fclose(fpx);
					}

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
							fprintf(fpCSV, ", %1.5f", dataOut[s]);
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
							fprintf(fpCSV, ", %1.5f", dataOut[s]);
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
							fprintf(fpCSV, ", %1.5f", dataOut[s]);
						}
					}
				}
			}

			// Write to disk
			if (fFirstRecordDim) {
				varDataOut->set_cur(t,0);
				varDataOut->put(&(dataOut[0]), 1, sShpRegionCount);
			} else {
				varDataOut->set_cur((long)0);
				varDataOut->put(&(dataOut[0]), sShpRegionCount);
			}

			if (fpCSV != NULL) {
				fprintf(fpCSV, "\n");
			}
		}

		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");
*/
	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

}

///////////////////////////////////////////////////////////////////////////////

