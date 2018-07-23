///////////////////////////////////////////////////////////////////////////////
///
///	\file    SpineARs_v2.cpp
///	\author  Paul Ullrich
///	\version July 17th, 2018
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

#include "Variable.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "SimpleGrid.h"

#include "DataVector.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <set>
#include <queue>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a thresholding operator.
///	</summary>
class ThresholdOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		GreaterThan,
		LessThan,
		GreaterThanEqualTo,
		LessThanEqualTo,
		EqualTo,
		NotEqualTo
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Distance,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Parse variable
		Variable var;
		int iLast = var.ParseFromString(varreg, strOp) + 1;
		m_varix = varreg.FindOrRegister(var);

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = NotEqualTo;
					} else {
						_EXCEPTION1("Threshold invalid operation \"%s\"",
							strSubStr.c_str());
					}

					iLast = i + 1;
					eReadMode = ReadMode_Value;

				// Read in value
				} else if (eReadMode == ReadMode_Value) {
					m_dValue = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in minimum count
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
							"\nRequired: \"<name>,<operation>"
							",<value>,<distance>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>,<distance>\"",
					strOp.c_str());
		}

		if (m_dDistance < 0.0) {
			_EXCEPTIONT("For threshold op, distance must be nonnegative");
		}

		// Output announcement
		std::string strDescription = var.ToString(varreg);
		if (m_eOp == GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == NotEqualTo) {
			strDescription += " is not equal to ";
		}

		char szBuffer[128];
		sprintf(szBuffer, "%f within %f degrees",
			m_dValue, m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for thresholding.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Distance to search for threshold value
	///	</summary>
	double m_dDistance;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field satisfies the threshold.
///	</summary>
template <typename real>
bool SatisfiesThreshold(
	const SimpleGrid & grid,
	const DataVector<real> & dataState,
	const int ix0,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const double dMaxDist
) {
	// Special case if dMaxDist is zero
	if (dMaxDist < 1.0e-12) {
		double dValue = dataState[ix0];

		if (op == ThresholdOp::GreaterThan) {
			if (dValue > dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThan) {
			if (dValue < dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThanEqualTo) {
			if (dValue <= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::EqualTo) {
			if (dValue == dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::NotEqualTo) {
			if (dValue != dTargetValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

	}

	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ix0);

	// Set of nodes that have already been visited
	std::set<int> setNodesVisited;

	// Latitude and longitude at the origin
	double dLat0 = grid.m_dLat[ix0];
	double dLon0 = grid.m_dLon[ix0];

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		double dR =
			sin(dLat0) * sin(dLatThis)
			+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0);

		if (dR >= 1.0) {
			dR = 0.0;
		} else if (dR <= -1.0) {
			dR = 180.0;
		} else {
			dR = 180.0 / M_PI * acos(dR);
		}
		if (dR != dR) {
			_EXCEPTIONT("NaN value detected");
		}

		if ((ix != ix0) && (dR > dMaxDist)) {
			continue;
		}

		// Value at this location
		double dValue = dataState[ix];

		// Apply operator
		if (op == ThresholdOp::GreaterThan) {
			if (dValue > dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThan) {
			if (dValue < dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::LessThanEqualTo) {
			if (dValue <= dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::EqualTo) {
			if (dValue == dTargetValue) {
				return true;
			}

		} else if (op == ThresholdOp::NotEqualTo) {
			if (dValue != dTargetValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		// Special case: zero distance
		if (dMaxDist == 0.0) {
			return false;
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse the list of input files.
///	</summary>
void ParseInputFiles(
	const std::string & strInputFile,
	std::vector<NcFile *> & vecFiles
) {
	int iLast = 0;
	for (int i = 0; i <= strInputFile.length(); i++) {
		if ((i == strInputFile.length()) ||
		    (strInputFile[i] == ';')
		) {
			std::string strFile =
				strInputFile.substr(iLast, i - iLast);

			NcFile * pNewFile = new NcFile(strFile.c_str());

			if (!pNewFile->is_valid()) {
				_EXCEPTION1("Cannot open input file \"%s\"",
					strFile.c_str());
			}

			vecFiles.push_back(pNewFile);
			iLast = i+1;
		}
	}

	if (vecFiles.size() == 0) {
		_EXCEPTION1("No input files found in \"%s\"",
			strInputFile.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

class SpineARsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	SpineARsParam() :
		fpLog(NULL),
		nLaplacianPoints(0),
		dLaplacianDist(0.0),
		dMinLaplacian(0.0),
		dMinAbsLat(0.0),
		fOutputLaplacian(false),
		fRegional(false),
		iVerbosityLevel(0),
		pvecThresholdOp(NULL)
	{ }

public:
	// Log
	FILE * fpLog;

	// Input integrated water vapor variable name
	//std::string strIWVVariable;

	// Number of points in Laplacian calculation
	int nLaplacianPoints;

	// Radius of the Laplacian (in degrees)
	double dLaplacianDist;

	// Minimum Laplacian
	double dMinLaplacian;

	// Minimum absolute latitude (in degrees)
	double dMinAbsLat;

	// Minimum pointwise integrated water vapor
	//double dMinIWV;

	// Output Laplacian
	bool fOutputLaplacian;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Verbosity level
	int iVerbosityLevel;

	// Vector of threshold operators
	std::vector<ThresholdOp> * pvecThresholdOp;
};

///////////////////////////////////////////////////////////////////////////////

void SpineARs(
	int iFile,
	const std::string & strInputFiles,
	const std::string & strOutputFile,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const SpineARsParam & param
) {
	// Set the Announce buffer
	if (param.fpLog == NULL) {
		_EXCEPTIONT("Invalid log buffer");
	}

	AnnounceSetOutputBuffer(param.fpLog);
	AnnounceOutputOnAllRanks();

	// Dereference pointers to operators
	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	// Unload data from the VariableRegistry
	varreg.UnloadAllGridData();

	// Define the SimpleGrid
	SimpleGrid grid;

	// Dimensions
	int nSize = 0;
	int nLon = 0;
	int nLat = 0;

	// Load in the benchmark file
	NcFileVector vecFiles;

	ParseInputFiles(strInputFiles, vecFiles);

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

		nSize = grid.GetSize();

	// No connectivity file; check for latitude/longitude dimension
	} else {

		NcDim * dimLat = vecFiles[0]->get_dim("lat");
		if (dimLat == NULL) {
			_EXCEPTIONT("No dimension \"lat\" found in input file");
		}

		NcDim * dimLon = vecFiles[0]->get_dim("lon");
		if (dimLon == NULL) {
			_EXCEPTIONT("No dimension \"lon\" found in input file");
		}

		NcVar * varLat = vecFiles[0]->get_var("lat");
		if (varLat == NULL) {
			_EXCEPTIONT("No variable \"lat\" found in input file");
		}

		NcVar * varLon = vecFiles[0]->get_var("lon");
		if (varLon == NULL) {
			_EXCEPTIONT("No variable \"lon\" found in input file");
		}

		nLat = dimLat->size();
		nLon = dimLon->size();

		DataVector<double> vecLat(nLat);
		varLat->get(vecLat, nLat);

		for (int j = 0; j < nLat; j++) {
			vecLat[j] *= M_PI / 180.0;
		}

		DataVector<double> vecLon(nLon);
		varLon->get(vecLon, nLon);

		for (int i = 0; i < nLon; i++) {
			vecLon[i] *= M_PI / 180.0;
		}

		// Generate the SimpleGrid
		grid.GenerateLatitudeLongitude(vecLat, vecLon, param.fRegional);
	}

	// Get time dimension
	NcDim * dimTime = vecFiles[0]->get_dim("time");
	NcVar * varTime = NULL;
	if (dimTime != NULL) {
		varTime = vecFiles[0]->get_var("time");
	}

	int nTime = 1;
	if (dimTime != NULL) {
		nTime = dimTime->size();
	}

	DataVector<double> dTime;
	dTime.Initialize(nTime);

	if (varTime != NULL) {
		if (varTime->type() == ncDouble) {
			varTime->get(dTime, nTime);

		} else if (varTime->type() == ncFloat) {
			DataVector<float> dTimeFloat;
			dTimeFloat.Initialize(nTime);

			varTime->get(dTimeFloat, nTime);
			for (int t = 0; t < nTime; t++) {
				dTime[t] = static_cast<double>(dTimeFloat[t]);
			}

		} else if (varTime->type() == ncInt) {
			DataVector<int> dTimeInt;
			dTimeInt.Initialize(nTime);

			varTime->get(dTimeInt, nTime);
			for (int t = 0; t < nTime; t++) {
				dTime[t] = static_cast<double>(dTimeInt[t]);
			}

		} else {
			_EXCEPTIONT("Variable \"time\" has an invalid type:\n"
				"Expected \"float\", \"double\" or \"int\"");
		}

	} else {
		for (int t = 0; t < nTime; t++) {
			dTime[t] = static_cast<double>(t);
		}
	}

	// Create reference to NetCDF input file
	NcFile & ncInput = *(vecFiles[0]);

	// Open the NetCDF output file
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open NetCDF file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Copy over latitude, longitude and time variables to output file
	NcDim * dimTimeOut = NULL;
	if ((dimTime != NULL) && (varTime != NULL)) {
		CopyNcVar(ncInput, ncOutput, "time", true);
		dimTimeOut = ncOutput.get_dim("time");
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

	} /*else if (nAddTimeDim != -1) {
		dimTimeOut = ncOutput.add_dim("time", 0);
		if (dimTimeOut == NULL) {
			_EXCEPTIONT("Error creating dimension \"time\" in output file");
		}

		NcVar * varTimeOut = ncOutput.add_var("time", ncDouble, dimTimeOut);
		if (varTimeOut == NULL) {
			_EXCEPTIONT("Error copying variable \"time\" to output file");
		}

		double dTime = static_cast<double>(nAddTimeDim);
		varTimeOut->set_cur((long)0);
		varTimeOut->put(&dTime, 1);

		if (strAddTimeDimUnits != "") {
			varTimeOut->add_att("units", strAddTimeDimUnits.c_str());
		}
		varTimeOut->add_att("long_name", "time");
		varTimeOut->add_att("calendar", "standard");
		varTimeOut->add_att("standard_name", "time");
	}*/

	// FIX for unstructured grids
	CopyNcVar(ncInput, ncOutput, "lat", true);
	CopyNcVar(ncInput, ncOutput, "lon", true);

	NcDim * dim0 = ncOutput.get_dim("lat");
	if (dim0 == NULL) {
		_EXCEPTIONT("Error copying variable \"lat\" to output file");
	}
	NcDim * dim1 = ncOutput.get_dim("lon");
	if (dim1 == NULL) {
		_EXCEPTIONT("Error copying variable \"lon\" to output file");
	}

	// Create output variables
	NcVar * varIWVtag = NULL;
	NcVar * varLaplacian = NULL;

	if (grid.m_nGridDim.size() == 1) {
		if (dimTimeOut != NULL) {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dimTimeOut,
				dim0);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dimTimeOut,
				dim0);

		} else {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dim0);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dim0);
		}

	} else if (grid.m_nGridDim.size() == 2) {
		if (dimTimeOut != NULL) {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dimTimeOut,
				dim0,
				dim1);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dimTimeOut,
				dim0,
				dim1);

		} else {
			varIWVtag = ncOutput.add_var(
				"ar_binary_tag",
				ncByte,
				dim0,
				dim1);

			varLaplacian = ncOutput.add_var(
				"laplacian",
				ncDouble,
				dim0,
				dim1);
		}

	} else {
		_EXCEPTIONT("Invalid grid dimension -- value must be 1 or 2");
	}

	AnnounceEndBlock("Done");

	// Tagged cell array
	DataVector<int> bIWVtag(grid.GetSize());

	// REMOVE: Laplacian
	DataVector<double> dLaplacian(grid.GetSize());

	// Loop through all times
	for (int t = 0; t < nTime; t ++) {

		// Announce
		char szBuffer[20];
		sprintf(szBuffer, "Time %i", t);
		AnnounceStartBlock(szBuffer);
/*
		// Load in data at this time slice
		AnnounceStartBlock("Reading data");
		if (dimTime != NULL) {
			if (grid.m_nGridDim.size() == 1) {
				varIWV->set_cur(t, 0);
				varIWV->get(&(dIWV[0]), 1, grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varIWV->set_cur(t, 0, 0);
				varIWV->get(&(dIWV[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}

		} else {
			if (grid.m_nGridDim.size() == 1) {
				varIWV->set_cur((long)0);
				varIWV->get(&(dIWV[0]), grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varIWV->set_cur(0, 0);
				varIWV->get(&(dIWV[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}
		}
		AnnounceEndBlock("Done");

		// Compute Laplacian
		AnnounceStartBlock("Applying Laplacian");
		opLaplacian.Apply(dIWV, dLaplacian);
		AnnounceEndBlock("Done");
*/
/*
		// DEBUG: Output Laplacian range
		double dMinLaplacian = dLaplacian[0];
		double dMaxLaplacian = dLaplacian[0];
		for (int i = 0; i < dLaplacian.GetRows(); i++) {
			if (dLaplacian[i] < dMinLaplacian) {
				dMinLaplacian = dLaplacian[i];
			}
			if (dLaplacian[i] > dMaxLaplacian) {
				dMaxLaplacian = dLaplacian[i];
			}
		}
		Announce("Laplacian Range: %1.5e %1.5e", dMinLaplacian, dMaxLaplacian);
*/

		AnnounceStartBlock("Build tagged cell array");
		bIWVtag.Zero();
		for (int i = 0; i < grid.GetSize(); i++) {
			if (fabs(grid.m_dLat[i]) < param.dMinAbsLat * M_PI / 180.0) {
				continue;
			}

/*
			if (dIWV[i] < dMinIWV) {
				continue;
			}
*/
/*
			if (dLaplacian[i] > -dMinLaplacian) {
				continue;
			}
*/
			bIWVtag[i] = 1;
		}
		AnnounceEndBlock("Done");

		// Eliminate based on threshold commands
		AnnounceStartBlock("Apply threshold commands");
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			// Load the search variable data
			Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataVector<float> & dataState = var.GetData();

			// Loop through data
			for (int i = 0; i < grid.GetSize(); i++) {
				if (bIWVtag[i] == 0) {
					continue;
				}

				// Determine if the threshold is satisfied
				bool fSatisfiesThreshold =
					SatisfiesThreshold<float>(
						grid,
						dataState,
						i,
						vecThresholdOp[tc].m_eOp,
						vecThresholdOp[tc].m_dValue,
						vecThresholdOp[tc].m_dDistance
					);

				if (!fSatisfiesThreshold) {
					bIWVtag[i] = 0;
				}
			}
		}
		AnnounceEndBlock("Done");

		// Output tagged cell array
		AnnounceStartBlock("Writing results");
		if (dimTimeOut != NULL) {
			if (varLaplacian != NULL) {
				if (grid.m_nGridDim.size() == 1) {
					varLaplacian->set_cur(t, 0);
					varLaplacian->put(&(dLaplacian[0]), 1, grid.m_nGridDim[0]);
				} else if (grid.m_nGridDim.size() == 2) {
					varLaplacian->set_cur(t, 0, 0);
					varLaplacian->put(&(dLaplacian[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
				} else {
					_EXCEPTION();
				}
			}

			if (grid.m_nGridDim.size() == 1) {
				varIWVtag->set_cur(t, 0);
				varIWVtag->put(&(bIWVtag[0]), 1, grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {

				varIWVtag->set_cur(t, 0, 0);
				varIWVtag->put(&(bIWVtag[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}

		} else {
			if (varLaplacian != NULL) {
				if (grid.m_nGridDim.size() == 1) {
					varLaplacian->set_cur((long)0);
					varLaplacian->put(&(dLaplacian[0]), grid.m_nGridDim[0]);
				} else if (grid.m_nGridDim.size() == 2) {
					varLaplacian->set_cur(0, 0);
					varLaplacian->put(&(dLaplacian[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
				} else {
					_EXCEPTION();
				}
			}

			if (grid.m_nGridDim.size() == 1) {
				varIWVtag->set_cur((long)0);
				varIWVtag->put(&(bIWVtag[0]), grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varIWVtag->set_cur(0, 0);
				varIWVtag->put(&(bIWVtag[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}
		}

		AnnounceEndBlock("Done");

		AnnounceEndBlock(NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);
#endif

	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {
	// Parameters for SpineARs
	SpineARsParam sarparam;

	// Input dat file
	std::string strInputFile;

	// Input list file
	std::string strInputFileList;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutput;

	// Output file list
	std::string strOutputFileList;

	// Output file list
	std::string strThresholdCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutput, "out", "");
		CommandLineString(strOutputFileList, "out_file_list", "");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		//CommandLineString(sarparam.strIWVVariable, "var", "");
		CommandLineInt(sarparam.nLaplacianPoints, "laplacianpoints", 8);
		CommandLineDoubleD(sarparam.dLaplacianDist, "laplaciandist", 10.0, "(degrees)");
		CommandLineDouble(sarparam.dMinLaplacian, "minlaplacian", 0.5e4);
		CommandLineDouble(sarparam.dMinAbsLat, "minabslat", 15.0);
		//CommandLineDouble(sarparam.dMinIWV, "minval", 20.0);
		//CommandLineInt(nMinArea, "minarea", 0);
		//CommandLineInt(nAddTimeDim, "addtimedim", -1);
		//CommandLineString(strAddTimeDimUnits, "addtimedimunits", "");
		CommandLineBool(sarparam.fOutputLaplacian, "laplacianout");
		CommandLineBool(sarparam.fRegional, "regional");
		CommandLineInt(sarparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(sarparam.iVerbosityLevel);

	// Check input
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}

	// Check output
	if ((strOutput.length() != 0) && (strOutputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--out) or (--out_data_list)"
			" may be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputFiles;

	if (strInputFile.length() != 0) {
		vecInputFiles.push_back(strInputFile);

	} else {
		std::ifstream ifInputFileList(strInputFileList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strInputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputFiles.push_back(strFileLine);
		}
	}

	// Load output file list
	std::vector<std::string> vecOutputFiles;

	if (strOutputFileList.length() != 0) {

		std::ifstream ifOutputFileList(strOutputFileList.c_str());
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strOutputFileList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputFiles.push_back(strFileLine);
		}

		if (vecOutputFiles.size() != vecInputFiles.size()) {
			_EXCEPTIONT("File --in_file_list must match --out_file_list");
		}
	}

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;
	sarparam.pvecThresholdOp = &vecThresholdOp;

	if (strThresholdCmd != "") {
		AnnounceStartBlock("Parsing threshold operations");

		int iLast = 0;
		for (int i = 0; i <= strThresholdCmd.length(); i++) {

			if ((i == strThresholdCmd.length()) ||
				(strThresholdCmd[i] == ';') ||
				(strThresholdCmd[i] == ':')
			) {
				std::string strSubStr =
					strThresholdCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecThresholdOp.size());
				vecThresholdOp.resize(iNextOp + 1);
				vecThresholdOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}
/*
	// Build Laplacian operator
	AnnounceStartBlock("Building Laplacian operator");

	SparseMatrix<double> opLaplacian;
	BuildLaplacianOperator(
		grid,
		nLaplacianPoints,
		dLaplacianDist,
		opLaplacian);

	AnnounceEndBlock("Done");
*/
#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
#endif

	AnnounceStartBlock("Begin search operation");
	if (vecInputFiles.size() != 1) {
		if (vecOutputFiles.size() != 0) {
			Announce("Output will be written following --out_file_list");
		} else if (strOutput == "") {
			Announce("Output will be written to outXXXXXX.dat");
		} else {
			Announce("Output will be written to %sXXXXXX.dat",
				strOutput.c_str());
		}
		Announce("Logs will be written to logXXXXXX.txt");
	}

	// Loop over all files to be processed
	for (int f = 0; f < vecInputFiles.size(); f++) {
#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
#endif
		// Generate output file name
		std::string strOutputFile;
		if (vecInputFiles.size() == 1) {
			sarparam.fpLog = stdout;

			if (strOutput == "") {
				strOutputFile = "out.dat";
			} else {
				strOutputFile = strOutput;
			}

		} else {
			char szFileIndex[32];
			sprintf(szFileIndex, "%06i", f);

			if (vecOutputFiles.size() != 0) {
				strOutputFile = vecOutputFiles[f];
			} else {
				if (strOutput == "") {
					strOutputFile =
						"out" + std::string(szFileIndex) + ".dat";
				} else {
					strOutputFile =
						strOutput + std::string(szFileIndex) + ".dat";
				}
			}

			std::string strLogFile = "log" + std::string(szFileIndex) + ".txt";
			sarparam.fpLog = fopen(strLogFile.c_str(), "w");
		}

		// Perform SpineARs
		SpineARs(
			f,
			vecInputFiles[f],
			strOutputFile,
			strConnectivity,
			varreg,
			sarparam);

		// Close the log file
		if (vecInputFiles.size() != 1) {
			fclose(sarparam.fpLog);
		}
	}

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

}

///////////////////////////////////////////////////////////////////////////////

