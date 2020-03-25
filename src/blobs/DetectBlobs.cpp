///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectBlobs.cpp
///	\author  Paul Ullrich
///	\version July 17th, 2018
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

#include "Variable.h"
#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "SimpleGrid.h"

#include "DataArray1D.h"
#include "DataArray2D.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include "../nodes/ThresholdOp.h"

#include <set>
#include <queue>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, int> Point;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field satisfies the threshold.
///	</summary>
template <typename real>
bool SatisfiesThreshold(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
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
///		A class storing a filtering operator.
///	</summary>
class FilterOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FilterOp() :
		m_varix(InvalidVariableIndex),
		m_eOp(ThresholdOp::GreaterThan),
		m_dValue(0.0),
		m_nCount(0)
	{ }

public:
	///	<summary>
	///		Parse a filter operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Value,
			ReadMode_Count,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					if (strSubStr == ">") {
						m_eOp = ThresholdOp::GreaterThan;
					} else if (strSubStr == "<") {
						m_eOp = ThresholdOp::LessThan;
					} else if (strSubStr == ">=") {
						m_eOp = ThresholdOp::GreaterThanEqualTo;
					} else if (strSubStr == "<=") {
						m_eOp = ThresholdOp::LessThanEqualTo;
					} else if (strSubStr == "=") {
						m_eOp = ThresholdOp::EqualTo;
					} else if (strSubStr == "!=") {
						m_eOp = ThresholdOp::NotEqualTo;
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
					eReadMode = ReadMode_Count;

				// Read in minimum count
				} else if (eReadMode == ReadMode_Count) {
					m_nCount = atoi(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
							"\nRequired: \"<name>,<operation>"
							",<value>,<count>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in filter op \"%s\""
					"\nRequired: \"<name>,<operation>,<value>,<count>\"",
					strOp.c_str());
		}

		if (m_nCount < 1) {
			_EXCEPTIONT("For filter op, count must be positive");
		}

		// Output announcement
		std::string strDescription = varreg.GetVariableString(m_varix);

		if (m_eOp == ThresholdOp::GreaterThan) {
			strDescription += " is greater than ";
		} else if (m_eOp == ThresholdOp::LessThan) {
			strDescription += " is less than ";
		} else if (m_eOp == ThresholdOp::GreaterThanEqualTo) {
			strDescription += " is greater than or equal to ";
		} else if (m_eOp == ThresholdOp::LessThanEqualTo) {
			strDescription += " is less than or equal to ";
		} else if (m_eOp == ThresholdOp::EqualTo) {
			strDescription += " is equal to ";
		} else if (m_eOp == ThresholdOp::NotEqualTo) {
			strDescription += " is not equal to ";
		}

		Announce("%s %f at at least %i points",
			strDescription.c_str(),
			m_dValue,
			m_nCount);
	}

public:
	///	<summary>
	///		Variable to use for filtering.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Operation that must be satisfied.
	///	</summary>
	ThresholdOp::Operation m_eOp;

	///	<summary>
	///		Filter threshold value.
	///	</summary>
	double m_dValue;

	///	<summary>
	///		Number of points that exceed the threshold.
	///	</summary>
	int m_nCount;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a output operator.
///	</summary>
class BlobOutputOp {

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
			ReadMode_Name,
			ReadMode_Invalid
		} eReadMode = ReadMode_Name;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in name
				if (eReadMode == ReadMode_Name) {
					m_strName = strSubStr;
					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nInsufficient entries in output op \"%s\""
							"\nRequired: \"<variable>,<name>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in output op \"%s\""
					"\nRequired: \"<variable>,<name>\"",
					strOp.c_str());
		}

		// Output announcement
		std::string strDescription = varreg.GetVariableString(m_varix);
		strDescription += " with name \"";
		strDescription += m_strName;
		strDescription += "\"";

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for output.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Name of variable in NetCDF output file.
	///	</summary>
	std::string m_strName;
};

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void ApplyFilters(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const ThresholdOp::Operation op,
	const double dTargetValue,
	const int nCount,
	DataArray1D<int> & bTag
) {
	_ASSERT(bTag.GetRows() == grid.GetSize());
	_ASSERT(bTag.GetRows() == dataState.GetRows());

	// Number of blobs
	int nBlobs = 0;
	int nBlobsFiltered = 0;

	// Visit all tagged points
	std::set<int> setNodesVisited;
	for (int i = 0; i < grid.GetSize(); i++) {

		// Verify point it tagged and we haven't visited before
		if (bTag[i] == 0) {
			continue;
		}
		if (setNodesVisited.find(i) != setNodesVisited.end()) {
			continue;
		}

		// Number of blobs
		nBlobs++;

		// New blob
		std::set<int> setCurrentBlob;

		// Number of points within blob that satisfy threshold
		int nThresholdPoints = 0;

		// Build connectivity
		std::queue<int> queueToVisit;
		queueToVisit.push(i);

		while (queueToVisit.size() != 0) {
			int iNext = queueToVisit.front();
			queueToVisit.pop();

			// Verify point is tagged and we haven't visited before
			if (bTag[iNext] == 0) {
				continue;
			}
			if (setNodesVisited.find(iNext) != setNodesVisited.end()) {
				continue;
			}
			setNodesVisited.insert(iNext);
			setCurrentBlob.insert(iNext);

			// Check if this point satisfies threshold
			bool fSatisfiesThreshold =
				SatisfiesThreshold<real>(
					grid,
					dataState,
					iNext,
					op,
					dTargetValue,
					0.0);

			if (fSatisfiesThreshold) {
				nThresholdPoints++;
			}

			// Insert all connected neighbors into "to visit" queue
			for (int n = 0; n < grid.m_vecConnectivity[iNext].size(); n++) {
				queueToVisit.push(grid.m_vecConnectivity[iNext][n]);
			}
		}

		// If not enough points satisfy the filter then eliminate this blob
		if (nThresholdPoints < nCount) {
			nBlobsFiltered++;
			for (auto it = setCurrentBlob.begin(); it != setCurrentBlob.end(); it++) {
				bTag[*it] = 0;
			}
		}
	}

	// Announce results
	Announce("Filter removed %i of %i blobs", nBlobsFiltered, nBlobs);
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

class DetectBlobsParam {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	DetectBlobsParam() :
		fpLog(NULL),
		dMinAbsLat(0.0),
		dMinLat(-90.0),
		dMaxLat(90.0),
		fRegional(false),
		iVerbosityLevel(0),
		strTagVar("binary_tag"),
		strLongitudeName("lon"),
		strLatitudeName("lat"),
		pvecThresholdOp(NULL),
		pvecFilterOp(NULL),
		pvecOutputOp(NULL)
	{ }

public:
	// Log
	FILE * fpLog;

	// Minimum absolute latitude (in degrees)
	double dMinAbsLat;

	// Minimum latitude (in degrees)
	double dMinLat;

	// Maximum latitude (in degrees)
	double dMaxLat;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Verbosity level
	int iVerbosityLevel;

	// Name of output variable for tag
	std::string strTagVar;

	// Name of longitude variabe
	std::string strLongitudeName;

	// Name of latitude variable
	std::string strLatitudeName;

	// Vector of threshold operators
	std::vector<ThresholdOp> * pvecThresholdOp;

	// Vector of filter operators
	std::vector<FilterOp> * pvecFilterOp;

	// Vector of output operators
	std::vector<BlobOutputOp> * pvecOutputOp;
};

///////////////////////////////////////////////////////////////////////////////

void DetectBlobs(
	int iFile,
	const std::string & strInputFiles,
	const std::string & strOutputFile,
	const std::string & strConnectivity,
	VariableRegistry & varreg,
	const DetectBlobsParam & param
) {
	// Set the Announce buffer
	if (param.fpLog == NULL) {
		_EXCEPTIONT("Invalid log buffer");
	}

	AnnounceSetOutputBuffer(param.fpLog);
	AnnounceOutputOnAllRanks();

	// Dereference pointers to operators
	_ASSERT(param.pvecThresholdOp != NULL);
	std::vector<ThresholdOp> & vecThresholdOp =
		*(param.pvecThresholdOp);

	// Dereference pointers to operators
	_ASSERT(param.pvecFilterOp != NULL);
	std::vector<FilterOp> & vecFilterOp =
		*(param.pvecFilterOp);

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

		// Get the longitude dimension
		NcDim * dimLon = vecFiles[0]->get_dim(param.strLongitudeName.c_str());
		if (dimLon == NULL) {
			_EXCEPTION1("Error accessing dimension \"%s\"", param.strLongitudeName.c_str());
		}

		// Get the longitude variable
		NcVar * varLon = vecFiles[0]->get_var(param.strLongitudeName.c_str());
		if (varLon == NULL) {
			_EXCEPTION1("Error accessing variable \"%s\"", param.strLongitudeName.c_str());
		}

		// Get the latitude dimension
		NcDim * dimLat = vecFiles[0]->get_dim(param.strLatitudeName.c_str());
		if (dimLat == NULL) {
			_EXCEPTION1("Error accessing dimension \"%s\"", param.strLatitudeName.c_str());
		}

		// Get the latitude variable
		NcVar * varLat = vecFiles[0]->get_var(param.strLatitudeName.c_str());
		if (varLat == NULL) {
			_EXCEPTION1("Error accessing variable \"%s\"", param.strLatitudeName.c_str());
		}

		nLat = dimLat->size();
		nLon = dimLon->size();

		DataArray1D<double> vecLat(nLat);
		varLat->get(vecLat, nLat);

		for (int j = 0; j < nLat; j++) {
			vecLat[j] *= M_PI / 180.0;
		}

		DataArray1D<double> vecLon(nLon);
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

	DataArray1D<double> dTime(nTime);

	if (varTime != NULL) {
		if (varTime->type() == ncDouble) {
			varTime->get(dTime, nTime);

		} else if (varTime->type() == ncFloat) {
			DataArray1D<float> dTimeFloat(nTime);

			varTime->get(dTimeFloat, nTime);
			for (int t = 0; t < nTime; t++) {
				dTime[t] = static_cast<double>(dTimeFloat[t]);
			}

		} else if (varTime->type() == ncInt) {
			DataArray1D<int> dTimeInt(nTime);

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

	// Create output variable
	CopyNcVarIfExists(ncInput, ncOutput, param.strLatitudeName, true);
	CopyNcVarIfExists(ncInput, ncOutput, param.strLongitudeName, true);

	NcDim * dim0 = NULL;
	NcDim * dim1 = NULL;

	NcVar * varTag = NULL;

	if (grid.m_nGridDim.size() == 1) {
		dim0 = ncOutput.get_dim("ncol");
		if (dim0 == NULL) {
			dim0 = ncOutput.add_dim("ncol", grid.GetSize());
		}
		if (dim0 == NULL) {
			_EXCEPTION1("Error creating dim \"ncol\" in file \"%s\"",
				strOutputFile.c_str());
		}
		if (dim0->size() != grid.GetSize()) {
			_EXCEPTION4("Dimension \"%s\" in file \"%s\" has inconsistent length (%i vs %i)",
				dim0->name(),
				strOutputFile.c_str(),
				dim0->size(),
				grid.GetSize());
		}

		if (dimTimeOut != NULL) {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dimTimeOut,
				dim0);

		} else {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dim0);
		}

	} else if (grid.m_nGridDim.size() == 2) {

		// Copy over latitude and longitude variable
		dim0 = ncOutput.get_dim(param.strLatitudeName.c_str());
		if (dim0 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				param.strLatitudeName.c_str());
		}
		dim1 = ncOutput.get_dim(param.strLongitudeName.c_str());
		if (dim1 == NULL) {
			_EXCEPTION1("Error copying variable \"%s\" to output file",
				param.strLongitudeName.c_str());
		}

		// Create output tag
		if (dimTimeOut != NULL) {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dimTimeOut,
				dim0,
				dim1);

		} else {
			varTag = ncOutput.add_var(
				param.strTagVar.c_str(),
				ncByte,
				dim0,
				dim1);
		}

	} else {
		_EXCEPTIONT("Invalid grid dimension -- value must be 1 or 2");
	}

	_ASSERT(varTag != NULL);

	AnnounceEndBlock("Done");

	// Tagged cell array
	DataArray1D<int> bTag(grid.GetSize());

	// Loop through all times
	for (int t = 0; t < nTime; t ++) {

		// Announce
		AnnounceStartBlock("Time %i", t);

		AnnounceStartBlock("Build tagged cell array");
		bTag.Zero();

		// Set all points within the specified latitude bounds to 1
		for (int i = 0; i < grid.GetSize(); i++) {
			if (fabs(grid.m_dLat[i]) < param.dMinAbsLat * M_PI / 180.0) {
				continue;
			}
			if (grid.m_dLat[i] < param.dMinLat * M_PI / 180.0) {
				continue;
			}
			if (grid.m_dLat[i] > param.dMaxLat * M_PI / 180.0) {
				continue;
			}

			bTag[i] = 1;
		}
		AnnounceEndBlock("Done");

		// Eliminate based on threshold commands
		AnnounceStartBlock("Apply threshold commands");
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			// Load the search variable data
			Variable & var = varreg.Get(vecThresholdOp[tc].m_varix);
			var.LoadGridData(varreg, vecFiles, grid, t);
			const DataArray1D<float> & dataState = var.GetData();

			// Loop through data
			for (int i = 0; i < grid.GetSize(); i++) {
				if (bTag[i] == 0) {
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
					bTag[i] = 0;
				}
			}
		}
		AnnounceEndBlock("Done");

		// Eliminate based on filter commands
		if (vecFilterOp.size() != 0) {
			AnnounceStartBlock("Apply filters");
			for (int fc = 0; fc < vecFilterOp.size(); fc++) {

				// Load the search variable data
				Variable & var = varreg.Get(vecFilterOp[fc].m_varix);
				var.LoadGridData(varreg, vecFiles, grid, t);
				const DataArray1D<float> & dataState = var.GetData();

				// Check the blobs against the filter
				ApplyFilters<float>(
					grid,
					dataState,
					vecFilterOp[fc].m_eOp,
					vecFilterOp[fc].m_dValue,
					vecFilterOp[fc].m_nCount,
					bTag);
			}
			AnnounceEndBlock("Done");
		}

		// Output tagged cell array
		AnnounceStartBlock("Writing results");
		if (dimTimeOut != NULL) {
			if (grid.m_nGridDim.size() == 1) {
				varTag->set_cur(t, 0);
				varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0]);

			} else if (grid.m_nGridDim.size() == 2) {
				varTag->set_cur(t, 0, 0);
				varTag->put(&(bTag[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

			} else {
				_EXCEPTION();
			}

			for (int oc = 0; oc < param.pvecOutputOp->size(); oc++) {
				if (grid.m_nGridDim.size() == 1) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dimTimeOut,
							dim0);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					var.LoadGridData(varreg, vecFiles, grid, t);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(t, 0);
					ncvar->put(&(dataState[0]), 1, grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dimTimeOut,
							dim0,
							dim1);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					var.LoadGridData(varreg, vecFiles, grid, t);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(t, 0, 0);
					ncvar->put(&(dataState[0]), 1, grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}

			}

		} else {
			if (grid.m_nGridDim.size() == 1) {
				varTag->set_cur((long)0);
				varTag->put(&(bTag[0]), grid.m_nGridDim[0]);
			} else if (grid.m_nGridDim.size() == 2) {
				varTag->set_cur(0, 0);
				varTag->put(&(bTag[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);
			} else {
				_EXCEPTION();
			}

			for (int oc = 0; oc < param.pvecOutputOp->size(); oc++) {
				if (grid.m_nGridDim.size() == 1) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dim0);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					var.LoadGridData(varreg, vecFiles, grid, t);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur((long)0);
					ncvar->put(&(dataState[0]), grid.m_nGridDim[0]);

				} else if (grid.m_nGridDim.size() == 2) {
					NcVar * ncvar =
						ncOutput.add_var(
							(*param.pvecOutputOp)[oc].m_strName.c_str(),
							ncFloat,
							dim0,
							dim1);

					Variable & var = varreg.Get((*param.pvecOutputOp)[oc].m_varix);
					var.LoadGridData(varreg, vecFiles, grid, t);
					const DataArray1D<float> & dataState = var.GetData();

					ncvar->set_cur(0, 0);
					ncvar->put(&(dataState[0]), grid.m_nGridDim[0], grid.m_nGridDim[1]);

				} else {
					_EXCEPTION();
				}
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
	// Parameters for DetectBlobs
	DetectBlobsParam dbparam;

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

	// Threshold commands
	std::string strThresholdCmd;

	// Filter commands
	std::string strFilterCmd;

	// Output commands
	std::string strOutputCmd;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strInputFileList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutput, "out", "");
		CommandLineString(strOutputFileList, "out_file_list", "");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strFilterCmd, "filtercmd", "", "[var,op,value,count]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,name;...]");
		CommandLineDouble(dbparam.dMinAbsLat, "minabslat", 0.0);
		CommandLineDouble(dbparam.dMinLat, "minlat", -90.0);
		CommandLineDouble(dbparam.dMaxLat, "maxlat", 90.0);
		CommandLineBool(dbparam.fRegional, "regional");
		CommandLineString(dbparam.strTagVar, "tagvar", "binary_tag");
		CommandLineString(dbparam.strLongitudeName, "lonname", "lon");
		CommandLineString(dbparam.strLatitudeName, "latname", "lat");
		CommandLineInt(dbparam.iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Set verbosity level
	AnnounceSetVerbosityLevel(dbparam.iVerbosityLevel);

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
	dbparam.pvecThresholdOp = &vecThresholdOp;

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

	// Parse the filter operator command string
	std::vector<FilterOp> vecFilterOp;
	dbparam.pvecFilterOp = &vecFilterOp;

	if (strFilterCmd != "") {
		AnnounceStartBlock("Parsing filter operations");

		int iLast = 0;
		for (int i = 0; i <= strFilterCmd.length(); i++) {

			if ((i == strFilterCmd.length()) ||
				(strFilterCmd[i] == ';') ||
				(strFilterCmd[i] == ':')
			) {
				std::string strSubStr =
					strFilterCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecFilterOp.size());
				vecFilterOp.resize(iNextOp + 1);
				vecFilterOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the threshold operator command string
	std::vector<BlobOutputOp> vecOutputOp;
	dbparam.pvecOutputOp = &vecOutputOp;

	if (strOutputCmd != "") {
		AnnounceStartBlock("Parsing output operations");

		int iLast = 0;
		for (int i = 0; i <= strOutputCmd.length(); i++) {

			if ((i == strOutputCmd.length()) ||
				(strOutputCmd[i] == ';') ||
				(strOutputCmd[i] == ':')
			) {
				std::string strSubStr =
					strOutputCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecOutputOp.size());
				vecOutputOp.resize(iNextOp + 1);
				vecOutputOp[iNextOp].Parse(varreg, strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

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
			dbparam.fpLog = stdout;

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
			dbparam.fpLog = fopen(strLogFile.c_str(), "w");
		}

		// Perform DetectBlobs
		DetectBlobs(
			f,
			vecInputFiles[f],
			strOutputFile,
			strConnectivity,
			varreg,
			dbparam);

		// Close the log file
		if (vecInputFiles.size() != 1) {
			fclose(dbparam.fpLog);
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

