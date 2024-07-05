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

/*
#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif
*/

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
#include <cfloat>

///////////////////////////////////////////////////////////////////////////////

static const int MaxHistogramGrids = 1000;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a histogram operator.
///	</summary>
class HistogramOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	HistogramOp() :
		m_varix(InvalidVariableIndex),
		m_dOffset(0.0),
		m_dBinWidth(0.0)
	{ }

public:
	///	<summary>
	///		Parse a nearbyblob operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Offset,
			ReadMode_BinWidth,
			ReadMode_Invalid
		} eReadMode = ReadMode_Offset;

		// Parse variable
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delimited
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in the offset
				if (eReadMode == ReadMode_Offset) {
					m_dOffset = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_BinWidth;

				// Read in the bin width
				} else if (eReadMode == ReadMode_BinWidth) {
					m_dBinWidth = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nExcess entries in --histogram op \"%s\""
							"\nRequired: \"<variable>,<offset>,<bin width>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in --histogram op \"%s\""
					"\nRequired: \"<variable>,<offset>,<bin width>]\"",
					strOp.c_str());
		}

		if (m_dBinWidth <= 0.0) {
			_EXCEPTION1("For --histogram, bin width (%2.6f) must be positive", m_dBinWidth);
		}

		// Output announcement
		Announce("Histogram of %s with offset %f and bin width %f",
			varreg.GetVariableString(m_varix).c_str(),
			m_dOffset,
			m_dBinWidth);
	}

public:
	///	<summary>
	///		Variable to use for the histogram.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Offset for the histogram (the leftmost edge of one of the bins).
	///	</summary>
	double m_dOffset;

	///	<summary>
	///		Bin width for the histogram.
	///	</summary>
	double m_dBinWidth;
};


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

	// Diagonal connectivity for RLL grids
	bool fDiagonalConnectivity;

	// Output grid type
	std::string strOutputGrid;

	// Output data file
	std::string strOutputData;

	// Start of time range to use
	std::string strTimeBegin;

	// End of time range to use
	std::string strTimeEnd;

	// List of variables to process
	std::string strVariables;

	// List of variables to output
	std::string strOutputVariables;

	// Data contains missing values
	bool fMissingData;

	// Scaling factor for max sample distance
	double dMaxDistanceScale;

	// List of operators
	std::string strOperators;

	// Histogram operators
	std::string strHistogramOp;

	// Remove the value from a nodefile column
	std::string strMinusNodefileColumn;

	// Grid spacing of output (Cartesian great-circle or radial distance)
	double dDeltaXDeg;

	// Resolution of output (Cartesian grid size or number of radial points)
	int nResolutionX;

	// Resolution of the output (azimuthal points)
	int nResolutionA;

	// Fixed longitude
	double dFixedLongitudeDeg;

	// Fixed latitude
	double dFixedLatitudeDeg;

	// Flip the composition in the y direction if the object is in the southern hemisphere.
	bool fFlipSouthernHemisphere;

	// Apply a time delta to the data
	std::string strApplyTimeDelta;

	// Maximum time difference
	std::string strMaxTimeDelta;

	// Output a snapshot series rather than the time-average
	bool fSnapshots;

	// Name of latitude dimension
	std::string strLatitudeName;

	// Name of longitude dimension
	std::string strLongitudeName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputNodeFile, "in_nodefile", "");
		//CommandLineString(strInputNodeFileList, "in_file_list", "");
		CommandLineStringD(strInputNodeFileType, "in_nodefile_type", "SN", "[DN|SN]");
		CommandLineString(strInputFormat, "in_fmt", "(auto)");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fDiagonalConnectivity, "diag_connect");
		CommandLineBool(fRegional, "regional");

		CommandLineStringD(strOutputGrid, "out_grid", "XY", "[XY|RAD|RLL]");
		CommandLineString(strOutputData, "out_data", "");

		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputVariables, "varout", "");
		CommandLineDouble(dMaxDistanceScale, "maxdistscale", 0.9);
		CommandLineBool(fMissingData, "missingdata");
		CommandLineBool(fSnapshots, "snapshots");
		CommandLineStringD(strOperators, "op", "mean", "[mean|min|max,...]");
		CommandLineStringD(strHistogramOp, "histogram", "", "[var,offset,binsize;...]");
		CommandLineString(strMinusNodefileColumn, "minus_nodefile_col", "");

		CommandLineDouble(dDeltaXDeg, "dx", 0.5);
		CommandLineInt(nResolutionX, "resx", 11);
		CommandLineInt(nResolutionA, "resa", 16);
		CommandLineDouble(dFixedLongitudeDeg, "fixlon", -999.);
		CommandLineDouble(dFixedLatitudeDeg, "fixlat", -999.);
		CommandLineBool(fFlipSouthernHemisphere, "flip_hemisphere");

		CommandLineString(strApplyTimeDelta, "apply_time_delta", "");
		CommandLineString(strMaxTimeDelta, "max_time_delta", "");

		CommandLineString(strLongitudeName, "lonname", "lon");
		CommandLineString(strLatitudeName, "latname", "lat");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registries
	VariableRegistry varreg;

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
	if ((strOperators.length() == 0) &&
		(strHistogramOp.length() == 0) &&
		(!fSnapshots)
	) {
		_EXCEPTIONT("At least one of --op, --histogram or --snapshots must be specified");
	}
	if (dFixedLongitudeDeg != -999.) {
		if (dFixedLatitudeDeg == -999.) {
			_EXCEPTIONT("--fixlon and --fixlat must be specified together");
		}
		if ((dFixedLongitudeDeg < -360.0) || (dFixedLongitudeDeg > 360.0)) {
			_EXCEPTION1("--fixlon %1.5f out of range [-360,+360]", dFixedLongitudeDeg);
		}
	}
	if (dFixedLatitudeDeg != -999.) {
		if (dFixedLongitudeDeg == -999.) {
			_EXCEPTIONT("--fixlon and --fixlat must be specified together");
		}
		if ((dFixedLatitudeDeg < -90.0) || (dFixedLongitudeDeg > 90.0)) {
			_EXCEPTION1("--fixlat %1.5f out of range [-90,+90]", dFixedLatitudeDeg);
		}
	}
	if (fSnapshots) {
		if ((dFixedLongitudeDeg != -999.) || (dFixedLatitudeDeg != -999.)) {
			_EXCEPTIONT("--fixlon and --fixlat cannot be combined with --snapshots");
		}
	}

	// Output grid options
	STLStringHelper::ToLower(strOutputGrid);
	if ((strOutputGrid != "xy") && (strOutputGrid != "rad") && (strOutputGrid != "rll")) {
		_EXCEPTIONT("Grid type of output (--out_grid) must be \"xy\", \"rad\" or \"rll\"");
	}
	if (dDeltaXDeg <= 0.0) {
		_EXCEPTIONT("Grid spacing of output (--dx) must be nonnegative");
	}
	if (nResolutionX < 1) {
		_EXCEPTIONT("Resolution of output (--resx) must be nonnegative");
	}
	if ((strOutputGrid == "rad") && (nResolutionA < 8)) {
		_EXCEPTIONT("Resolution of output (--resa) must be >= 8");
	}
	if ((strOutputGrid == "rll") && (dFixedLongitudeDeg == -999.)) {
		_EXCEPTIONT("Grid \"rll\" may only be used with fixed coordinate composites");
	}
	if (fFlipSouthernHemisphere && (strOutputGrid != "xy")) {
		_EXCEPTIONT("--flip_hemisphere may only be used with --out_grid \"xy\"");
	}

	// Input file type
	NodeFile::PathType iftype;
	if (strInputNodeFileType == "DN") {
		iftype = NodeFile::PathTypeDN;
	} else if (strInputNodeFileType == "SN") {
		iftype = NodeFile::PathTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_nodefile_type, expected \"SN\" or \"DN\"");
	}

	// Convert degrees to radians
	double dDeltaXRad = dDeltaXDeg * M_PI / 180.0;
	double dFixedLongitudeRad = -999.;
	double dFixedLatitudeRad = -999.;

	if (dFixedLongitudeDeg != -999.) {
		dFixedLongitudeRad = dFixedLongitudeDeg * M_PI / 180.0;
	}
	if (dFixedLatitudeDeg != -999.) {
		dFixedLatitudeRad = dFixedLatitudeDeg * M_PI / 180.0;
	}

	// NodeFile
	NodeFile nodefile;

	// Parse --in_fmt string
	ColumnDataHeader cdhInput;
	cdhInput.Parse(strInputFormat);

	int iLonColIx = (-1);
	int iLatColIx = (-1);

	// Parse --var argument
	std::vector<std::string> vecVariableNamesIn;
	VariableIndexVector vecVarIxIn;
	if (strVariables != "") {
		std::string strVariablesTemp = strVariables;
		for (;;) {
			VariableIndex varixIn;
			int iLast = varreg.FindOrRegisterSubStr(strVariablesTemp, &varixIn) + 1;

			if (fMissingData) {
				Variable & var = varreg.Get(varixIn);
				if (var.IsOp()) {
					_EXCEPTIONT("When using --missingdata flag, variable --var cannot be an operator (not implemented)");
				}
			}

			vecVarIxIn.push_back(varixIn);
			vecVariableNamesIn.push_back(varreg.GetVariableString(varixIn));

			if (iLast >= strVariablesTemp.length()) {
				break;
			}
			strVariablesTemp = strVariablesTemp.substr(iLast);
		}

		_ASSERT(vecVariableNamesIn.size() == vecVarIxIn.size());
	}
 
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

	// Parse --op argument
	bool fCompositeMean = false;
	bool fCompositeMin = false;
	bool fCompositeMax = false;

	if (strOperators != "") {
		int i = 0;
		int iLast = 0;
		for (;;) {
			if ((i == strOperators.length()) || (strOperators[i] == ',')) {
				std::string strOp = strOperators.substr(iLast, i - iLast);
				STLStringHelper::ToLower(strOp);
				if (strOp == "mean") {
					fCompositeMean = true;
				} else if (strOp == "min") {
					fCompositeMin = true;
				} else if (strOp == "max") {
					fCompositeMax = true;
				} else {
					_EXCEPTION1("Invalid --op argument \"%s\"", strOp.c_str());
				}
				iLast = i+1;

				if (i == strOperators.length()) {
					break;
				}
			}

			i++;
		}
	}

	// Parse --histogram argument
	const int NoHistogram = (-1);

	std::vector<int> iVarHistogramOpIx;
	iVarHistogramOpIx.resize(vecVarIxIn.size(), NoHistogram);

	std::vector<HistogramOp> vecHistogramOps;
	if (strHistogramOp != "") {
		int i = 0;
		int iLast = 0;
		for (;;) {
			if ((i == strHistogramOp.length()) || (strHistogramOp[i] == ';')) {
				HistogramOp opHistogram;
				opHistogram.Parse(varreg, strHistogramOp.substr(iLast, i - iLast));

				int v = 0;
				for (; v < vecVarIxIn.size(); v++) {
					if (opHistogram.m_varix == vecVarIxIn[v]) {
						break;
					}
				}
				if (v == vecVarIxIn.size()) {
					_EXCEPTION1("Variable \"%s\" must appear in --var argument",
						varreg.GetVariableString(opHistogram.m_varix).c_str());
				}
				if (iVarHistogramOpIx[v] != NoHistogram) {
					_EXCEPTION1("Variable \"%s\" can only appear once in --histogram",
						varreg.GetVariableString(opHistogram.m_varix).c_str());
				}

				iVarHistogramOpIx[v] = vecHistogramOps.size();
				vecHistogramOps.push_back(opHistogram);

				iLast = i + 1;
				if (i == strHistogramOp.length()) {
					break;
				}
			}

			i++;
		}
	}

	// Parse --apply_time_delta
	bool fAddTimeDelta = true;
	Time timeDelta;
	if (strApplyTimeDelta != "") {
		if (strApplyTimeDelta[0] == '-') {
			fAddTimeDelta = false;
			strApplyTimeDelta = strApplyTimeDelta.substr(1);
		}
		timeDelta.FromFormattedString(strApplyTimeDelta);
	}

	// Parse --max_time_delta
	Time timeMaxDelta(Time::CalendarNone);
	if (strMaxTimeDelta != "") {
		timeMaxDelta.FromFormattedString(strMaxTimeDelta);
		if (timeMaxDelta.GetTimeType() != Time::TypeDelta) {
			_EXCEPTIONT("--max_time_delta must be specified as a time difference ([#d][#h][#m][#s])");
		}
	}

	// Parse --minus_nodefile_col and get column indices for each output variable
	std::vector<int> vecMinusNodefileColIx;
	if (strMinusNodefileColumn != "") {
		if ((dFixedLatitudeRad != -999.) || (dFixedLongitudeRad != -999.)) {
			_EXCEPTIONT("Cannot combine --minus_nodefile_col with --fixlon and --fixlat");
		}
		std::vector<std::string> vecMinusNodefileCol;
		STLStringHelper::ParseVariableList(strMinusNodefileColumn, vecMinusNodefileCol);
		if ((vecMinusNodefileCol.size() != vecVariableNamesIn.size()) &&
		    (vecMinusNodefileCol.size() != 1)
		) {
			_EXCEPTIONT("--minus_nodefile_col must be a single name or a list of same length as --var");
		}
		for (std::string & strNodefileCol : vecMinusNodefileCol) {
			int ixCol = cdhInput.GetIndexFromString(strNodefileCol);
			if (ixCol == ColumnDataHeader::InvalidIndex) {
				_EXCEPTION1("Column name from --minus_nodefile_col \"%s\" must appear in --in_fmt",
					strNodefileCol.c_str());
			}
			vecMinusNodefileColIx.push_back(ixCol);
			if ((vecMinusNodefileCol.size() == 1) && (vecVariableNamesIn.size() != 1)) {
				for (int n = 0; n < vecVariableNamesIn.size()-1; n++) {
					vecMinusNodefileColIx.push_back(ixCol);
				}
			}
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
	if (autocurator.GetTimeCount() == 0) {
		_EXCEPTIONT("No time slices found among input file(s)");
	}
	AnnounceEndBlock("Done");

	// If --regional then also --missingdata
	if (fRegional) {
		fMissingData = true;
	}

	// Get the time units
	NcType eNcTimeType = autocurator.GetNcTimeType();
	std::string strNcTimeUnits = autocurator.GetNcTimeUnits();

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
			strLatitudeName,
			strLongitudeName,
			fRegional,
			fDiagonalConnectivity);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
		AnnounceEndBlock("Done");
	}

	// Maximum distance permitted at each grid node
	DataArray1D<double> dMaxDistance;
	if (fRegional) {
		dMaxDistance.Allocate(grid.m_vecConnectivity.size());
		for (int i = 0; i < grid.m_vecConnectivity.size(); i++) {
			for (int j = 0; j < grid.m_vecConnectivity[i].size(); j++) {
				double dDist =
					GreatCircleDistance_Rad(
						grid.m_dLon[i],
						grid.m_dLat[i],
						grid.m_dLon[grid.m_vecConnectivity[i][j]],
						grid.m_dLat[grid.m_vecConnectivity[i][j]]);

				if (dDist > dMaxDistance[i]) {
					dMaxDistance[i] = dDist;
				}
			}
		}
		for (int i = 0; i < dMaxDistance.GetRows(); i++) {
			dMaxDistance[i] *= dMaxDistanceScale;
		}
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

	if ((vecInputNodeFiles.size() > 1) && (fSnapshots)) {
		_EXCEPTIONT("Only one nodefile allowed with --snapshotseries");
	}

	// Open output file
	AnnounceStartBlock("Preparing output file");
	NcFile ncoutfile(strOutputData.c_str(), NcFile::Replace);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output datafile \"%s\"",
			strOutputData.c_str());
	}

	// Write dimensions
	long lOutputDimSize0 = (-1);
	long lOutputDimSize1 = (-1);

	if (strOutputGrid == "xy") {
		NcDim * dimX = ncoutfile.add_dim("x", nResolutionX);
		NcDim * dimY = ncoutfile.add_dim("y", nResolutionX);

		DataArray1D<double> dX(nResolutionX);
		for (int i = 0; i < nResolutionX; i++) {
			dX[i] = dDeltaXDeg * (
				static_cast<double>(i)
				- 0.5 * static_cast<double>(nResolutionX-1));
		}

		NcVar * varX = ncoutfile.add_var("x", ncDouble, dimX);
		varX->put(&(dX[0]), nResolutionX);
		varX->add_att("name", "stereographic x coordinate");
		varX->add_att("units", "degrees_east");

		NcVar * varY = ncoutfile.add_var("y", ncDouble, dimY);
		varY->put(&(dX[0]), nResolutionX);
		varY->add_att("name", "stereographic y coordinate");
		varY->add_att("units", "degrees_north");

		lOutputDimSize0 = nResolutionX;
		lOutputDimSize1 = nResolutionX;

	} else if (strOutputGrid == "rad") {
		NcDim * dimX = ncoutfile.add_dim("az", nResolutionA);
		NcDim * dimY = ncoutfile.add_dim("r", nResolutionX);

		DataArray1D<double> dAz(nResolutionA);
		for (int i = 0; i < nResolutionA; i++) {
			dAz[i] = 360.0 * static_cast<double>(i)
				/ static_cast<double>(nResolutionA);
		}

		DataArray1D<double> dR(nResolutionX);
		for (int i = 0; i < nResolutionX; i++) {
			dR[i] = dDeltaXDeg * (static_cast<double>(i) + 0.5);
		}

		NcVar * varAz = ncoutfile.add_var("az", ncDouble, dimX);
		varAz->put(&(dAz[0]), nResolutionA);
		varAz->add_att("name", "stereographic azimuth angle");
		varAz->add_att("units", "degrees");

		NcVar * varR = ncoutfile.add_var("r", ncDouble, dimY);
		varR->put(&(dR[0]), nResolutionX);
		varR->add_att("name", "stereographic great circle distance");
		varR->add_att("units", "degrees");

		lOutputDimSize0 = nResolutionX;
		lOutputDimSize1 = nResolutionA;

	} else if (strOutputGrid == "rll") {
		NcDim * dimX = ncoutfile.add_dim("lon", nResolutionX);
		NcDim * dimY = ncoutfile.add_dim("lat", nResolutionX);

		double dHalfWidthDeg = 0.5 * static_cast<double>(nResolutionX) * dDeltaXDeg;

		DataArray1D<double> dLonDeg(nResolutionX);
		for (int i = 0; i < nResolutionX; i++) {
			dLonDeg[i] = dFixedLongitudeDeg
				- dHalfWidthDeg
				+ dDeltaXDeg * (static_cast<double>(i) + 0.5);
		}

		DataArray1D<double> dLatDeg(nResolutionX);
		for (int j = 0; j < nResolutionX; j++) {
			dLatDeg[j] = dFixedLatitudeDeg
				- dHalfWidthDeg
				+ dDeltaXDeg * (static_cast<double>(j) + 0.5);
		}

		NcVar * varLon = ncoutfile.add_var("lon", ncDouble, dimX);
		varLon->put(&(dLonDeg[0]), nResolutionX);
		varLon->add_att("name", "longitude");
		varLon->add_att("units", "degrees_east");

		NcVar * varLat = ncoutfile.add_var("lat", ncDouble, dimY);
		varLat->put(&(dLatDeg[0]), nResolutionX);
		varLon->add_att("name", "latitude");
		varLat->add_att("units", "degrees_north");

		lOutputDimSize0 = nResolutionX;
		lOutputDimSize1 = nResolutionX;

	} else {
		_EXCEPTIONT("Invalid grid");
	}
	AnnounceEndBlock("Done");

	// Vector of output data
	int nDataInstances = 0;

	DataArray1D<float> dOutputDataSnapshot;
	
	std::vector<float> vecFillValueFloat;
	vecFillValueFloat.resize(vecVarIxIn.size(), -FLT_MAX);

	std::vector< DataArray1D<float> > vecOutputDataMean;
	vecOutputDataMean.resize(vecVarIxIn.size());

	std::vector< DataArray1D<float> > vecOutputDataMin;
	vecOutputDataMin.resize(vecVarIxIn.size());

	std::vector< DataArray1D<float> > vecOutputDataMax;
	vecOutputDataMax.resize(vecVarIxIn.size());

	// Vector of output data for histograms
	int nHistogramGrids = 0;
	typedef std::map<int, DataArray1D<int> *> HistogramMap;
	std::vector<HistogramMap> vecmapHistograms;

	// Vector of output NcDim * for each variable
	std::vector< std::vector<NcDim *> > vecOutputNcDim;
	vecOutputNcDim.resize(vecVarIxIn.size());

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

		for (int i = 0; i < nodefile.m_cdh.size(); i++) {
			if ((nodefile.m_cdh[i] == "lon") || (nodefile.m_cdh[i] == "longitude")) {
				iLonColIx = i;
			}
			if ((nodefile.m_cdh[i] == "lat") || (nodefile.m_cdh[i] == "latitude")) {
				iLatColIx = i;
			}
		}

		if ((iLonColIx == (-1)) || (iLatColIx == (-1))) {
			_EXCEPTIONT("--in_fmt or CSV header must contain \"lon\" and \"lat\" columns");
		}

		if (strApplyTimeDelta != "") {
			nodefile.ApplyTimeDelta(timeDelta, fAddTimeDelta);
		}

		nodefile.GenerateTimeToPathNodeMap();
		AnnounceEndBlock("Done");

		if (mapTimeToPathNode.size() == 0) {
			_EXCEPTIONT("No time slices found in file");
		}

		// Get number of paths and number of snapshots
		std::vector< DataArray1D<size_t> > vecPathCountNonMissing;
		if (fMissingData) {
			vecPathCountNonMissing.resize(vecVarIxIn.size());
		}

		size_t sPathCount = 0;
		size_t sSnapshotCount = 0;
		{
			const PathVector & pathvec = nodefile.GetPathVector();
			sPathCount = pathvec.size();
			for (size_t p = 0; p < sPathCount; p++) {
				sSnapshotCount += pathvec[p].size();
			}
			if (sSnapshotCount == 0) {
				_EXCEPTIONT("No pathnode entries found in nodefile; cannot build composite");
			}
			Announce("%lu path snapshots found", sSnapshotCount);
		}

		// Variables used for outputing snapshots
		NcDim * dimSnapshot = NULL;

		std::vector<NcVar *> vecvarSnapshots;
		vecvarSnapshots.resize(vecVarIxIn.size());

		std::vector<bool> vecSnapshotsFillValue;
		vecSnapshotsFillValue.resize(vecVarIxIn.size(), false);

		if (fSnapshots) {
			dimSnapshot = ncoutfile.add_dim("snapshot", sSnapshotCount);
			if (dimSnapshot == NULL) {
				_EXCEPTION1("Unable to create dimension \"snapshot\" in file %s",
					strOutputData.c_str());
			}

			// Create auxiliary variables
			DataArray1D<int> dataPathId(sSnapshotCount);
			DataArray1D<double> dataPathLonDeg(sSnapshotCount);
			DataArray1D<double> dataPathLatDeg(sSnapshotCount);
			DataArray1D<double> dataPathTimeDouble(sSnapshotCount);
			size_t ipathnode = 0;
			for (size_t p = 0; p < pathvec.size(); p++) {
				for (size_t n = 0; n < pathvec[p].size(); n++) {
					const PathNode & pn = pathvec[p][n];
					dataPathId[ipathnode] = p;
					dataPathLonDeg[ipathnode] = RadToDeg(grid.m_dLon[pn.m_gridix]);
					dataPathLatDeg[ipathnode] = RadToDeg(grid.m_dLat[pn.m_gridix]);

					dataPathTimeDouble[ipathnode] =
						pn.m_time.GetCFCompliantUnitsOffsetDouble(strNcTimeUnits);

					ipathnode++;
				}
			}

			// Output auxiliary variables
			NcVar * varSnapPathId = ncoutfile.add_var("snap_pathid", ncInt, dimSnapshot);
			if (varSnapPathId == NULL) {
				_EXCEPTION1("Unable to create variable \"snap_pathid\" in file %s",
					strOutputData.c_str());
			}
			varSnapPathId->put(&(dataPathId[0]), sSnapshotCount);

			NcVar * varSnapLon = ncoutfile.add_var("snap_lon", ncDouble, dimSnapshot);
			if (varSnapLon == NULL) {
				_EXCEPTION1("Unable to create variable \"snap_lon\" in file %s",
					strOutputData.c_str());
			}
			varSnapLon->put(&(dataPathLonDeg[0]), sSnapshotCount);
			varSnapLon->add_att("units", "degrees_east");

			NcVar * varSnapLat = ncoutfile.add_var("snap_lat", ncDouble, dimSnapshot);
			if (varSnapLat == NULL) {
				_EXCEPTION1("Unable to create variable \"snap_lat\" in file %s",
					strOutputData.c_str());
			}
			varSnapLat->put(&(dataPathLatDeg[0]), sSnapshotCount);
			varSnapLat->add_att("units", "degrees_north");

			NcVar * varSnapTime = ncoutfile.add_var("snap_time", eNcTimeType, dimSnapshot);
			if (varSnapTime == NULL) {
				_EXCEPTION1("Unable to create variable \"snap_time\" in file %s",
					strOutputData.c_str());
			}
			varSnapTime->add_att("units", strNcTimeUnits.c_str());

			if (eNcTimeType == ncInt) {
				DataArray1D<double> dataPathTimeInt(sSnapshotCount);
				for (size_t i = 0; i < sSnapshotCount; i++) {
					dataPathTimeInt[i] = static_cast<int>(dataPathTimeDouble[i]);
				}
				varSnapTime->put(&(dataPathTimeInt[0]), sSnapshotCount);

			} else if (eNcTimeType == ncFloat) {
				DataArray1D<double> dataPathTimeFloat(sSnapshotCount);
				for (size_t i = 0; i < sSnapshotCount; i++) {
					dataPathTimeFloat[i] = static_cast<float>(dataPathTimeDouble[i]);
				}
				varSnapTime->put(&(dataPathTimeFloat[0]), sSnapshotCount);

			} else if (eNcTimeType == ncDouble) {
				varSnapTime->put(&(dataPathTimeDouble[0]), sSnapshotCount);
			}
		}

		// Loop through all Times in the NodeFile
		for (auto iter = mapTimeToPathNode.begin(); iter != mapTimeToPathNode.end(); iter++) {

			// Generate a NcFileVector at this Time
			AnnounceStartBlock("Time %s", iter->first.ToString().c_str());

			NcFileVector vecncDataFiles;
			if (strMaxTimeDelta == "") {
				autocurator.FindFilesAtTime(
					iter->first,
					vecncDataFiles);
			} else {
				autocurator.FindFilesNearTime(
					iter->first,
					vecncDataFiles,
					timeMaxDelta);
			}

			if (vecncDataFiles.size() == 0) {
				_EXCEPTION1("Time (%s) does not exist in input data fileset",
					iter->first.ToString().c_str());
			}

			// If this is the first time through the loop load the aux dimension info
			// and generate the output data structures and file.
			if (!fOutputInitialized) {

				AnnounceStartBlock("Initializing processing queue and output variables");

				varreg.ClearProcessingQueue();

				_ASSERT(lOutputDimSize0 > 0);
				_ASSERT(lOutputDimSize1 > 0);

				// Loop through all variables
				for (int v = 0; v < vecVarIxIn.size(); v++) {

					// Get auxiliary dimension info and verify consistency
					DimInfoVector vecAuxDimInfo;

					varreg.AppendVariableToProcessingQueue(
						vecncDataFiles,
						grid,
						vecVarIxIn[v],
						&vecAuxDimInfo);

					size_t s2DVariables = vecAuxDimInfo.GetTotalSize();

					if ((s2DVariables != 1) && (vecHistogramOps.size() != 0)) {
						_EXCEPTIONT("Not implemented: Histograms only available for 2D variables");
					}

					// Generate output variables
					if (strOutputGrid == "xy") {
						vecAuxDimInfo.push_back(DimInfo("y", lOutputDimSize0));
						vecAuxDimInfo.push_back(DimInfo("x", lOutputDimSize1));

					} else if (strOutputGrid == "rad") {
						vecAuxDimInfo.push_back(DimInfo("r", lOutputDimSize0));
						vecAuxDimInfo.push_back(DimInfo("az", lOutputDimSize1));

					} else if (strOutputGrid == "rll") {
						vecAuxDimInfo.push_back(DimInfo("lat", lOutputDimSize0));
						vecAuxDimInfo.push_back(DimInfo("lon", lOutputDimSize1));
					}

					// Initialize data storage for output
					if (fSnapshots) {
						dOutputDataSnapshot.Allocate(lOutputDimSize0 * lOutputDimSize1);
					}
					if (fMissingData) {
						vecPathCountNonMissing[v].Allocate(vecAuxDimInfo.GetTotalSize());
					}
					if (fCompositeMean) {
						vecOutputDataMean[v].Allocate(vecAuxDimInfo.GetTotalSize());
					}
					if (fCompositeMin) {
						vecOutputDataMin[v].Allocate(vecAuxDimInfo.GetTotalSize());
						for (int i = 0; i < vecOutputDataMin[v].GetRows(); i++) {
							vecOutputDataMin[v][i] = FLT_MAX;
						}
					}
					if (fCompositeMax) {
						vecOutputDataMax[v].Allocate(vecAuxDimInfo.GetTotalSize());
						for (int i = 0; i < vecOutputDataMin[v].GetRows(); i++) {
							vecOutputDataMax[v][i] = -FLT_MAX;
						}
					}

					// Copy auxiliary dimension variables from input to output
					vecOutputNcDim[v].resize(vecAuxDimInfo.size());
					for (int d = 0; d < vecAuxDimInfo.size(); d++) {
						vecOutputNcDim[v][d] =
							ncoutfile.get_dim(vecAuxDimInfo[d].name.c_str());
						if (vecOutputNcDim[v][d] == NULL) {
							vecOutputNcDim[v][d] = ncoutfile.add_dim(
								vecAuxDimInfo[d].name.c_str(),
								vecAuxDimInfo[d].size);

							if (vecOutputNcDim[v][d] == NULL) {
								_EXCEPTION2("Error adding dimension \"%s\" (%li) to output file",
									vecAuxDimInfo[d].name.c_str(),
									vecAuxDimInfo[d].size);
							}

						} else {
							if (vecOutputNcDim[v][d]->size() != vecAuxDimInfo[d].size) {
								std::string strVarName =
									varreg.GetVariableString(vecVarIxIn[v]);

								_EXCEPTION4("Dimension size mismatch when initializing variable \"%s\": Expected dimension \"%s\" to have size \"%li\" (found \"%li\")",
									strVarName.c_str(),
									vecAuxDimInfo[d].name.c_str(),
									vecAuxDimInfo[d].size,
									vecOutputNcDim[v][d]->size());
							}
						}
					}

					// Initialize snapshot variables
					if (fSnapshots) {
						std::vector<NcDim *> vecSnapshotDims;
						vecSnapshotDims.push_back(dimSnapshot);
						for (int d = 0; d < vecOutputNcDim[v].size(); d++) {
							vecSnapshotDims.push_back(vecOutputNcDim[v][d]);
						}

						std::string strVarName = std::string("snap_");
						strVarName += vecVariableNamesOut[v];

						NcVar * varSnapshots =
							ncoutfile.add_var(
								strVarName.c_str(),
								ncFloat,
								vecSnapshotDims.size(),
								const_cast<const NcDim**>(&(vecSnapshotDims[0])));

						if (varSnapshots == NULL) {
							_EXCEPTION1("Unable to create variable \"%s\" in output file",
								strVarName.c_str());
						}

						vecvarSnapshots[v] = varSnapshots;
					}
				}

				// Allocate storage for histograms
				if (vecHistogramOps.size() != 0) {
					vecmapHistograms.resize(vecVarIxIn.size());
				}

				// Done
				fOutputInitialized = true;

				AnnounceEndBlock("Done");
			}

			// Generate the SimpleGrid for each node
			AnnounceStartBlock("Building composites");
			const PathNodeIndexVector & vecPathNodes = iter->second;

			// Loop through all Variables
			varreg.ResetProcessingQueue();
			while(varreg.AdvanceProcessingQueue()) {
			//for (int v = 0; v < vecVarIxIn.size(); v++) {
				size_t v = varreg.GetProcessingQueueVarPos();

				_ASSERT(v < vecVarIxIn.size());

				// Load the data for the search variable
				Variable & var = varreg.GetProcessingQueueVariable();
				var.LoadGridData(varreg, vecncDataFiles, grid);
				const DataArray1D<float> & dataState = var.GetData();
				_ASSERT(dataState.GetRows() == grid.GetSize());

				vecFillValueFloat[v] = var.GetFillValueFloat();

				if (fSnapshots && fMissingData && !vecSnapshotsFillValue[v]) {
					vecSnapshotsFillValue[v] = true;
					vecvarSnapshots[v]->add_att("_FillValue", vecFillValueFloat[v]);
				}

				Announce("%s", var.ToString(varreg).c_str());

				// Get the offset from the beginning of the array
				VariableAuxIndex auxix = varreg.GetProcessingQueueAuxIx();

				size_t sOffset = varreg.GetProcessingQueueOffset() * lOutputDimSize0 * lOutputDimSize1;

				//std::cout << sOffset << "/" << vecOutputDataMean[v].GetRows() << std::endl;

				/////////////////////////////////
				// PathNode centered composite
				for (int p = 0; p < vecPathNodes.size(); p++) {
					const Path & path = pathvec[vecPathNodes[p].first];
					const PathNode & pathnode = path[vecPathNodes[p].second];

					double dPathNodeLonRad =
						DegToRad(pathnode.GetColumnDataAsDouble(iLonColIx));
					double dPathNodeLatRad =
						DegToRad(pathnode.GetColumnDataAsDouble(iLatColIx));

					if (dPathNodeLatRad > 0.5 * M_PI + HighTolerance) {
						_EXCEPTION3("Latitude out of range (%2.14f) of pathnode (%i,%i)",
							RadToDeg(dPathNodeLatRad),
							vecPathNodes[p].first,
							vecPathNodes[p].second);
					}
/*
					int ixOrigin = static_cast<int>(pathnode.m_gridix);

					double dPathNodeLon = grid.m_dLon[ixOrigin];
					double dPathNodeLat = grid.m_dLat[ixOrigin];

					_ASSERT((ixOrigin >= 0) && (ixOrigin < grid.GetSize()));
*/
					// Fixed point composites
					if ((dFixedLatitudeRad != -999.) || (dFixedLongitudeRad != -999.)) {
						_ASSERT(dFixedLatitudeRad != -999.);
						_ASSERT(dFixedLongitudeRad != -999.);

						dPathNodeLonRad = dFixedLongitudeRad;
						dPathNodeLatRad = dFixedLatitudeRad;
					}

					// Obtain reference value from nodefile
					float dNodefileRefValue = 0.0;
					if (vecMinusNodefileColIx.size() != 0) {
						dNodefileRefValue =
							pathnode.GetColumnDataAsDouble(vecMinusNodefileColIx[v]);
					}

					// Generate the SimpleGrid for this pathnode
					SimpleGrid gridNode;
					if (strOutputGrid == "xy") {
						gridNode.GenerateRectilinearStereographic(
							dPathNodeLonRad,
							dPathNodeLatRad,
							nResolutionX,
							dDeltaXRad,
							fFlipSouthernHemisphere,
							false);

					} else if (strOutputGrid == "rad") {
						gridNode.GenerateRadialStereographic(
							dPathNodeLonRad,
							dPathNodeLatRad,
							nResolutionX,
							nResolutionA,
							dDeltaXRad,
							false);

					} else if (strOutputGrid == "rll") {
						double dHalfWidth =
							0.5 * static_cast<double>(nResolutionX) * dDeltaXRad;

						gridNode.GenerateRegionalLatitudeLongitude(
							dPathNodeLatRad - dHalfWidth,
							dPathNodeLatRad + dHalfWidth,
							dPathNodeLonRad - dHalfWidth,
							dPathNodeLonRad + dHalfWidth,
							nResolutionX,
							nResolutionX,
							fDiagonalConnectivity);
					}

					// Only calculate the mean
					if (fCompositeMean && !fCompositeMin && !fCompositeMax && !fSnapshots && !fMissingData && (vecHistogramOps.size() == 0)) {
						for (size_t i = 0; i < gridNode.GetSize(); i++) {
							int ixGridIn =
								grid.NearestNode(
									gridNode.m_dLon[i],
									gridNode.m_dLat[i]);
/*
							if (i == 0) {
								printf("%6f %1.5f %1.5f %1.5f %1.5f\n", dataState[ixGridIn],
									gridNode.m_dLon[i] * 180.0 / M_PI,
									gridNode.m_dLat[i] * 180.0 / M_PI,
									grid.m_dLon[ixGridIn] * 180.0 / M_PI,
									grid.m_dLat[ixGridIn] * 180.0 / M_PI);
							}
*/
							vecOutputDataMean[v][sOffset+i] +=
								dataState[ixGridIn];
						}

						if (dNodefileRefValue != 0.0) {
							for (size_t i = 0; i < gridNode.GetSize(); i++) {
								vecOutputDataMean[v][sOffset+i] -= dNodefileRefValue;
							}
						}

					// Calculate some subset of mean, min, max
					} else {
						for (size_t i = 0; i < gridNode.GetSize(); i++) {
							int ixGridIn =
								grid.NearestNode(
									gridNode.m_dLon[i],
									gridNode.m_dLat[i]);

							// Verify sample point isn't beyond sample range
							if (fRegional) {
								double dDist =
									GreatCircleDistance_Rad(
										gridNode.m_dLon[i],
										gridNode.m_dLat[i],
										grid.m_dLon[ixGridIn],
										grid.m_dLat[ixGridIn]);

								if (dDist > dMaxDistance[ixGridIn]) {
									if (fSnapshots) {
										dOutputDataSnapshot[i] = vecFillValueFloat[v];
									}

									continue;
								}
							}

							// State value minus reference (adjusted if state is not missing)
							float dStateValueMinusRef = dataState[ixGridIn];

							// Check for non-missing data
							if (fMissingData &&
							    (dataState[ixGridIn] != vecFillValueFloat[v]) &&
							    (!std::isnan(dataState[ixGridIn]))
							) {
								vecPathCountNonMissing[v][sOffset+i]++;
								dStateValueMinusRef -= dNodefileRefValue;
							} else if (!fMissingData) {
								dStateValueMinusRef -= dNodefileRefValue;
							}

							if (fSnapshots) {
								dOutputDataSnapshot[i] = dStateValueMinusRef;
							}

							if (fMissingData &&
							    ((dataState[ixGridIn] == vecFillValueFloat[v]) || std::isnan(dataState[ixGridIn]))
							) {
								continue;
							}

							if (fCompositeMean) {
								vecOutputDataMean[v][sOffset+i] += dStateValueMinusRef;
							}
							if (fCompositeMin) {
								if (dStateValueMinusRef < vecOutputDataMin[v][sOffset+i]) {
									vecOutputDataMin[v][sOffset+i] = dStateValueMinusRef;
								}
							}
							if (fCompositeMax) {
								if (dStateValueMinusRef > vecOutputDataMax[v][sOffset+i]) {
									vecOutputDataMax[v][sOffset+i] = dStateValueMinusRef;
								}
							}

							// Build histograms
							if (vecHistogramOps.size() != 0) {
								size_t sHistIx = varreg.GetProcessingQueueOffset();
								for (int hop = 0; hop < vecHistogramOps.size(); hop++) {
									int iBin =
										static_cast<int>(
											dStateValueMinusRef
											- vecHistogramOps[hop].m_dOffset
										) / vecHistogramOps[hop].m_dBinWidth;

									HistogramMap::iterator iter =
										vecmapHistograms[sHistIx].find(iBin);

									DataArray1D<int> * pdata = NULL;
									if (iter == vecmapHistograms[sHistIx].end()) {
										nHistogramGrids++;
										pdata = new DataArray1D<int>(gridNode.GetSize());
										vecmapHistograms[sHistIx].insert(
											HistogramMap::value_type(iBin, pdata));
									} else {
										pdata = iter->second;
									}
									if (nHistogramGrids > MaxHistogramGrids) {
										_EXCEPTION1("Sanity check failed: NodeFileCompose limits number of histogram grids to %i", MaxHistogramGrids);
									}
									(*pdata)[sOffset+i]++;
								}
							}
						}
					}

					// Output snapshots
					if (fSnapshots) {
						if ((pathnode.m_fileix < 0) ||
						    (pathnode.m_fileix >= sSnapshotCount)
						) {
							_EXCEPTION2("pathnode file index out of range (%lu/%lu)",
								pathnode.m_fileix, sSnapshotCount);
						}

						VariableAuxIndex lSize;
						lSize.push_back(1);
						for (size_t d = 0; d < auxix.size(); d++) {
							lSize.push_back(1);
						}
						lSize.push_back(lOutputDimSize0);
						lSize.push_back(lOutputDimSize1);

						VariableAuxIndex lPos = auxix;
						lPos.insert(lPos.begin(), pathnode.m_fileix);
						lPos.push_back(0);
						lPos.push_back(0);
/*
						for (size_t d = 0; d < lSize.size(); d++) {
							std::cout << lSize[d] << ",";
						}
						std::cout << std::endl;
						for (size_t d = 0; d < lPos.size(); d++) {
							std::cout << lPos[d] << ",";
						}
						std::cout << std::endl;
*/
						vecvarSnapshots[v]->set_cur(&(lPos[0]));

						vecvarSnapshots[v]->put(
							&(dOutputDataSnapshot[0]),
							&(lSize[0]));
					}

					// Fixed point composites only use the time, not the location
					if ((dFixedLatitudeRad != -999.) || (dFixedLongitudeRad != -999.)) {
						break;
					}
				}
			}

			AnnounceEndBlock(NULL);
			AnnounceEndBlock(NULL);
		}

		// Average all Variables
		if ((dFixedLatitudeRad != -999.) || (dFixedLongitudeRad != -999.)) {
			if (fMissingData) {
				for (int v = 0; v < vecOutputDataMean.size(); v++) {
					for (int i = 0; i < vecOutputDataMean[v].GetRows(); i++) {
						if (vecPathCountNonMissing[v][i] == 0) {
							vecOutputDataMean[v][i] = vecFillValueFloat[v];
						} else {
							vecOutputDataMean[v][i] /= 
								static_cast<float>(vecPathCountNonMissing[v][i]);
						}
					}
				}

			} else {
				for (int v = 0; v < vecOutputDataMean.size(); v++) {
					for (int i = 0; i < vecOutputDataMean[v].GetRows(); i++) {
						vecOutputDataMean[v][i] /=
							static_cast<float>(sPathCount);
					}
				}
			}
		} else {
			if (fMissingData) {
				for (int v = 0; v < vecOutputDataMean.size(); v++) {
					for (int i = 0; i < vecOutputDataMean[v].GetRows(); i++) {
						if (vecPathCountNonMissing[v][i] == 0) {
							vecOutputDataMean[v][i] = vecFillValueFloat[v];
						} else {
							vecOutputDataMean[v][i] /=
								static_cast<float>(vecPathCountNonMissing[v][i]);
						}
					}
				}
			
			} else {
				for (int v = 0; v < vecOutputDataMean.size(); v++) {
					for (int i = 0; i < vecOutputDataMean[v].GetRows(); i++) {
						vecOutputDataMean[v][i] /=
							static_cast<float>(sSnapshotCount);
					}
				}
			}
		}

		// Write output variables
		{
			if (!fOutputInitialized) {
				_EXCEPTIONT("Failed to initialize output (no time slices found)");
			}

			AnnounceStartBlock("Writing output");
			if (fCompositeMean) {
				_ASSERT(vecVarIxIn.size() == vecOutputDataMean.size());
			}
			if (fCompositeMin) {
				_ASSERT(vecVarIxIn.size() == vecOutputDataMin.size());
			}
			if (fCompositeMax) {
				_ASSERT(vecVarIxIn.size() == vecOutputDataMax.size());
			}

			long lOutputDimSize0;
			long lOutputDimSize1;

			if (strOutputGrid == "xy") {
				lOutputDimSize0 = nResolutionX;
				lOutputDimSize1 = nResolutionX;
			} else if (strOutputGrid == "rad") {
				lOutputDimSize0 = nResolutionX;
				lOutputDimSize1 = nResolutionA;
			} else if (strOutputGrid == "rll") {
				lOutputDimSize0 = nResolutionX;
				lOutputDimSize1 = nResolutionX;
			}

			//varreg.ResetProcessingQueue();
			//while(varreg.AdvanceProcessingQueue()) {
			for (int v = 0; v < vecVarIxIn.size(); v++) {

				Variable & var = varreg.Get(vecVarIxIn[v]);
				if (vecVariableNamesIn[v] == vecVariableNamesOut[v]) {
					AnnounceStartBlock("Writing \"%s\"",
						vecVariableNamesIn[v].c_str());

				} else {
					AnnounceStartBlock("Writing \"%s\" -> \"%s\"",
						vecVariableNamesIn[v].c_str(),
						vecVariableNamesOut[v].c_str());
				}

				std::vector<long> lSize(vecOutputNcDim[v].size());
				for (size_t d = 0; d < vecOutputNcDim[v].size(); d++) {
					lSize[d] = vecOutputNcDim[v][d]->size();
				}

				// Mean of composite
				if (fCompositeMean) {
					Announce("mean");
					std::string strVarNameMean = vecVariableNamesOut[v];

					NcVar * pvar =
						ncoutfile.add_var(
							strVarNameMean.c_str(),
							ncFloat,
							vecOutputNcDim[v].size(),
							const_cast<const NcDim**>(&(vecOutputNcDim[v][0])));

					if (pvar == NULL) {
						_EXCEPTION1("Unable to add variable \"%s\" to output file",
							strVarNameMean.c_str());
					}

					if (fMissingData) {
						pvar->add_att("_FillValue", vecFillValueFloat[v]);
					}

					pvar->put(&(vecOutputDataMean[v][0]), &(lSize[0]));
				}

				// Min of composite
				if (fCompositeMin) {
					Announce("min");
					std::string strVarNameMin = vecVariableNamesOut[v] + "_min";

					NcVar * pvar =
						ncoutfile.add_var(
							strVarNameMin.c_str(),
							ncFloat,
							vecOutputNcDim[v].size(),
							const_cast<const NcDim**>(&(vecOutputNcDim[v][0])));

					if (pvar == NULL) {
						_EXCEPTION1("Unable to add variable \"%s\" to output file",
							strVarNameMin.c_str());
					}

					if (fMissingData) {
						pvar->add_att("_FillValue", vecFillValueFloat[v]);
					}

					pvar->put(&(vecOutputDataMin[v][0]), &(lSize[0]));
				}

				// Max of composite
				if (fCompositeMax) {
					Announce("max");
					std::string strVarNameMax = vecVariableNamesOut[v] + "_max";

					NcVar * pvar =
						ncoutfile.add_var(
							strVarNameMax.c_str(),
							ncFloat,
							vecOutputNcDim[v].size(),
							const_cast<const NcDim**>(&(vecOutputNcDim[v][0])));

					if (pvar == NULL) {
						_EXCEPTION1("Unable to add variable \"%s\" to output file",
							strVarNameMax.c_str());
					}

					if (fMissingData) {
						pvar->add_att("_FillValue", vecFillValueFloat[v]);
					}

					pvar->put(&(vecOutputDataMax[v][0]), &(lSize[0]));
				}

				// Histogram of composite
				if (iVarHistogramOpIx[v] != NoHistogram) {
					const HistogramOp & opHist = vecHistogramOps[iVarHistogramOpIx[v]];
					std::string strVarNameHist = vecVariableNamesOut[v] + "_hist";

					int nBins = vecmapHistograms[v].size();

					Announce("histogram (%i bins)", nBins);

					char szBuffer[128];
					snprintf(szBuffer, 128, "hist%i", v);

					NcDim * dimHist = ncoutfile.add_dim(szBuffer, nBins);
					if (dimHist == NULL) {
						_EXCEPTION1("Unable to add dimension \"%s\" to output file",
							szBuffer);
					}

					NcVar * varHist = ncoutfile.add_var(szBuffer, ncDouble, dimHist);
					if (dimHist == NULL) {
						_EXCEPTION1("Unable to add variable \"%s\" to output file",
							szBuffer);
					}

					DataArray1D<double> dBins(nBins);
					auto iter = vecmapHistograms[v].begin();
					for (int b = 0; iter != vecmapHistograms[v].end(); iter++, b++) {
						dBins[b] =
							opHist.m_dOffset
							+ opHist.m_dBinWidth * (static_cast<double>(iter->first) + 0.5);
					}

					varHist->put(&(dBins[0]), (long)nBins);

					// Add the histogram dimension to the NcDim array for this variable
					std::vector<NcDim *> vecHistNcDims = vecOutputNcDim[v];
					_ASSERT(vecHistNcDims.size() >= grid.DimCount());

					int iHistDimPos = vecHistNcDims.size() - grid.DimCount();
					vecHistNcDims.insert(vecHistNcDims.begin() + iHistDimPos, dimHist);

					// Insert the new variable
					NcVar * pvar =
						ncoutfile.add_var(
							strVarNameHist.c_str(),
							ncInt,
							vecHistNcDims.size(),
							const_cast<const NcDim**>(&(vecHistNcDims[0])));

					if (pvar == NULL) {
						_EXCEPTION1("Unable to add variable \"%s\" to output file",
							strVarNameHist.c_str());
					}

					// Write the data
					iter = vecmapHistograms[v].begin();
					for (int b = 0; iter != vecmapHistograms[v].end(); iter++, b++) {
						DataArray1D<int> & dataHist = *(iter->second);
						int nNonZeros = 0;
						for (int i = 0; i < dataHist.GetRows(); i++) {
							if (dataHist[i] > 0) {
								nNonZeros++;
							}
						}

						pvar->set_cur(b, 0, 0);
						pvar->put(
							&(dataHist[0]),
							1,
							lOutputDimSize0,
							lOutputDimSize1);
					}
					AnnounceEndBlock(NULL);
				}

				AnnounceEndBlock("Done");
			}
			AnnounceEndBlock("Done");
		}

		// Clear the histogram data
		{
			for (int v = 0; v < vecmapHistograms.size(); v++) {
				for (auto iter = vecmapHistograms[v].begin();
				     iter != vecmapHistograms[v].end();
					 iter++
				) {
					delete iter->second;
				}
			}
		}

		// Done processing this nodefile
		AnnounceEndBlock("Done");

		// Done
		AnnounceBanner();
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


