///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculatePosthocOutput.cpp
///	\author  Paul Ullrich
///	\version August 14, 2018
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
#include "Variable.h"
#include "AutoCurator.h"
#include "DataMatrix.h"

#include "netcdfcpp.h"

#include <fstream>
#include <queue>
#include <set>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

struct VelocityVector {
	double dU;
	double dV;
};

///////////////////////////////////////////////////////////////////////////////

struct PathNode {
	Time m_time;
	int m_gridix;
	std::vector<std::string> m_vecData;
	double m_dVelocityLat;
	double m_dVelocityLon;
};

///////////////////////////////////////////////////////////////////////////////

struct Path {
	Time m_timeStart;
	std::vector<PathNode> m_vecPathNodes;
};

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<Path> PathVector;

///////////////////////////////////////////////////////////////////////////////

class RadialProfileCalculator {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RadialProfileCalculator() :
		m_varixU(0),
		m_varixV(0),
		m_nBins(0),
		m_dBinWidth(0.0)
	{ }

	///	<summary>
	///		Parse the command-line argument for parameters.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_UVar,
			ReadMode_VVar,
			ReadMode_Bins,
			ReadMode_BinWidth,
			ReadMode_Invalid
		} eReadMode = ReadMode_UVar;

		// Loop through string
		int iLast = 0;

		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in amount
				if (eReadMode == ReadMode_UVar) {
					std::string strUVar = strSubStr;

					Variable varU;
					varU.ParseFromString(varreg, strUVar);
					m_varixU = varreg.FindOrRegister(varU);

					iLast = i + 1;
					eReadMode = ReadMode_VVar;

				} else if (eReadMode == ReadMode_VVar) {
					std::string strVVar = strSubStr;

					Variable varV;
					varV.ParseFromString(varreg, strVVar);
					m_varixV = varreg.FindOrRegister(varV);

					iLast = i + 1;
					eReadMode = ReadMode_Bins;

				// Read in distance
				} else if (eReadMode == ReadMode_Bins) {
					m_nBins = atoi(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_BinWidth;

				// Read in min/max distance
				} else if (eReadMode == ReadMode_BinWidth) {
					m_dBinWidth = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nToo many entries in --radial_wind_profile \"%s\""
						"\nRequired: \"<u var>,<v var>,<bins>,<bin_width>\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in --radial_wind_profile \"%s\""
					"\nRequired: \"<u var>,<v var>,<bins>,<bin_width>\"",
					strOp.c_str());
		}
		if (m_nBins <= 0) {
			_EXCEPTIONT("\nNonpositive value of <bins> argument given");
		}
		if (m_dBinWidth <= 0.0) {
			_EXCEPTIONT("\nNonpositive value of <bin_width> argument given");
		}
		if (static_cast<double>(m_nBins) * m_dBinWidth > 180.0) {
			_EXCEPTIONT("\n<bins> x <bin_width> must be no larger than 180 (degrees)");
		}
	}

	///	<summary>
	///		Parse the command-line argument for parameters.
	///	</summary>
	void Apply(
		VariableRegistry & varreg,
		NcFileVector & vecFiles,
		SimpleGrid & grid,
		int iTime,
		const PathNode & pathnode,
		std::vector<double> & dProfile
	) {
		const int ix0 = pathnode.m_gridix;

		// Load the zonal wind data
		Variable & varU = varreg.Get(m_varixU);
		varU.LoadGridData(varreg, vecFiles, grid, iTime);
		const DataVector<float> & dataStateU = varU.GetData();

		// Load the meridional wind data
		Variable & varV = varreg.Get(m_varixV);
		varV.LoadGridData(varreg, vecFiles, grid, iTime);
		const DataVector<float> & dataStateV = varV.GetData();

		// Verify that dRadius is less than 180.0
		double dRadius = m_dBinWidth * static_cast<double>(m_nBins);

		if ((dRadius < 0.0) || (dRadius > 180.0)) {
			_EXCEPTIONT("Radius must be in the range [0.0, 180.0]");
		}

		// Check grid index
		if (ix0 >= grid.m_vecConnectivity.size()) {
			_EXCEPTION2("Grid index (%i) out of range (< %i)",
				ix0, static_cast<int>(grid.m_vecConnectivity.size()));
		}

		// Central lat/lon and Cartesian coord
		double dLon0 = grid.m_dLon[ix0];
		double dLat0 = grid.m_dLat[ix0];

		double dX0 = cos(dLon0) * cos(dLat0);
		double dY0 = sin(dLon0) * cos(dLat0);
		double dZ0 = sin(dLat0);

		// Allocate bins
		std::vector< std::vector<double> > dVelocities;
		dVelocities.resize(m_nBins);

		// Queue of nodes that remain to be visited
		std::queue<int> queueNodes;
		for (int n = 0; n < grid.m_vecConnectivity[ix0].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix0][n]);
		}

		// Set of nodes that have already been visited
		std::set<int> setNodesVisited;

		// Loop through all latlon elements
		while (queueNodes.size() != 0) {
			int ix = queueNodes.front();
			queueNodes.pop();

			if (setNodesVisited.find(ix) != setNodesVisited.end()) {
				continue;
			}

			setNodesVisited.insert(ix);

			// Don't perform calculation on central node
			if (ix == ix0) {
				continue;
			}

			// lat/lon and Cartesian coords of this point
			double dLat = grid.m_dLat[ix];
			double dLon = grid.m_dLon[ix];

			double dX = cos(dLon) * cos(dLat);
			double dY = sin(dLon) * cos(dLat);
			double dZ = sin(dLat);

			// Great circle distance to this element (in degrees)
			double dR =
				sin(dLat0) * sin(dLat)
				+ cos(dLat0) * cos(dLat) * cos(dLon - dLon0);

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

			if (dR >= dRadius) {
				continue;
			}

			// Velocities at this location
			double dUlon = dataStateU[ix];
			double dUlat = dataStateV[ix];

			// Cartesian velocities at this location
			double dUx = - sin(dLat) * cos(dLon) * dUlat - sin(dLon) * dUlon;
			double dUy = - sin(dLat) * sin(dLon) * dUlat + cos(dLon) * dUlon;
			double dUz = cos(dLat) * dUlat;

			double dUtot = sqrt(dUlon * dUlon + dUlat * dUlat);

			// Calculate local radial vector from central lat/lon
			// i.e. project \vec{X} - \vec{X}_0 to the surface of the
			//      sphere and normalize to unit length.
			double dRx = dX - dX0;
			double dRy = dY - dY0;
			double dRz = dZ - dZ0;

			double dDot = dRx * dX + dRy * dY + dRz * dZ;

			dRx -= dDot * dX;
			dRy -= dDot * dY;
			dRz -= dDot * dZ;

			double dMag = sqrt(dRx * dRx + dRy * dRy + dRz * dRz);

			dRx /= dMag;
			dRy /= dMag;
			dRz /= dMag;

			// Calculate local azimuthal velocity vector
			double dAx = dY * dRz - dZ * dRy;
			double dAy = dZ * dRx - dX * dRz;
			double dAz = dX * dRy - dY * dRx;
/*
			dDot = dAx * dAx + dAy * dAy + dAz * dAz;
			if (fabs(dDot - 1.0) > 1.0e-10) {
				std:: cout << dDot << std::endl;
				_EXCEPTIONT("Logic error");
			}
*/
			// Calculate radial velocity
			//double dUr = dUx * dRx + dUy * dRy + dUz * dRz;

			// Calculate azimuthal velocity
			double dUa = dUx * dAx + dUy * dAy + dUz * dAz;

			//printf("%1.5e %1.5e :: %1.5e %1.5e\n", dUlon, dUlat, dUr, dUa);

			// Determine bin
			int iBin = static_cast<int>(dR / m_dBinWidth);
			if (iBin >= m_nBins) {
				_EXCEPTIONT("Logic error");
			}

			dVelocities[iBin].push_back(dUa);

			if (iBin < m_nBins-1) {
				dVelocities[iBin+1].push_back(dUa);
			}

			// Add all neighbors of this point
			for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
				queueNodes.push(grid.m_vecConnectivity[ix][n]);
			}
		}

		// Construct radial profile
		dProfile.resize(m_nBins);

		dProfile[0] = 0.0;

		for (int i = 1; i < m_nBins; i++) {
			double dAvg = 0.0;
			if (dVelocities[i].size() != 0) {
				for (int j = 0; j < dVelocities[i].size(); j++) {
					dAvg += dVelocities[i][j];
				}
				dAvg /= static_cast<double>(dVelocities[i].size());
			}
			dProfile[i] = dAvg;
		}
	}

public:
	///	<summary>
	///		Variable to use for zonal velocity.
	///	</summary>
	VariableIndex m_varixU;

	///	<summary>
	///		Variable to use for meridional velocity.
	///	</summary>
	VariableIndex m_varixV;

	///	<summary>
	///		Number of bins to use in calculation.
	///	</summary>
	int m_nBins;

	///	<summary>
	///		Maximum radius (in degrees) to use in calculation.
	///	</summary>
	double m_dBinWidth;

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Polynomial cubic Hermite interpolation function using Fritsch-Carlson
///		method.
///	</summary>
///	<reference>
///		https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
///	</reference>
///	<param name="n">
///		Number of points to use in interpolant.
///	</param>
///	<param name="x">
///		Array of size n with independent coordinate values of interpolant.
///		The values of x must be monotone increasing or behavior is undefined.
///	</param>
///	<param name="y">
///		Array of size n with dependent coordinate values of interpolant.
///	</param>
///	<param name="c">
///		Output array of cubic polynomial coefficients of size (n-1) x 4.
///		The interpolation function is then defined by
///		  xdiff = x - x[i]
///		  f_i(x) = c[i][0]
///		         + xdiff * (c[i][1] + xdiff * (c[i][2] + xdiff * c[i][3]))
///	</param>
///	<param name="work">
///		Work array of size 3*n.
///	</param>
///	<returns>
///		(-1) if the value of n is less than 2.
///		(0) if the calculation completed successfully.
///	</returns>
int pchip_fc_coeff(
	int n,
	double * const x,
	double * const y,
	double * c,
	double * work
) {
	if (n < 2) {
		return (-1);
	}

	double * dx = work;
	double * dy = work+n;
	double * m = work+2*n;

	for (int i = 0; i < n-1; i++) {
		dx[i] = x[i+1] - x[i];
		dy[i] = y[i+1] - y[i];
		m[i] = dy[i] / dx[i];
		c[i*4] = y[i];
	}
	m[n-1] = dy[n-2] / dx[n-2];

	c[1] = m[0];
	for (int i = 0; i < n-2; i++) {
		double mi = m[i];
		double mn = m[i+1];
		if (mi * mn <= 0.0) {
			c[4*(i+1)+1] = 0.0;
		} else {
			double dc = dx[i] + dx[i+1];
			c[4*(i+1)+1] = 3.0 * dc / ((dc + dx[i+1]) / mi + (dc + dx[i]) / mn);
		}
	}

	for (int i = 0; i < n-2; i++) {
		double c1 = c[4*i+1];
		double invdx = 1.0/dx[i];
		double dc = c1 + c[4*(i+1)+1] - 2.0 * m[i];

		c[4*i+2] = (m[i] - c1 - dc) * invdx;
		c[4*i+3] = dc * invdx * invdx;
	}
	{
		int i = n-2;
		double c1 = c[4*i+1];
		double invdx = 1.0/dx[i];
		double dc = c1 + m[n-1] - 2.0 * m[i];

		c[4*i+2] = (m[i] - c1 - dc) * invdx;
		c[4*i+3] = dc * invdx * invdx;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

void CalculateStormVelocity(
	const SimpleGrid & grid,
	Path & path
) {
	const double MetersPerRadian = 111.325 * 1000.0 * 180.0 / M_PI;

	int nPathNodes = path.m_vecPathNodes.size();
	if (path.m_vecPathNodes.size() < 2) {
		_EXCEPTIONT("Path must contain at least two nodes");
	}

	DataVector<double> dX(nPathNodes);
	DataVector<double> dY(nPathNodes);
	DataMatrix<double> dC(nPathNodes-1, 4);
	DataVector<double> dWork(nPathNodes*3);

	// Calculate meridional velocity
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];

		dX[i] = pathnode.m_time - path.m_timeStart;

		if (pathnode.m_gridix >= grid.m_dLat.GetRows()) {
			_EXCEPTION2("Grid index of path (%i) out of range "
				"(only %i nodes in grid)",
				pathnode.m_gridix,
				grid.m_dLat.GetRows());
		}
		dY[i] = grid.m_dLat[pathnode.m_gridix];
	}

	const double dFinalDeltaT = dX[nPathNodes-1] - dX[nPathNodes-2];

	pchip_fc_coeff(nPathNodes, &(dX[0]), &(dY[0]), &(dC[0][0]), &(dWork[0]));

	for (int i = 0; i < nPathNodes-1; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];
		pathnode.m_dVelocityLat = dC[i][1];
	}
	path.m_vecPathNodes[nPathNodes-1].m_dVelocityLat =
		dC[nPathNodes-2][1]
		+ 2.0 * dC[nPathNodes-2][2] * dFinalDeltaT
		+ 3.0 * dC[nPathNodes-2][3] * dFinalDeltaT * dFinalDeltaT;

	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];
		pathnode.m_dVelocityLat *= MetersPerRadian;
	}

	// Calculate zonal velocity
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];
		dY[i] = grid.m_dLon[pathnode.m_gridix];
	}

	pchip_fc_coeff(nPathNodes, &(dX[0]), &(dY[0]), &(dC[0][0]), &(dWork[0]));

	for (int i = 0; i < nPathNodes-1; i++) {
		path.m_vecPathNodes[i].m_dVelocityLon = dC[i][1];
	}
	path.m_vecPathNodes[nPathNodes-1].m_dVelocityLon =
		dC[nPathNodes-2][1]
		+ 2.0 * dC[nPathNodes-2][2] * dFinalDeltaT
		+ 3.0 * dC[nPathNodes-2][3] * dFinalDeltaT * dFinalDeltaT;

	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];
		pathnode.m_dVelocityLon *=
			MetersPerRadian * cos(grid.m_dLat[pathnode.m_gridix]);
	}

/*
	for (int i = 0; i < nPathNodes-1; i++) {
		printf("%1.5e\n", path.m_vecPathNodes[i].m_dVelocityLat);
		printf("%1.5e %1.5e :: %1.5e %1.5e %1.5e %1.5e\n", dX[i], dY[i], dC[i][0], dC[i][1], dC[i][2], dC[i][3]); //path.m_vecPathNodes[i].m_dVelocityLat);
	}
	printf("%1.5e\n", path.m_vecPathNodes[nPathNodes-1].m_dVelocityLat);
	printf("%1.5e %1.5e\n", dX[nPathNodes-1], dY[nPathNodes-1]);
*/
/*
	for (int i = 0; i < nPathNodes; i++) {
		PathNode & pathnode = path.m_vecPathNodes[i];
		printf("%1.5e %1.5e\n", pathnode.m_dVelocityLat, pathnode.m_dVelocityLon);
	}
*/
	//_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

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

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();

try {

	// Input text file
	std::string strInputFile;

	// Input list of text files
	std::string strInputFileList;

	// Input file type
	std::string strInputFileType;

	// Input data file
	std::string strInputData;

	// Input list of data files
	std::string strInputDataList;

	// Connectivity file
	std::string strConnectivity;

	// Data is regional
	bool fRegional;

	// Output file
	std::string strOutputFile;

	// Append output to input file
	bool fOutputAppend;

	// Append trajectory velocities
	bool fAppendTrajectoryVelocity;

	// Variables for calculating radial wind profile
	std::string strRadialWindProfileVars;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_file", "");
		//CommandLineString(strInputFileList, "in_file_list", "");
		CommandLineStringD(strInputFileType, "in_file_type", "SN", "[DCU|SN]");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineBool(fRegional, "regional");

		CommandLineString(strOutputFile, "out_file", "");
		//CommandLineBool(fOutputAppend, "out_append");

		CommandLineBool(fAppendTrajectoryVelocity, "append_traj_vel");
		CommandLineStringD(strRadialWindProfileVars, "radial_wind_profile", "", "(U,V,bins,bin_width[,opts])");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Create Variable registry
	VariableRegistry varreg;

	// Create autocurator
	AutoCurator autocurator;

	// Check arguments
	if ((strInputData.length() == 0) && (strInputDataList.length() == 0)) {
		_EXCEPTIONT("No input data file (--in_data) or (--in_data_list)"
			" specified");
	}
	if ((strInputData.length() != 0) && (strInputDataList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_data) or (--in_data_list)"
			" may be specified");
	}
	if ((strInputFile.length() == 0) && (strInputFileList.length() == 0)) {
		_EXCEPTIONT("No input file (--in_file) or (--in_file_list)"
			" specified");
	}
	if ((strInputFile.length() != 0) && (strInputFileList.length() != 0)) {
		_EXCEPTIONT("Only one of (--in_file) or (--in_file_list)"
			" may be specified");
	}

	// Input file type
	enum InputFileType {
		InputFileTypeDCU,
		InputFileTypeSN
	} iftype;
	if (strInputFileType == "DCU") {
		iftype = InputFileTypeDCU;
	} else if (strInputFileType == "SN") {
		iftype = InputFileTypeSN;
	} else {
		_EXCEPTIONT("Invalid --in_file_type, expected \"SN\" or \"DCU\"");
	}

	// Radial wind profile can only be calculated with StitchNodes output
	if ((fAppendTrajectoryVelocity) && (iftype != InputFileTypeSN)) {
		_EXCEPTIONT("--append_traj_vel can only be used with --in_file_type SN");
	}

	if ((strRadialWindProfileVars != "") && (iftype != InputFileTypeSN)) {
		_EXCEPTIONT("--radial_wind_profile can only be used with --in_file_type SN");
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Curate input data
	AnnounceStartBlock("Autocurating data");
	if (strInputData.length() != 0) {
		autocurator.IndexFiles(strInputData);

	} else {
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
			autocurator.IndexFiles(strFileLine);
		}
	}
	AnnounceEndBlock("Done");

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

	// No connectivity file; check for latitude/longitude dimension
	} else {
		const std::vector<std::string> & vecFiles = autocurator.GetFilenames();

		if (vecFiles.size() < 1) {
			_EXCEPTIONT("No data files specified");
		}

		NcFile ncFile(vecFiles[0].c_str());
		if (!ncFile.is_valid()) {
			_EXCEPTION1("Unable to open NetCDF file \"%s\"", vecFiles[0].c_str());
		}

		grid.GenerateLatitudeLongitude(&ncFile, fRegional);

		if (grid.m_nGridDim.size() != 2) {
			_EXCEPTIONT("Logic error when generating connectivity");
		}
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

	// Parse radial input parameters
	RadialProfileCalculator radprofcalc;
	if (strRadialWindProfileVars != "") {
		radprofcalc.Parse(varreg, strRadialWindProfileVars);
	}

	// Loop through all files
	std::string strBuffer;

	std::vector<int> coord;
	coord.resize(grid.m_nGridDim.size());

	// Vector of Path information loaded from StitchNodes file
	PathVector vecPaths;

	// A map from Times to file lines, used for StitchNodes formatted output
	typedef std::map<Time, std::vector< std::pair<int, int> > > TimeToPathNodeMap;
	TimeToPathNodeMap mapTimeToPathNode;

	// Loop over all files
	for (int f = 0; f < vecInputFiles.size(); f++) {

		AnnounceStartBlock("Processing input (%s)", vecInputFiles[f].c_str());

		std::ifstream ifInput(vecInputFiles[f]);
		if (!ifInput.is_open()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFiles[f].c_str());
		}

		AnnounceStartBlock("Reading input file");
		int iLine = 1;
		for (;;) {

			int nCount = 0;
			Time time;

			// Read header lines
			{
				getline(ifInput, strBuffer);
				if (ifInput.eof()) {
					break;
				}

				std::istringstream iss(strBuffer);

				// DetectCyclonesUnstructured output
				if (iftype == InputFileTypeDCU) {
					int iYear;
					int iMonth;
					int iDay;
					int iHour;

					iss >> iYear;
					iss >> iMonth;
					iss >> iDay;
					iss >> nCount;
					iss >> iHour;

					if (iss.eof()) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					time = Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						autocurator.GetCalendarType());

				// StitchNodes output
				} else if (iftype == InputFileTypeSN) {
					std::string strStart;
					int iYear;
					int iMonth;
					int iDay;
					int iHour;

					iss >> strStart;
					iss >> nCount;
					iss >> iYear;
					iss >> iMonth;
					iss >> iDay;
					iss >> iHour;

					if (iss.bad()) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					vecPaths.resize(vecPaths.size() + 1);
					vecPaths[vecPaths.size()-1].m_timeStart =
						Time(
							iYear,
							iMonth-1,
							iDay-1,
							3600 * iHour,
							0,
							autocurator.GetCalendarType());

					vecPaths[vecPaths.size()-1].m_vecPathNodes.resize(nCount);
				}

				iLine++;
			}

			// Read contents under each header line
			for (int i = 0; i < nCount; i++) {

				getline(ifInput, strBuffer);
				if (ifInput.eof()) {
					break;
				}

				std::istringstream iss(strBuffer);

				for (int n = 0; n < grid.m_nGridDim.size(); n++) {
					iss >> coord[n];
					if ((coord[n] < 0) || (coord[n] >= grid.m_nGridDim[n])) {
						_EXCEPTION1("Coordinate index out of range on line %i", coord[n]);
					}
				}
				if (iss.eof()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, vecInputFiles[f].c_str());
				}

				std::string strBuf;
				std::vector<std::string> vecDelimitedOutput;
				for (;;) {
					iss >> strBuf;
					vecDelimitedOutput.push_back(strBuf);
					if (iss.eof()) {
						break;
					}
				}

				int nOutputSize = vecDelimitedOutput.size();

				// StitchNodes requires us to reorder the data 
				if (iftype == InputFileTypeSN) {
					PathNode & pathnode =
						vecPaths[vecPaths.size()-1].m_vecPathNodes[i];

					if (nOutputSize < 4) {
						_EXCEPTION2("Format error on line %i of \"%s\"",
							iLine, vecInputFiles[f].c_str());
					}

					int iYear = std::stoi(vecDelimitedOutput[nOutputSize-4]);
					int iMonth = std::stoi(vecDelimitedOutput[nOutputSize-3]);
					int iDay = std::stoi(vecDelimitedOutput[nOutputSize-2]);
					int iHour = std::stoi(vecDelimitedOutput[nOutputSize-1]);

					time = Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						autocurator.GetCalendarType());

					TimeToPathNodeMap::iterator iter =
						mapTimeToPathNode.find(time);
					if (iter == mapTimeToPathNode.end()) {
						iter = mapTimeToPathNode.insert(
							TimeToPathNodeMap::value_type(
								time, std::vector< std::pair<int,int> >())).first;
					}
					iter->second.push_back(
						std::pair<int,int>(
							static_cast<int>(vecPaths.size()-1),i));

					pathnode.m_time = time;

					if (coord.size() == 1) {
						pathnode.m_gridix = coord[0];
					} else if (coord.size() == 2) {
						pathnode.m_gridix = coord[0] + grid.m_nGridDim[1] * coord[1];
					} else {
						_EXCEPTIONT("Undefined behavior for SimpleGrid dimensionality > 2");
					}
					pathnode.m_vecData = vecDelimitedOutput;
				}

				iLine++;
			}

			// Calculate velocity at each point
			if (iftype == InputFileTypeSN) {

				// Calculate trajectory velocity
				Path & path = vecPaths[vecPaths.size()-1];
				CalculateStormVelocity(grid, path);

				// Append velocity to file
				if (fAppendTrajectoryVelocity) {
					for (int i = 0; i < path.m_vecPathNodes.size(); i++) {
						PathNode & pathnode = path.m_vecPathNodes[i];
						if (pathnode.m_vecData.size() < 4) {
							_EXCEPTIONT("Logic error");
						}

						char buf[100];
						sprintf(buf, "%3.6f", path.m_vecPathNodes[i].m_dVelocityLon);
						pathnode.m_vecData.insert(
							pathnode.m_vecData.end() - 4, buf);
						sprintf(buf, "%3.6f", path.m_vecPathNodes[i].m_dVelocityLat);
						pathnode.m_vecData.insert(
							pathnode.m_vecData.end() - 4, buf);
					}
				}
			}
		}
		AnnounceEndBlock("Done");

		// Perform the radial wind profile calculation
		if (strRadialWindProfileVars != "") {
			TimeToPathNodeMap::iterator iterPathNode =
				mapTimeToPathNode.begin();
			for (; iterPathNode != mapTimeToPathNode.end(); iterPathNode++) {
				const Time & time = iterPathNode->first;

				// Unload data from the VariableRegistry
				varreg.UnloadAllGridData();

				// Open NetCDF files with data at this time
				NcFileVector vecncDataFiles;
				int iTime;
				autocurator.Find(time, vecncDataFiles, iTime);
				if (vecncDataFiles.size() == 0) {
					_EXCEPTION1("Time (%s) does not exist in input data fileset",
						time.ToString().c_str());
				}

				// Loop through all PathNodes which need calculating at this time
				for (int i = 0; i < iterPathNode->second.size(); i++) {
					int iPath = iterPathNode->second[i].first;
					int iPathNode = iterPathNode->second[i].second;

					PathNode & pathnode =
						vecPaths[iPath].m_vecPathNodes[iPathNode];

					std::vector<double> dProfile;
					radprofcalc.Apply(
						varreg,
						vecncDataFiles,
						grid,
						iTime,
						pathnode,
						dProfile);

					// Write profile to string
					std::string strProfile = "\"";
					for (int i = 0; i < dProfile.size(); i++) {
						char buf[100];
						sprintf(buf, "%3.6f", dProfile[i]);
						strProfile += buf;
						if (i != dProfile.size()-1) {
							strProfile += ",";
						}
					}
					strProfile += "\"";

					pathnode.m_vecData.insert(
						pathnode.m_vecData.end() - 4, strProfile);
				}
			}
		}

		// Output
		if (strOutputFile != "") {
			AnnounceStartBlock("Writing output");
			FILE * fpOutput = fopen(strOutputFile.c_str(),"w");

			if (iftype == InputFileTypeSN) {
				for (int p = 0; p < vecPaths.size(); p++) {
					Path & path = vecPaths[p];
					fprintf(fpOutput, "start\t%i\t%i\t%i\t%i\t%i\n",
						static_cast<int>(path.m_vecPathNodes.size()),
						path.m_timeStart.GetYear(),
						path.m_timeStart.GetMonth(),
						path.m_timeStart.GetDay(),
						path.m_timeStart.GetSecond() / 3600);

					for (int i = 0; i < vecPaths[p].m_vecPathNodes.size(); i++) {
						PathNode & pathnode = path.m_vecPathNodes[i];
						fprintf(fpOutput, "\t%i", pathnode.m_gridix);

						for (int j = 0; j < pathnode.m_vecData.size(); j++) {
							fprintf(fpOutput, "\t%s", pathnode.m_vecData[j].c_str());
						}
						fprintf(fpOutput, "\n");
					}
				}

			} else {
				_EXCEPTIONT("Sorry, not yet implemented!");
			}
			AnnounceEndBlock("Done");
		}
		AnnounceEndBlock("Done");
	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif
}

