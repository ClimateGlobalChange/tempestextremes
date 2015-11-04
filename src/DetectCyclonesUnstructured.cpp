///////////////////////////////////////////////////////////////////////////////
///
///	\file    DetectCyclonesUnstructured.cpp
///	\author  Paul Ullrich
///	\version September 18, 2015
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
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
#include "DataVector.h"
#include "DataMatrix.h"
#include "TimeObj.h"

#include "kdtree.h"

#include "netcdfcpp.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A data structure describing the grid, including coordinates of
///		each data point and graph connectivity of elements.
///	</summary>
class SimpleGrid {

public:
	///	<summary>
	///		Generate the unstructured grid information for a
	///		longitude-latitude grid.
	///	</summary>
	void GenerateLatitudeLongitude(
		const DataVector<double> & vecLat,
		const DataVector<double> & vecLon,
		bool fRegional
	) {
		int nLat = vecLat.GetRows();
		int nLon = vecLon.GetRows();

		m_dLat.Initialize(nLon * nLat);
		m_dLon.Initialize(nLon * nLat);
		m_vecConnectivity.resize(nLon * nLat);

		m_nGridDim.resize(2);
		m_nGridDim[0] = nLat;
		m_nGridDim[1] = nLon;

		int ixs = 0;
		for (int j = 0; j < nLat; j++) {
		for (int i = 0; i < nLon; i++) {

			// Vectorize coordinates
			m_dLat[ixs] = vecLat[j];
			m_dLon[ixs] = vecLon[i];

			// Connectivity in each compass direction
			if (j != 0) {
				m_vecConnectivity[ixs].push_back((j-1) * nLon + i);
			}
			if (j != nLat-1) {
				m_vecConnectivity[ixs].push_back((j+1) * nLon + i);
			}

			if ((!fRegional) ||
			    ((i != 0) && (i != nLon-1))
			) {
				m_vecConnectivity[ixs].push_back(
					j * nLon + ((i + 1) % nLon));
				m_vecConnectivity[ixs].push_back(
					j * nLon + ((i + nLon - 1) % nLon));
			}

			ixs++;
		}
		}

	}

	///	<summary>
	///		Load the grid information from a file.
	///	</summary>
	void FromFile(
		const std::string & strGridInfoFile
	) {
		std::ifstream fsGrid(strGridInfoFile.c_str());
		if (!fsGrid.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				strGridInfoFile.c_str());
		}

		size_t nFaces;
		fsGrid >> nFaces;

		m_nGridDim.resize(1);
		m_nGridDim[0] = nFaces;

		m_dLon.Initialize(nFaces);
		m_dLat.Initialize(nFaces);
		m_vecConnectivity.resize(nFaces);

		for (size_t f = 0; f < nFaces; f++) {
			size_t sNeighbors;
			char cComma;
			fsGrid >> m_dLon[f];
			fsGrid >> cComma;
			fsGrid >> m_dLat[f];
			fsGrid >> cComma;
			fsGrid >> sNeighbors;
			fsGrid >> cComma;

			// Convert to radians
			m_dLon[f] *= M_PI / 180.0;
			m_dLat[f] *= M_PI / 180.0;

			// Load connectivity
			m_vecConnectivity[f].resize(sNeighbors);
			for (size_t n = 0; n < sNeighbors; n++) {
				fsGrid >> m_vecConnectivity[f][n];
				if (n != sNeighbors-1) {
					fsGrid >> cComma;
				}
				m_vecConnectivity[f][n]--;
			}
			if (fsGrid.eof()) {
				if (f != nFaces-1) {
					_EXCEPTIONT("Premature end of file");
				}
			}
		}
	}

	///	<summary>
	///		Get the size of the SimpleGrid (number of points).
	///	</summary>
	size_t GetSize() const {
		return (m_vecConnectivity.size());
	}

public:
	///	<summary>
	///		Longitude of each grid point.
	///	</summary>
	DataVector<double> m_dLon;

	///	<summary>
	///		Latitude of each grid point.
	///	</summary>
	DataVector<double> m_dLat;

	///	<summary>
	///		Connectivity of each grid point.
	///	</summary>
	std::vector< std::vector<int> > m_vecConnectivity;

	///	<summary>
	///		Grid dimensions.
	///	</summary>
	std::vector<size_t> m_nGridDim;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing a parsed variable name.
///	</summary>
class Variable {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Variable() :
		m_fOp(false),
		m_strName(),
		m_nSpecifiedDim(0)
	{ }

	///	<summary>
	///		Constructor from variable name.
	///	</summary>
	Variable(
		const std::string & strName
	) :
		m_fOp(false),
		m_strName(strName),
		m_nSpecifiedDim(0)
	{ }

public:
	///	<summary>
	///		Parse the variable information from a string.  Return the index
	///		of the first character after the variable information.
	///	</summary>
	int ParseFromString(
		const std::string & strIn
	) {
		m_fOp = false;
		m_strName = "";
		m_nSpecifiedDim = 0;
		m_varArg.clear();

		bool fDimMode = false;
		std::string strDim;

		if (strIn.length() >= 1) {
			if (strIn[0] == '_') {
				m_fOp = true;
			}
		}

		for (int n = 0; n <= strIn.length(); n++) {
			// Reading the variable name
			if (!fDimMode) {
				if (n == strIn.length()) {
					m_strName = strIn;
					return n;
				}
				if (strIn[n] == ',') {
					m_strName = strIn.substr(0, n);
					return n;
				}
				if (strIn[n] == '(') {
					m_strName = strIn.substr(0, n);
					fDimMode = true;
					continue;
				}
				if (strIn[n] == ')') {
					m_strName = strIn.substr(0, n);
					return n;
				}

			// Reading in dimensions
			} else if (!m_fOp) {
				if (m_nSpecifiedDim == 4) {
					_EXCEPTIONT("Only 4 dimensions / arguments may "
						"be specified");
				}
				if (n == strIn.length()) {
					_EXCEPTION1("Variable dimension list must be terminated"
						" with ): %s", strIn.c_str());
				}
				if (strIn[n] == ',') {
					if (strDim.length() == 0) {
						_EXCEPTIONT("Invalid dimension index in variable");
					}
					m_iDim[m_nSpecifiedDim] = atoi(strDim.c_str());
					m_nSpecifiedDim++;
					strDim = "";

				} else if (strIn[n] == ')') {
					if (strDim.length() == 0) {
						_EXCEPTIONT("Invalid dimension index in variable");
					}
					m_iDim[m_nSpecifiedDim] = atoi(strDim.c_str());
					m_nSpecifiedDim++;
					return (n+1);

				} else {
					strDim += strIn[n];
				}

			// Reading in arguments
			} else {
				if (m_nSpecifiedDim == 4) {
					_EXCEPTIONT("Only 4 dimensions / arguments may "
						"be specified");
				}
				if (n == strIn.length()) {
					_EXCEPTION1("Op argument list must be terminated"
						" with ): %s", strIn.c_str());
				}
				m_varArg.resize(m_nSpecifiedDim+1);
				n += m_varArg[m_nSpecifiedDim]
					.ParseFromString(strIn.substr(n));
				m_nSpecifiedDim++;

				if (strIn[n] == ')') {
					return (n+1);
				}
			}
		}

		_EXCEPTION1("Malformed variable string \"%s\"", strIn.c_str());
	}

	///	<summary>
	///		Get a string representation of this variable.
	///	</summary>
	std::string ToString() const {
		char szBuffer[20];
		std::string strOut = m_strName;
		if (m_nSpecifiedDim == 0) {
			return strOut;
		}
		strOut += "(";
		for (int d = 0; d < m_nSpecifiedDim; d++) {
			if (m_fOp) {
				strOut += m_varArg[d].ToString();
			} else {
				sprintf(szBuffer, "%i", m_iDim[d]);
				strOut += szBuffer;
			}
			if (d != m_nSpecifiedDim-1) {
				strOut += ",";
			} else {
				strOut += ")";
			}
		}
		return strOut;
	}

	///	<summary>
	///		Get this variable in the given NcFile.
	///	</summary>
	NcVar * GetFromNetCDF(
		NcFile & ncFile,
		int iTime = (-1)
	) const {
		if (m_fOp) {
			_EXCEPTION1("Cannot call GetFromNetCDF() on operator \"%s\"",
				m_strName.c_str());
		}

		NcVar * var = ncFile.get_var(m_strName.c_str());
		if (var == NULL) {
			_EXCEPTION1("Variable \"%s\" not found in input file",
				m_strName.c_str());
		}

		int nVarDims = var->num_dims();
		if ((nVarDims > 0) && (iTime != (-1))) {
			if (strcmp(var->get_dim(0)->name(), "time") != 0) {
				iTime = (-1);
			}
		}
		if (iTime != (-1)) {
			if ((nVarDims != m_nSpecifiedDim + 2) &&
				(nVarDims != m_nSpecifiedDim + 3)
			) {
				_EXCEPTION1("Dimension size inconsistency in \"%s\"",
					m_strName.c_str());
			}
		} else {
			if ((nVarDims != m_nSpecifiedDim + 1) &&
				(nVarDims != m_nSpecifiedDim + 2)
			) {
				_EXCEPTION1("Dimension size inconsistency in \"%s\"",
					m_strName.c_str());
			}
		}

		int nSetDims = 0;
		long iDim[7];
		memset(&(iDim[0]), 0, sizeof(int) * 7);

		if (iTime != (-1)) {
			iDim[0] = iTime;
			nSetDims++;
		}

		if (m_nSpecifiedDim != 0) {
			for (int i = 0; i < m_nSpecifiedDim; i++) {
				iDim[nSetDims] = m_iDim[i];
				nSetDims++;
			}
			iDim[nSetDims  ] = 0;
			iDim[nSetDims+1] = 0;
		}

		var->set_cur(&(iDim[0]));

		return var;
	}

	///	<summary>
	///		Load a data block from the NcFile.
	///	</summary>
	template <typename real>
	void LoadGridData(
		NcFile & ncFile,
		const SimpleGrid & grid,
		DataVector<real> & data,
		int iTime = (-1)
	) {
		// Get the data directly from a variable
		if (!m_fOp) {
			// Get pointer to variable
			NcVar * var = GetFromNetCDF(ncFile, iTime);
			if (var == NULL) {
				_EXCEPTION1("Variable \"%s\" not found in NetCDF file",
					m_strName.c_str());
			}

			// Check grid dimensions
			int nVarDims = var->num_dims();
			if (nVarDims < grid.m_nGridDim.size()) {
				_EXCEPTION1("Variable \"%s\" has insufficient dimensions",
					m_strName.c_str());
			}

			int nSize = 0;
			int nLat = 0;
			int nLon = 0;

			long nDataSize[7];
			for (int i = 0; i < 7; i++) {
				nDataSize[i] = 1;
			}

			// Latitude/longitude grid
			if (grid.m_nGridDim.size() == 2) {
				nLat = grid.m_nGridDim[0];
				nLon = grid.m_nGridDim[1];

				int nVarDimX0 = var->get_dim(nVarDims-2)->size();
				int nVarDimX1 = var->get_dim(nVarDims-1)->size();

				if (nVarDimX0 != nLat) {
					_EXCEPTION1("Dimension mismatch with variable"
						" \"%s\" on \"lat\"",
						m_strName.c_str());
				}
				if (nVarDimX1 != nLon) {
					_EXCEPTION1("Dimension mismatch with variable"
						" \"%s\" on \"lon\"",
						m_strName.c_str());
				}

				nDataSize[nVarDims-2] = nLat;
				nDataSize[nVarDims-1] = nLon;

			// Unstructured grid
			} else if (grid.m_nGridDim.size() == 1) {
				nSize = grid.m_nGridDim[0];

				int nVarDimX0 = var->get_dim(nVarDims-1)->size();

				if (nVarDimX0 != nSize) {
					_EXCEPTION1("Dimension mismatch with variable"
						" \"%s\" on \"ncol\"",
						m_strName.c_str());
				}

				nDataSize[nVarDims-1] = nSize;
			}

			// Load the data
			var->get(&(data[0]), &(nDataSize[0]));
			return;
		}

		// Evaluate the vector magnitude operator
		if (m_strName == "_VECMAG") {
			if (m_varArg.size() != 2) {
				_EXCEPTION1("_VECMAG expects two arguments: %i given",
					m_varArg.size());
			}
			DataVector<real> dataLeft(data.GetRows());
			DataVector<real> dataRight(data.GetRows());

			m_varArg[0].LoadGridData<real>(
				ncFile, grid, dataLeft, iTime);
			m_varArg[1].LoadGridData<real>(
				ncFile, grid, dataRight, iTime);

			for (int i = 0; i < data.GetRows(); i++) {
				data[i] =
					sqrt(dataLeft[i] * dataLeft[i]
						+ dataRight[i] * dataRight[i]);
			}

		// Evaluate the minus operator
		} else if (m_strName == "_DIFF") {
			if (m_varArg.size() != 2) {
				_EXCEPTION1("_VECMAG expects two arguments: %i given",
					m_varArg.size());
			}
			DataVector<real> dataLeft(data.GetRows());
			DataVector<real> dataRight(data.GetRows());

			m_varArg[0].LoadGridData<real>(
				ncFile, grid, dataLeft, iTime);
			m_varArg[1].LoadGridData<real>(
				ncFile, grid, dataRight, iTime);

			for (int i = 0; i < data.GetRows(); i++) {
				data[i] = dataLeft[i] - dataRight[i];
			}

		} else {
			_EXCEPTION1("Unexpected operator \"%s\"", m_strName.c_str());
		}
	}

public:
	///	<summary>
	///		Flag indicating this is an operator.
	///	</summary>
	bool m_fOp;

	///	<summary>
	///		Variable name.
	///	</summary>
	std::string m_strName;

	///	<summary>
	///		Number of dimensions specified.
	///	</summary>
	int m_nSpecifiedDim;

	///	<summary>
	///		Specified dimension values.
	///	</summary>
	int m_iDim[4];

	///	<summary>
	///		Specified operator arguments.
	///	</summary>
	std::vector<Variable> m_varArg;

};

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
		int iLast = m_var.ParseFromString(strOp) + 1;

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
		std::string strDescription = m_var.ToString();
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
	Variable m_var;

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
///		A class describing a general closed contour operation.
///	</summary>
class ClosedContourOp {

public:
	///	<summary>
	///		Parse a closed contour operation string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Amount,
			ReadMode_Distance,
			ReadMode_MinMaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Amount;

		// Get variable information
		int iLast = m_var.ParseFromString(strOp) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in amount
				if (eReadMode == ReadMode_Amount) {
					m_dDeltaAmount = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in distance
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_MinMaxDist;

				// Read in min/max distance
				} else if (eReadMode == ReadMode_MinMaxDist) {
					m_dMinMaxDist = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nToo many entries in closed contour op \"%s\""
						"\nRequired: \"<name>,<amount>,<distance>,"
						"<minmaxdist>\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in closed contour op \"%s\""
					"\nRequired: \"<name>,<amount>,<distance>,"
					"<minmaxdist>\"",
					strOp.c_str());
		}

		// Output announcement
		char szBuffer[128];

		std::string strDescription;
		strDescription += m_var.ToString();

		if (m_dDeltaAmount == 0.0) {
			_EXCEPTIONT("For closed contour op, delta amount must be non-zero");
		}
		if (m_dDistance <= 0.0) {
			_EXCEPTIONT("For closed contour op, distance must be positive");
		}
		if (m_dMinMaxDist < 0.0) {
			_EXCEPTIONT("For closed contour op, min/max dist must be nonnegative");
		}

		if (m_dDeltaAmount < 0.0) {
			Announce("%s decreases by %f over %f degrees"
				   " (max search %f deg)",
				m_var.ToString().c_str(),
				-m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);

		} else {
			Announce("%s increases by %f over %f degrees"
					" (min search %f deg)",
				m_var.ToString().c_str(),
				m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);
		}
	}

public:
	///	<summary>
	///		Variable to use for closed contour op.
	///	</summary>
	Variable m_var;

	///	<summary>
	///		Threshold amount.  If positive this represents a minimum
	///		increase.  If negative this represents a minimum decrease.
	///	</summary>
	double m_dDeltaAmount;

	///	<summary>
	///		Threshold distance.
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		Distance to search for min or max.
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class storing an output operator.
///	</summary>
class OutputOp {

public:
	///	<summary>
	///		Possible operations.
	///	</summary>
	enum Operation {
		Max,
		Min,
		Avg,
		MaxDist,
		MinDist
	};

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Op,
			ReadMode_Distance,
			ReadMode_Invalid
		} eReadMode = ReadMode_Op;

		// Get variable information
		int iLast = m_var.ParseFromString(strOp) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in operation
				if (eReadMode == ReadMode_Op) {
					if (strSubStr == "max") {
						m_eOp = Max;
					} else if (strSubStr == "min") {
						m_eOp = Min;
					} else if (strSubStr == "avg") {
						m_eOp = Avg;
					} else if (strSubStr == "maxdist") {
						m_eOp = MaxDist;
					} else if (strSubStr == "mindist") {
						m_eOp = MinDist;
					} else {
						_EXCEPTION1("Output invalid operation \"%s\"",
							strSubStr.c_str());
					}

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
							"\nRequired: \"<name>,<operation>,<distance>\"",
							strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in threshold op \"%s\""
					"\nRequired: \"<name>,<operation>,<distance>\"",
					strOp.c_str());
		}

		if (m_dDistance < 0.0) {
			_EXCEPTIONT("For output op, distance must be nonnegative");
		}

		// Output announcement
		std::string strDescription;

		if (m_eOp == Max) {
			strDescription += "Maximum of ";
		} else if (m_eOp == Min) {
			strDescription += "Minimum of ";
		} else if (m_eOp == Avg) {
			strDescription += "Average of ";
		} else if (m_eOp == MaxDist) {
			strDescription += "Distance to maximum of ";
		} else if (m_eOp == MinDist) {
			strDescription += "Distance to minimum of ";
		}

		char szBuffer[128];

		sprintf(szBuffer, "%s", m_var.ToString().c_str());
		strDescription += szBuffer;

		sprintf(szBuffer, " within %f degrees", m_dDistance);
		strDescription += szBuffer;

		Announce("%s", strDescription.c_str());
	}

public:
	///	<summary>
	///		Variable to use for output.
	///	</summary>
	Variable m_var;

	///	<summary>
	///		Operation.
	///	</summary>
	Operation m_eOp;

	///	<summary>
	///		Distance to use when applying operation.
	///	</summary>
	double m_dDistance;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all minima in the given DataMatrix.
///	</summary>
template <typename real>
void FindAllLocalMinima(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	std::set<int> & setMinima
) {
	int sFaces = grid.m_vecConnectivity.size();
	for (int f = 0; f < sFaces; f++) {
		
		bool fMinimum = true;

		real dValue = data[f];
		int sNeighbors = grid.m_vecConnectivity[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[grid.m_vecConnectivity[f][n]] < dValue) {
				fMinimum = false;
				break;
			}
		}

		if (fMinimum) {
			setMinima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the locations of all maxima in the given DataMatrix.
///	</summary>
template <typename real>
void FindAllLocalMaxima(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	std::set<int> & setMaxima
) {
	int sFaces = grid.m_vecConnectivity.size();
	for (int f = 0; f < sFaces; f++) {
		
		bool fMaximum = true;

		real dValue = data[f];
		int sNeighbors = grid.m_vecConnectivity[f].size();
		for (int n = 0; n < sNeighbors; n++) {
			if (data[grid.m_vecConnectivity[f][n]] > dValue) {
				fMaximum = false;
				break;
			}
		}

		if (fMaximum) {
			setMaxima.insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the minimum/maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalMinMax(
	const SimpleGrid & grid,
	bool fMinimum,
	const DataVector<real> & data,
	int ix0,
	double dMaxDist,
	int & ixExtremum,
	real & dMaxValue,
	float & dRMax
) {
	// Verify that dMaxDist is less than 180.0
	if (dMaxDist > 180.0) {
		_EXCEPTIONT("MaxDist must be less than 180.0");
	}

	// Initialize the maximum to the central location
	ixExtremum = ix0;
	dMaxValue = data[ix0];
	dRMax = 0.0;

	// Queue of nodes that remain to be visited
	std::queue<int> queueNodes;
	queueNodes.push(ixExtremum);

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

		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		// Great circle distance to this element
		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		if (fMinimum) {
			if (data[ix] < dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}

		} else {
			if (data[ix] > dMaxValue) {
				ixExtremum = ix;
				dMaxValue = data[ix];
				dRMax = dR;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a pair of Date values.
///	</summary>
void ParseDate(
	int nDate,
	int nDateSec,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	nDateYear  = nDate / 10000;
	nDateMonth = (nDate % 10000) / 100;
	nDateDay   = (nDate % 100);
	nDateHour  = nDateSec / 3600;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Parse a time value.
///	</summary>
void ParseTimeDouble(
	const std::string & strTimeUnits,
	const std::string & strTimeCalendar,
	double dTime,
	int & nDateYear,
	int & nDateMonth,
	int & nDateDay,
	int & nDateHour
) {
	// Get calendar type
	Time::CalendarType cal;
	if ((strTimeCalendar.length() >= 6) &&
		(strncmp(strTimeCalendar.c_str(), "noleap", 6) == 0)
	) {
		cal = Time::CalendarNoLeap;

	} else if (
		(strTimeCalendar.length() >= 8) &&
		(strncmp(strTimeCalendar.c_str(), "standard", 8) == 0)
	) {
		cal = Time::CalendarStandard;

	} else {
		_EXCEPTION1("Unknown calendar type \"%s\"", strTimeCalendar.c_str());
	}
/*
	Time time(Time::CalendarStandard);
	time.FromFormattedString("1800-01-01 00:00:00");
	printf("%1.15e %i\n", 3600.0 * 1577832.0, (int)(3600.0 * 1577832.0));
	time.AddHours(1577832);

	Announce("Time (YMDS): %i %i %i %i",
			time.GetYear(),
			time.GetMonth(),
			time.GetDay(),
			time.GetSecond());

	_EXCEPTION();
*/
	// Time format is "days since ..."
	if ((strTimeUnits.length() >= 11) &&
	    (strncmp(strTimeUnits.c_str(), "days since ", 11) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(11);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		int nDays = static_cast<int>(dTime);
		time.AddDays(nDays);

		int nSeconds = static_cast<int>(fmod(dTime, 1.0) * 86400.0);
		time.AddSeconds(nSeconds);

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());


		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	// Time format is "hours since ..."
	} else if (
	    (strTimeUnits.length() >= 12) &&
	    (strncmp(strTimeUnits.c_str(), "hours since ", 12) == 0)
	) {
		std::string strSubStr = strTimeUnits.substr(12);
		Time time(cal);
		time.FromFormattedString(strSubStr);

		time.AddHours(static_cast<int>(dTime));

		Announce("Time (YMDS): %i %i %i %i",
				time.GetYear(),
				time.GetMonth(),
				time.GetDay(),
				time.GetSecond());

		nDateYear = time.GetYear();
		nDateMonth = time.GetMonth();
		nDateDay = time.GetDay();
		nDateHour = time.GetSecond() / 3600;

		//printf("%s\n", strSubStr.c_str());

	} else {
		_EXCEPTIONT("Unknown \"time::units\" format");
	}
	//_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if the given field has a closed contour about this point.
///	</summary>
template <typename real>
bool HasClosedContour(
	const SimpleGrid & grid,
	const DataVector<real> & dataState,
	const int ix0,
	double dDeltaAmt,
	double dDeltaDist,
	double dMinMaxDist
) {
	// Verify arguments
	if (dDeltaAmt == 0.0) {
		_EXCEPTIONT("Closed contour amount must be non-zero");
	}
	if (dDeltaDist <= 0.0) {
		_EXCEPTIONT("Closed contour distance must be positive");
	}

	// Find min/max near point
	int ixOrigin;

	if (dMinMaxDist == 0.0) {
		ixOrigin = ix0;

	// Find a local minimum / maximum
	} else {
		real dValue;
		float dR;

		FindLocalMinMax<real>(
			grid,
			(dDeltaAmt > 0.0),
			dataState,
			ix0,
			dMinMaxDist,
			ixOrigin,
			dValue,
			dR);
	}

	//printf("%lu %lu : %lu %lu : %1.5f %1.5f\n", ix0 % grid.m_nGridDim[1], ix0 / grid.m_nGridDim[1], ixOrigin % grid.m_nGridDim[1], ixOrigin / grid.m_nGridDim[1], dataState[ix0], dataState[ixOrigin]);

	// Set of visited nodes
	std::set<int> setNodesVisited;

	// Set of nodes to visit
	std::queue<int> queueToVisit;
	queueToVisit.push(ixOrigin);

	// Reference value
	real dRefValue = dataState[ixOrigin];

	const double dLat0 = grid.m_dLat[ixOrigin];
	const double dLon0 = grid.m_dLon[ixOrigin];

	Announce(2, "Checking (%lu) : (%1.5f %1.5f)",
		ixOrigin, dLat0, dLon0);


	// Build up nodes
	while (queueToVisit.size() != 0) {
		int ix = queueToVisit.front();
		queueToVisit.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		// Great circle distance to this element
		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		Announce(2, "-- (%lu) : (%1.5f %1.5f) : dx %1.5f",
			ix, dLatThis, dLonThis, dR);

		//printf("-- (%lu %lu) %1.5f %1.5f\n", ix % grid.m_nGridDim[1], ix / grid.m_nGridDim[1], dR, dataState[ix] - dRefValue);

		// Check great circle distance
		if (dR > dDeltaDist) {
			Announce(2, "Failed criteria; returning");
			AnnounceEndBlock(2, NULL);
			return false;
		}

		// Verify sufficient increase in value
		if (dDeltaAmt > 0.0) {
			if (dataState[ix] - dRefValue >= dDeltaAmt) {
				continue;
			}

		// Verify sufficient decrease in value
		} else {
			if (dRefValue - dataState[ix] >= -dDeltaAmt) {
				continue;
			}
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueToVisit.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	// Report success with criteria
	Announce(2, "Passed criteria; returning");
	AnnounceEndBlock(2, NULL);
	return true;
}

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

		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

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
///		Find the maximum value of a field near the given point.
///	</summary>
///	<param name="dMaxDist">
///		Maximum distance from the initial point in degrees.
///	</param>
template <typename real>
void FindLocalAverage(
	const SimpleGrid & grid,
	const DataVector<real> & data,
	int ix0,
	double dMaxDist,
	float & dAverage
) {
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

	// Number of points
	float dSum = 0.0;
	int nCount = 0;

	// Loop through all latlon elements
	while (queueNodes.size() != 0) {
		int ix = queueNodes.front();
		queueNodes.pop();

		if (setNodesVisited.find(ix) != setNodesVisited.end()) {
			continue;
		}

		setNodesVisited.insert(ix);

		double dLatThis = grid.m_dLat[ix];
		double dLonThis = grid.m_dLon[ix];

		// Great circle distance to this element
		double dR = 180.0 / M_PI * acos(sin(dLat0) * sin(dLatThis)
				+ cos(dLat0) * cos(dLatThis) * cos(dLonThis - dLon0));

		if (dR > dMaxDist) {
			continue;
		}

		// Check for new local extremum
		dSum += data[ix];
		nCount++;

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	dAverage = dSum / static_cast<float>(nCount);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::verbose_nonfatal);

try {
	// Gravitational constant
	const float ParamGravity = 9.80616;

	// Radius of the Earth
	const float ParamEarthRadius = 6.37122e6;

	// Input dat file
	std::string strInputFile;

	// Connectivity file
	std::string strConnectivity;

	// Output file
	std::string strOutputFile;

	// Variable to search for the minimum
	std::string strSearchByMin;

	// Variable to search for the maximum
	std::string strSearchByMax;

	// Maximum latitude for detection
	double dMaxLatitude;

	// Minimum latitude for detection
	double dMinLatitude;

	// File containing information on topography
	std::string strTopoFile;

	// Maximum topographic height for a detection
	double dMaxTopoHeight;

	// Merge distance
	double dMergeDist;

	// Closed contour commands
	std::string strClosedContourCmd;

	// Closed contour commands
	std::string strNoClosedContourCmd;

	// Threshold commands
	std::string strThresholdCmd;

	// Output commands
	std::string strOutputCmd;

	// Time stride
	int nTimeStride;

	// Regional (do not wrap longitudinal boundaries)
	bool fRegional;

	// Output header
	bool fOutputHeader;

	// Verbosity level
	int iVerbosityLevel;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in_data", "");
		CommandLineString(strConnectivity, "in_connect", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineStringD(strSearchByMin, "searchbymin", "", "(default PSL)");
		CommandLineString(strSearchByMax, "searchbymax", "");
		CommandLineDoubleD(dMaxLatitude, "maxlat", 0.0, "(degrees)");
		CommandLineDoubleD(dMinLatitude, "minlat", 0.0, "(degrees)");
		CommandLineString(strTopoFile, "topofile", "");
		CommandLineDoubleD(dMaxTopoHeight, "maxtopoht", 0.0, "(m)");
		CommandLineDoubleD(dMergeDist, "mergedist", 0.0, "(degrees)");
		CommandLineStringD(strClosedContourCmd, "closedcontourcmd", "", "[var,dist,delta,minmaxdist;...]");
		CommandLineStringD(strNoClosedContourCmd, "noclosedcontourcmd", "", "[var,dist,delta,minmaxdist;...]");
		CommandLineStringD(strThresholdCmd, "thresholdcmd", "", "[var,op,value,dist;...]");
		CommandLineStringD(strOutputCmd, "outputcmd", "", "[var,op,dist;...]");
		CommandLineInt(nTimeStride, "timestride", 1);
		CommandLineBool(fRegional, "regional");
		CommandLineBool(fOutputHeader, "out_header");
		CommandLineInt(iVerbosityLevel, "verbosity", 0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Set verbosity level
	AnnounceSetVerbosityLevel(iVerbosityLevel);

	// Check input
	if (strInputFile.length() == 0) {
		_EXCEPTIONT("No input data file (--in_data) specified");
	}

	// Check output
	if (strOutputFile.length() == 0) {
		_EXCEPTIONT("No output file (--out) specified");
	}

	// Only one of search by min or search by max should be specified
	if ((strSearchByMin == "") && (strSearchByMax == "")) {
		strSearchByMin = "PSL";
	}
	if ((strSearchByMin != "") && (strSearchByMax != "")) {
		_EXCEPTIONT("Only one of --searchbymin or --searchbymax can"
			" be specified");
	}

	bool fSearchByMinima = false;
	Variable varSearchBy;
	if (strSearchByMin != "") {
		varSearchBy.ParseFromString(strSearchByMin);
		fSearchByMinima = true;
	}
	if (strSearchByMax != "") {
		varSearchBy.ParseFromString(strSearchByMax);
		fSearchByMinima = false;
	}

	// Parse the closed contour command string
	std::vector<ClosedContourOp> vecClosedContourOp;

	if (strClosedContourCmd != "") {
		AnnounceStartBlock("Parsing closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strClosedContourCmd.length(); i++) {

			if ((i == strClosedContourCmd.length()) ||
				(strClosedContourCmd[i] == ';') ||
				(strClosedContourCmd[i] == ':')
			) {
				std::string strSubStr =
					strClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecClosedContourOp.size());
				vecClosedContourOp.resize(iNextOp + 1);
				vecClosedContourOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the no closed contour command string
	std::vector<ClosedContourOp> vecNoClosedContourOp;

	if (strNoClosedContourCmd != "") {
		AnnounceStartBlock("Parsing no closed contour operations");

		int iLast = 0;
		for (int i = 0; i <= strNoClosedContourCmd.length(); i++) {

			if ((i == strNoClosedContourCmd.length()) ||
				(strNoClosedContourCmd[i] == ';') ||
				(strNoClosedContourCmd[i] == ':')
			) {
				std::string strSubStr =
					strNoClosedContourCmd.substr(iLast, i - iLast);
			
				int iNextOp = (int)(vecNoClosedContourOp.size());
				vecNoClosedContourOp.resize(iNextOp + 1);
				vecNoClosedContourOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the threshold operator command string
	std::vector<ThresholdOp> vecThresholdOp;

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
				vecThresholdOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Parse the output operator command string
	std::vector<OutputOp> vecOutputOp;

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
				vecOutputOp[iNextOp].Parse(strSubStr);

				iLast = i + 1;
			}
		}

		AnnounceEndBlock("Done");
	}

	// Check minimum latitude and maximum latitude
	if (dMaxLatitude < 0.0) {
		_EXCEPTIONT("--maxlat must be nonnegative");
	}
	if (dMinLatitude < 0.0) {
		_EXCEPTIONT("--minlat must be nonnegative");
	}

	dMaxLatitude *= M_PI / 180.0;
	dMinLatitude *= M_PI / 180.0;

	// Load the netcdf file
	NcFile ncInput(strInputFile.c_str());
	if (!ncInput.is_valid()) {
		_EXCEPTION1("Cannot open input file \"%s\"", strInputFile.c_str());
	}

	// Define the SimpleGrid
	SimpleGrid grid;

	// Dimensions
	int nSize = 0;
	int nLon = 0;
	int nLat = 0;

	// Check for connectivity file
	if (strConnectivity != "") {
		grid.FromFile(strConnectivity);

		nSize = grid.GetSize();

	// No connectivity file; check for latitude/longitude dimension
	} else {

		NcDim * dimLat = ncInput.get_dim("lat");
		if (dimLat == NULL) {
			_EXCEPTIONT("No dimension \"lat\" found in input file");
		}

		NcDim * dimLon = ncInput.get_dim("lon");
		if (dimLon == NULL) {
			_EXCEPTIONT("No dimension \"lon\" found in input file");
		}

		NcVar * varLat = ncInput.get_var("lat");
		if (varLat == NULL) {
			_EXCEPTIONT("No variable \"lat\" found in input file");
		}

		NcVar * varLon = ncInput.get_var("lon");
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
		grid.GenerateLatitudeLongitude(vecLat, vecLon, fRegional);
	}

	// Get time dimension
	NcDim * dimTime = ncInput.get_dim("time");
	if (dimTime == NULL) {
		_EXCEPTIONT("No dimension \"time\" found in input file");
	}

	NcVar * varTime = ncInput.get_var("time");
	if (varTime == NULL) {
		_EXCEPTIONT("No variable \"time\" found in input file");
	}

	int nTime = dimTime->size();

	DataVector<double> dTime;
	dTime.Initialize(nTime);

	varTime->get(dTime, nTime);

	// Search variable data
	DataVector<float> dataSearch(grid.GetSize());

	// Load topography data (if requested)
	NcVar * varPHIS = NULL;

	DataVector<float> dataPHIS(grid.GetSize());

	if (strTopoFile != "") {
		NcFile ncTopo(strTopoFile.c_str());
		if (!ncTopo.is_valid()) {
			_EXCEPTION1("Unable to open file \"%s\"", strTopoFile.c_str());
		}

		Variable varPHIS("PHIS");
		varPHIS.LoadGridData<float>(ncTopo, grid, dataPHIS);
	}

	if ((strTopoFile == "") && (dMaxTopoHeight != 0.0)) {
		_EXCEPTIONT("No topography file specified; required for --maxtopoht");
	}

	// Open output file
	FILE * fpOutput = fopen(strOutputFile.c_str(), "w");
	if (fpOutput == NULL) {
		_EXCEPTION1("Could not open output file \"%s\"",
			strOutputFile.c_str());
	}

	if (fOutputHeader) {
		fprintf(fpOutput, "#day\tmonth\tyear\tcount\thour\n");
		fprintf(fpOutput, "#\t#\ti\tj\tpsl_lon\tpsl_lat\twind_max\tr_wind_max\tpsl_min\n");
	}

	// Loop through all times
	for (int t = 0; t < nTime; t += nTimeStride) {
	//for (int t = 0; t < 1; t++) {

		char szStartBlock[128];
		sprintf(szStartBlock, "Time %i", t);
		AnnounceStartBlock(szStartBlock);

		// Load the data for the search variable
		varSearchBy.LoadGridData<float>(ncInput, grid, dataSearch, t);

		// Tag all minima
		std::set<int> setCandidates;

		if (strSearchByMin != "") {
			FindAllLocalMinima<float>(grid, dataSearch, setCandidates);
		} else {
			FindAllLocalMaxima<float>(grid, dataSearch, setCandidates);
		}

		// Total number of candidates
		int nTotalCandidates = setCandidates.size();

		int nRejectedLatitude = 0;
		int nRejectedTopography = 0;
		int nRejectedMerge = 0;

		DataVector<int> vecRejectedClosedContour(vecClosedContourOp.size());
		DataVector<int> vecRejectedNoClosedContour(vecNoClosedContourOp.size());
		DataVector<int> vecRejectedThreshold(vecThresholdOp.size());

		// Eliminate based on maximum latitude
		if (dMaxLatitude != 0.0) {
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];

				if (fabs(dLat) <= dMaxLatitude) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedLatitude++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on minimum latitude
		if (dMinLatitude != 0.0) {
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];

				if (fabs(dLat) >= dMinLatitude) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedLatitude++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on maximum topographic height
		if (dMaxTopoHeight != 0.0) {
			std::set<int> setNewCandidates;

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				float dPHIS = dataPHIS[*iterCandidate];

				double dTopoHeight =
					static_cast<double>(dPHIS / ParamGravity);

				if (dTopoHeight <= dMaxTopoHeight) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					nRejectedTopography++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on merge distance
		if (dMergeDist != 0.0) {
			std::set<int> setNewCandidates;

			// Calculate spherical distance
			double dSphDist = 2.0 * sin(0.5 * dMergeDist / 180.0 * M_PI);

			// Create a new KD Tree containing all nodes
			kdtree * kdMerge = kd_create(3);

			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];
				double dLon = grid.m_dLon[*iterCandidate];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				kd_insert3(kdMerge, dX, dY, dZ, (void*)(&(*iterCandidate)));
			}

			// Loop through all candidates find set of nearest neighbors
			iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {
				double dLat = grid.m_dLat[*iterCandidate];
				double dLon = grid.m_dLon[*iterCandidate];

				double dX = cos(dLon) * cos(dLat);
				double dY = sin(dLon) * cos(dLat);
				double dZ = sin(dLat);

				// Find all neighbors within dSphDist
				kdres * kdresMerge =
					kd_nearest_range3(kdMerge, dX, dY, dZ, dSphDist);

				// Number of neighbors
				int nNeighbors = kd_res_size(kdresMerge);
				if (nNeighbors == 0) {
					setNewCandidates.insert(*iterCandidate);

				} else {
					double dValue =
						static_cast<double>(dataSearch[*iterCandidate]);

					bool fExtrema = true;
					for (;;) {
						int * ppr = (int *)(kd_res_item_data(kdresMerge));

						if (fSearchByMinima) {
							if (static_cast<double>(dataSearch[*ppr]) < dValue) {
								fExtrema = false;
								break;
							}

						} else {
							if (static_cast<double>(dataSearch[*ppr]) > dValue) {
								fExtrema = false;
								break;
							}
						}

						int iHasMore = kd_res_next(kdresMerge);
						if (!iHasMore) {
							break;
						}
					}

					if (fExtrema) {
						setNewCandidates.insert(*iterCandidate);
					} else {
						nRejectedMerge++;
					}
				}

				kd_res_free(kdresMerge);
			}

			// Destroy the KD Tree
			kd_free(kdMerge);

			// Update set of pressure minima
			setCandidates = setNewCandidates;
		}

		// Eliminate based on thresholds
		for (int tc = 0; tc < vecThresholdOp.size(); tc++) {

			std::set<int> setNewCandidates;

			// Load the search variable data
			DataVector<float> dataState(grid.GetSize());
			
			vecThresholdOp[tc].m_var.LoadGridData<float>(
				ncInput, grid, dataState, t);

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if the threshold is satisfied
				bool fSatisfiesThreshold =
					SatisfiesThreshold<float>(
						grid,
						dataState,
						*iterCandidate,
						vecThresholdOp[tc].m_eOp,
						vecThresholdOp[tc].m_dValue,
						vecThresholdOp[tc].m_dDistance
					);

				// If not rejected, add to new pressure minima array
				if (fSatisfiesThreshold) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					vecRejectedThreshold[tc]++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on closed contours
		for (int ccc = 0; ccc < vecClosedContourOp.size(); ccc++) {
			std::set<int> setNewCandidates;

			// Load the search variable data
			DataVector<float> dataState(grid.GetSize());
			
			vecClosedContourOp[ccc].m_var.LoadGridData<float>(
				ncInput, grid, dataState, t);

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if a closed contour is present
				bool fHasClosedContour =
					HasClosedContour<float>(
						grid,
						dataState,
						*iterCandidate,
						vecClosedContourOp[ccc].m_dDeltaAmount,
						vecClosedContourOp[ccc].m_dDistance,
						vecClosedContourOp[ccc].m_dMinMaxDist
					);

				// If not rejected, add to new pressure minima array
				if (fHasClosedContour) {
					setNewCandidates.insert(*iterCandidate);
				} else {
					vecRejectedClosedContour[ccc]++;
				}
			}

			setCandidates = setNewCandidates;
		}

		// Eliminate based on no closed contours
		for (int ccc = 0; ccc < vecNoClosedContourOp.size(); ccc++) {
			std::set<int> setNewCandidates;

			// Load the search variable data
			DataVector<float> dataState(grid.GetSize());
			
			vecNoClosedContourOp[ccc].m_var.LoadGridData<float>(
				ncInput, grid, dataState, t);

			// Loop through all pressure minima
			std::set<int>::const_iterator iterCandidate
				= setCandidates.begin();

			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				// Determine if a closed contour is present
				bool fHasClosedContour =
					HasClosedContour<float>(
						grid,
						dataState,
						*iterCandidate,
						vecNoClosedContourOp[ccc].m_dDeltaAmount,
						vecNoClosedContourOp[ccc].m_dDistance,
						vecNoClosedContourOp[ccc].m_dMinMaxDist
					);

				// If a closed contour is present, reject this candidate
				if (fHasClosedContour) {
					vecRejectedNoClosedContour[ccc]++;
				} else {
					setNewCandidates.insert(*iterCandidate);
				}
			}

			setCandidates = setNewCandidates;
		}

		Announce("Total candidates: %i", setCandidates.size());
		Announce("Rejected (  latitude): %i", nRejectedLatitude);
		Announce("Rejected (topography): %i", nRejectedTopography);
		Announce("Rejected (    merged): %i", nRejectedMerge);

		for (int tc = 0; tc < vecRejectedThreshold.GetRows(); tc++) {
			Announce("Rejected (thresh. %s): %i",
					vecThresholdOp[tc].m_var.m_strName.c_str(),
					vecRejectedThreshold[tc]);
		}

		for (int ccc = 0; ccc < vecRejectedClosedContour.GetRows(); ccc++) {
			Announce("Rejected (contour %s): %i",
					vecClosedContourOp[ccc].m_var.m_strName.c_str(),
					vecRejectedClosedContour[ccc]);
		}

		for (int ccc = 0; ccc < vecRejectedNoClosedContour.GetRows(); ccc++) {
			Announce("Rejected (nocontour %s): %i",
					vecNoClosedContourOp[ccc].m_var.m_strName.c_str(),
					vecRejectedNoClosedContour[ccc]);
		}

		// Write results to file
		{
			// Parse time information
			//NcVar * varDate = ncInput.get_var("date");
			//NcVar * varDateSec = ncInput.get_var("datesec");

			int nDateYear;
			int nDateMonth;
			int nDateDay;
			int nDateHour;

			NcAtt * attTimeUnits = varTime->get_att("units");
			if (attTimeUnits == NULL) {
				_EXCEPTIONT("Variable \"time\" has no \"units\" attribute");
			}

			std::string strTimeUnits = attTimeUnits->as_string(0);


			std::string strTimeCalendar = "noleap";
			NcAtt * attTimeCalendar = varTime->get_att("calendar");
			if (attTimeUnits != NULL) {
				strTimeCalendar = attTimeCalendar->as_string(0);
			}

			ParseTimeDouble(
				strTimeUnits,
				strTimeCalendar,
				dTime[t],
				nDateYear,
				nDateMonth,
				nDateDay,
				nDateHour);

			// Write time information
			fprintf(fpOutput, "%i\t%i\t%i\t%i\t%i\n",
				nDateYear,
				nDateMonth,
				nDateDay,
				static_cast<int>(setCandidates.size()),
				nDateHour);

			// Write candidate information
			int iCandidateCount = 0;

			// Apply output operators
			DataMatrix<float> dOutput(setCandidates.size(), vecOutputOp.size());
			for (int outc = 0; outc < vecOutputOp.size(); outc++) {

				// Load the search variable data
				DataVector<float> dataState(grid.GetSize());
			
				vecOutputOp[outc].m_var.LoadGridData<float>(
					ncInput, grid, dataState, t);

				// Loop through all pressure minima
				std::set<int>::const_iterator iterCandidate
					= setCandidates.begin();

				iCandidateCount = 0;
				for (; iterCandidate != setCandidates.end(); iterCandidate++) {

					int ixExtremum;
					float dValue;
					float dRMax;

					if (vecOutputOp[outc].m_eOp == OutputOp::Max) {
						FindLocalMinMax<float>(
							grid,
							false,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dValue;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::MaxDist) {
						FindLocalMinMax<float>(
							grid,
							false,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dRMax;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::Min) {
						FindLocalMinMax<float>(
							grid,
							true,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dValue;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::MinDist) {
						FindLocalMinMax<float>(
							grid,
							true,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							ixExtremum,
							dValue,
							dRMax);

						dOutput[iCandidateCount][outc] = dRMax;

					} else if (vecOutputOp[outc].m_eOp == OutputOp::Avg) {
						FindLocalAverage<float>(
							grid,
							dataState,
							*iterCandidate,
							vecOutputOp[outc].m_dDistance,
							dValue);

						dOutput[iCandidateCount][outc] = dValue;

					} else {
						_EXCEPTIONT("Invalid Output operator");
					}

					iCandidateCount++;
				}
			}

			// Output all candidates
			iCandidateCount = 0;

			std::set<int>::const_iterator iterCandidate = setCandidates.begin();
			for (; iterCandidate != setCandidates.end(); iterCandidate++) {

				if (grid.m_nGridDim.size() == 1) {
					fprintf(fpOutput, "\t%i", *iterCandidate);

				} else if (grid.m_nGridDim.size() == 2) {
					fprintf(fpOutput, "\t%i\t%i",
						(*iterCandidate) % static_cast<int>(grid.m_nGridDim[1]),
						(*iterCandidate) / static_cast<int>(grid.m_nGridDim[1]));
				}

				fprintf(fpOutput, "\t%3.6f\t%3.6f",
					grid.m_dLon[*iterCandidate] * 180.0 / M_PI,
					grid.m_dLat[*iterCandidate] * 180.0 / M_PI);

				for (int outc = 0; outc < vecOutputOp.size(); outc++) {
					fprintf(fpOutput, "\t%3.6e",
						dOutput[iCandidateCount][outc]);
				}

				fprintf(fpOutput, "\n");

				iCandidateCount++;
			}
		}

		AnnounceEndBlock("Done");
	}

	fclose(fpOutput);

	ncInput.close();

	AnnounceEndBlock("Done");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

