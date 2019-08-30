///////////////////////////////////////////////////////////////////////////////
///
///	\file    Variable.cpp
///	\author  Paul Ullrich
///	\version July 22, 2018
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
#include "STLStringHelper.h"

#include <set>
#include <cctype>

///////////////////////////////////////////////////////////////////////////////
// VariableRegistry
///////////////////////////////////////////////////////////////////////////////

VariableRegistry::VariableRegistry() {
	m_domDataOp.Add("_VECMAG");
	m_domDataOp.Add("_ABS");
	m_domDataOp.Add("_SIGN");
	m_domDataOp.Add("_ALLPOS");
	m_domDataOp.Add("_SUM");
	m_domDataOp.Add("_AVG");
	m_domDataOp.Add("_DIFF");
	m_domDataOp.Add("_MULT");
	m_domDataOp.Add("_DIV");
	m_domDataOp.Add("_LAT");
	m_domDataOp.Add("_F");
}

///////////////////////////////////////////////////////////////////////////////

VariableRegistry::~VariableRegistry() {
}

///////////////////////////////////////////////////////////////////////////////

int VariableRegistry::FindOrRegister(
	const Variable & var
) {
	for (int i = 0; i < m_vecVariables.size(); i++) {
		if (m_vecVariables[i] == var) {
			return i;
		}
	}
	m_vecVariables.push_back(var);
	return (m_vecVariables.size()-1);
}

///////////////////////////////////////////////////////////////////////////////

Variable & VariableRegistry::Get(
	VariableIndex varix
) {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}
	return m_vecVariables[varix];
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::UnloadAllGridData() {
	for (int i = 0; i < m_vecVariables.size(); i++) {
		m_vecVariables[i].UnloadGridData();
	}
}

///////////////////////////////////////////////////////////////////////////////

DataOp * VariableRegistry::GetDataOp(
	const std::string & strName
) {
	DataOp * pdo = m_domDataOp.Find(strName);
	if (pdo != NULL) {
		return pdo;
	}

	return m_domDataOp.Add(strName);
}

///////////////////////////////////////////////////////////////////////////////
// Variable
///////////////////////////////////////////////////////////////////////////////

bool Variable::operator==(
	const Variable & var
) {
	if (m_fOp != var.m_fOp) {
		return false;
	}
	if (m_strName != var.m_strName) {
		return false;
	}
	if (m_nSpecifiedDim != var.m_nSpecifiedDim) {
		return false;
	}
	for (int i = 0; i < m_nSpecifiedDim; i++) {
		if (m_iDim[i] != var.m_iDim[i]) {
			return false;
		}
	}
	for (int i = 0; i < m_varArg.size(); i++) {
		if (m_varArg[i] != var.m_varArg[i]) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

int Variable::ParseFromString(
	VariableRegistry & varreg,
	const std::string & strIn
) {
	m_fOp = false;
	m_strName = "";
	m_nSpecifiedDim = 0;
	m_varArg.clear();

	bool fParamMode = false;
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
				if (fParamMode) {
					_EXCEPTIONT("Unbalanced curly brackets in variable");
				}
				m_strName = strIn;
				return n;
			}

			// Items in curly brackets are included in variable name
			if (fParamMode) {
				if (strIn[n] == '(') {
					_EXCEPTIONT("Unexpected \'(\' in variable");
				}
				if (strIn[n] == ')') {
					_EXCEPTIONT("Unexpected \')\' in variable");
				}
				if (strIn[n] == '{') {
					_EXCEPTIONT("Unexpected \'{\' in variable");
				}
				if (strIn[n] == '}') {
					fParamMode = false;
				}
				continue;
			}
			if (strIn[n] == '{') {
				fParamMode = true;
				continue;
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

			// No arguments
			if (strIn[n] == ')') {
				return (n+1);
			}

			// Check for floating point argument
			if (isdigit(strIn[n]) || (strIn[n] == '.') || (strIn[n] == '-')) {
				int nStart = n;
				for (; n <= strIn.length(); n++) {
					if (n == strIn.length()) {
						_EXCEPTION1("Op argument list must be terminated"
							" with ): %s", strIn.c_str());
					}
					if ((strIn[n] == ',') || (strIn[n] == ')')) {
						break;
					}
				}

				std::string strFloat = strIn.substr(nStart, n-nStart);
				if (!STLStringHelper::IsFloat(strFloat)) {
					_EXCEPTION2("Invalid floating point number at position %i in: %s",
						nStart, strIn.c_str());
				}
				m_strArg.push_back(strFloat);
				m_varArg.push_back(InvalidVariableIndex);
				m_nSpecifiedDim++;

			// Check for string argument
			} else if (strIn[n] == '\"') {
				int nStart = n;
				for (; n <= strIn.length(); n++) {
					if (n == strIn.length()) {
						_EXCEPTION1("String must be terminated with \": %s",
							strIn.c_str());
					}
					if (strIn[n] == '\"') {
						break;
					}
				}
				if (n >= strIn.length()-1) {
					_EXCEPTION1("Op argument list must be terminated"
						" with ): %s", strIn.c_str());
				}
				if ((strIn[n+1] != ',') && (strIn[n+1] != ')')) {
					_EXCEPTION2("Invalid character in argument list after "
						"string at position %i in: %s",
						n+1, strIn.c_str());
				}

				m_strArg.push_back(strIn.substr(nStart+1,n-nStart-1));
				m_varArg.push_back(InvalidVariableIndex);
				m_nSpecifiedDim++;

			// Check for variable
			} else {
				Variable var;

				n += var.ParseFromString(varreg, strIn.substr(n));

				m_strArg.push_back("");
				m_varArg.push_back(varreg.FindOrRegister(var));
				m_nSpecifiedDim++;
			}

			if (strIn[n] == ')') {
				return (n+1);
			}
		}
	}

	_EXCEPTION1("Malformed variable string \"%s\"", strIn.c_str());
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::ToString(
	VariableRegistry & varreg
) const {
	char szBuffer[20];
	std::string strOut = m_strName;
	if (m_nSpecifiedDim == 0) {
		return strOut;
	}
	strOut += "(";
	for (int d = 0; d < m_nSpecifiedDim; d++) {
		if (m_fOp) {
			if (m_varArg[d] != InvalidVariableIndex) {
				Variable & var = varreg.Get(m_varArg[d]);
				strOut += var.ToString(varreg);
			} else {
				strOut += m_strArg[d];
			}
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

///////////////////////////////////////////////////////////////////////////////

NcVar * Variable::GetFromNetCDF(
	NcFileVector & vecFiles,
	int iTime
) {
	if (m_fOp) {
		_EXCEPTION1("Cannot call GetFromNetCDF() on operator \"%s\"",
			m_strName.c_str());
	}

	// Find the NcVar in all open NcFiles
	NcVar * var = NULL;
	for (int i = 0; i < vecFiles.size(); i++) {
		var = vecFiles[i]->get_var(m_strName.c_str());
		if (var != NULL) {
			break;
		}
	}
	if (var == NULL) {
		_EXCEPTION1("Variable \"%s\" not found in input files",
			m_strName.c_str());
	}

	// If the first dimension of the variable is not "time" then
	// ignore any time index that has been specified.
	int nVarDims = var->num_dims();
	if ((nVarDims > 0) && (iTime != (-1))) {
		if (strcmp(var->get_dim(0)->name(), "time") != 0) {
			iTime = (-1);
			m_fNoTimeInNcFile = true;
		} else {
			if (var->get_dim(0)->size() == 1) {
				iTime = 0;
				m_fNoTimeInNcFile = true;
			} else if (iTime >= var->get_dim(0)->size()) {
				_EXCEPTIONT("Requested time index larger than available in input files:\n"
					"Likely mismatch in time dimension lengths among files");
			}
		}
	}

	// Verify correct dimensionality
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

///////////////////////////////////////////////////////////////////////////////

void Variable::LoadGridData(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	int iTime
) {
	// Check if data already loaded
	if (iTime == (-2)) {
		_EXCEPTIONT("Invalid time index");
	}
	if (iTime == m_iTime) {
		if (m_data.GetRows() != grid.GetSize()) {
			_EXCEPTIONT("Logic error");
		}
		return;
	}
	if ((m_fNoTimeInNcFile) && (m_iTime != (-2))) {
		if (m_data.GetRows() != grid.GetSize()) {
			_EXCEPTIONT("Logic error");
		}
		return;
	}

	//std::cout << "Loading " << ToString(varreg) << " " << iTime << std::endl;

	// Allocate data
	m_data.Allocate(grid.GetSize());
	m_iTime = iTime;

	// Get the data directly from a variable
	if (!m_fOp) {
		// Get pointer to variable
		NcVar * var = GetFromNetCDF(vecFiles, iTime);
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
		var->get(&(m_data[0]), &(nDataSize[0]));

		NcError err;
		if (err.get_err() != NC_NOERR) {
			_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
		}

		return;
	}

	{
		// Get the associated operator
		DataOp * pop = varreg.GetDataOp(m_strName);
		if (pop == NULL) {
			_EXCEPTION1("Unexpected operator \"%s\"", m_strName.c_str());
		}

		// Build argument list
		std::vector<DataArray1D<float> const *> vecArgData;
		for (int i = 0; i < m_varArg.size(); i++) {
			if (m_varArg[i] != InvalidVariableIndex) {
				Variable & var = varreg.Get(m_varArg[i]);
				var.LoadGridData(varreg, vecFiles, grid, iTime);

				vecArgData.push_back(&var.GetData());
			} else {
				vecArgData.push_back(NULL);
			}
		}

		// Apply the DataOp
		pop->Apply(grid, m_strArg, vecArgData, m_data);
	}
/*
	// Evaluate the mean operator
	} else if (m_strName == "_MEAN") {
		if (m_varArg.size() != 2) {
			_EXCEPTION1("_MEAN expects two arguments: %i given",
				m_varArg.size());
		}

		// Obtain field and distance
		Variable & varField = varreg.Get(m_varArg[0]);
		Variable & varDist = varreg.Get(m_varArg[1]);

		varField.LoadGridData(varreg, vecFiles, grid, iTime);

		// Load distance (in degrees) and convert to radians
		double dDist = atof(varDist.m_strName.c_str());

		if ((dDist < 0.0) || (dDist > 360.0)) {
			_EXCEPTION1("Distance argument in _MEAN out of range\n"
				"Expected [0,360], found %1.3e", dDist);
		}

		// Calculate mean of field
		m_data.Zero();

		if (grid.m_vecConnectivity.size() != m_data.GetRows()) {
			_EXCEPTIONT("Invalid grid connectivity array");
		}

		for (int i = 0; i < m_data.GetRows(); i++) {
			std::set<int> setNodesVisited;
			std::set<int> setNodesToVisit;
			setNodesToVisit.insert(i);

			double dLat0 = grid.m_dLat[i];
			double dLon0 = grid.m_dLon[i];

			while (setNodesToVisit.size() != 0) {

				// Next node to explore
				int j = *(setNodesToVisit.begin());

				setNodesToVisit.erase(setNodesToVisit.begin());
				setNodesVisited.insert(j);

				// Update the mean
				m_data[i] += varField.m_data[j];

				// Find additional neighbors to explore
				for (int k = 0; k < grid.m_vecConnectivity[j].size(); k++) {
					int l = grid.m_vecConnectivity[j][k];

					// Check if already visited
					if (setNodesVisited.find(l) != setNodesVisited.end()) {
						continue;
					}

					// Test that this node satisfies the distance criteria
					double dLat1 = grid.m_dLat[l];
					double dLon1 = grid.m_dLon[l];

					double dR =
						sin(dLat0) * sin(dLat1)
						+ cos(dLat0) * cos(dLat1) * cos(dLon1 - dLon0);

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

					if (dR > dDist) {
						continue;
					}

					// Add node to visit
					setNodesToVisit.insert(l);
				}
			}

			// Average data to obtain mean
			if (setNodesVisited.size() == 0) {
				_EXCEPTIONT("Logic error");
			}

			m_data[i] /= static_cast<float>(setNodesVisited.size());
		}

	} else {
		_EXCEPTION1("Unexpected operator \"%s\"", m_strName.c_str());
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void Variable::UnloadGridData() {

	// Force data to be loaded within this structure
	m_iTime = (-2);
}

///////////////////////////////////////////////////////////////////////////////

