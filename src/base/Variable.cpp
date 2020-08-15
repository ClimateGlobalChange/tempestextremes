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
#include "Units.h"

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
	for (int v = 0; v < m_vecVariables.size(); v++) {
		delete m_vecVariables[v];
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::InsertUniqueOrDelete(
	Variable * pvar,
	VariableIndex * pvarix
) {
	for (int v = 0; v < m_vecVariables.size(); v++) {
		if (*pvar == *(m_vecVariables[v])) {
			delete pvar;
			if (pvarix != NULL) {
				*pvarix = v;
			}
			return;
		}
	}

	m_vecVariables.push_back(pvar);
	if (pvarix != NULL) {
		*pvarix = m_vecVariables.size()-1;
	}
}

///////////////////////////////////////////////////////////////////////////////

int VariableRegistry::FindOrRegisterSubStr(
	const std::string & strIn,
	VariableIndex * pvarix
) {
	// Does not exist
	Variable * pvar = new Variable();
	_ASSERT(pvar != NULL);

	if (pvarix != NULL) {
		*pvarix = InvalidVariableIndex;
	}

	// Parse the string
	bool fParamMode = false;
	bool fDimMode = false;
	std::string strDim;

	if (strIn.length() >= 1) {
		if (strIn[0] == '_') {
			pvar->m_fOp = true;
		}
	}

	for (int n = 0; n <= strIn.length(); n++) {

		// Reading the variable name
		if (!fDimMode) {
			if (n == strIn.length()) {
				if (fParamMode) {
					_EXCEPTIONT("Unbalanced curly brackets in variable");
				}
				pvar->m_strName = strIn;
				InsertUniqueOrDelete(pvar, pvarix);
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
				pvar->m_strName = strIn.substr(0, n);
				InsertUniqueOrDelete(pvar, pvarix);
				return n;
			}
			if (strIn[n] == '(') {
				pvar->m_strName = strIn.substr(0, n);
				fDimMode = true;
				continue;
			}
			if (strIn[n] == ')') {
				pvar->m_strName = strIn.substr(0, n);
				InsertUniqueOrDelete(pvar, pvarix);
				return n;
			}

		// Reading in dimensions
		} else if (!pvar->m_fOp) {
			if (pvar->m_strArg.size() >= MaxVariableArguments) {
				_EXCEPTION1("Sanity check fail: Only %i dimensions / arguments may "
					"be specified", MaxVariableArguments);
			}
			if (n == strIn.length()) {
				_EXCEPTION1("Variable dimension list must be terminated"
					" with ): %s", strIn.c_str());
			}
			if ((strIn[n] == ',') || (strIn[n] == ')')) {
				if (strDim.length() == 0) {
					_EXCEPTIONT("Invalid dimension index in variable");
				}
				pvar->m_strArg.push_back(strDim);
				//pvar->m_lArg.push_back(std::stol(strDim));
				pvar->m_varArg.push_back(InvalidVariableIndex);
				strDim = "";

				if (strIn[n] == ')') {
					InsertUniqueOrDelete(pvar, pvarix);
					return (n+1);
				}

			} else {
				strDim += strIn[n];
			}

		// Reading in arguments
		} else {
			if (pvar->m_strArg.size() >= MaxVariableArguments) {
				_EXCEPTION1("Sanity check fail: Only %i dimensions / arguments may "
					"be specified", MaxVariableArguments);
			}
			if (n == strIn.length()) {
				_EXCEPTION1("Op argument list must be terminated"
					" with ): %s", strIn.c_str());
			}

			// No arguments
			if (strIn[n] == ')') {
				InsertUniqueOrDelete(pvar, pvarix);
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
				pvar->m_strArg.push_back(strFloat);
				pvar->m_varArg.push_back(InvalidVariableIndex);
				//pvar->m_lArg.push_back(Variable::InvalidArgument);

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

				pvar->m_strArg.push_back(strIn.substr(nStart+1,n-nStart-1));
				pvar->m_varArg.push_back(InvalidVariableIndex);
				//pvar->m_lArg.push_back(Variable::InvalidArgument);

			// Check for variable
			} else {
				VariableIndex varix;
				n += FindOrRegisterSubStr(strIn.substr(n), &varix);

				pvar->m_strArg.push_back("");
				pvar->m_varArg.push_back(varix);
				//pvar->m_lArg.push_back(Variable::InvalidArgument);
			}

			if (strIn[n] == ')') {
				InsertUniqueOrDelete(pvar, pvarix);
				return (n+1);
			}
		}
	}

	_EXCEPTION1("Malformed variable string \"%s\"", strIn.c_str());
}

///////////////////////////////////////////////////////////////////////////////

VariableIndex VariableRegistry::FindOrRegister(
	const std::string & strIn
) {
	VariableIndex varix;
	int iFinalStringPos = FindOrRegisterSubStr(strIn, &varix);

	if (iFinalStringPos != strIn.length()) {
		_EXCEPTION1("Malformed variable: Extra characters found at end of string \"%s\"",
			strIn.c_str());
	}

	_ASSERT((varix >= 0) && (varix < m_vecVariables.size()));

	return varix;
}

///////////////////////////////////////////////////////////////////////////////

Variable & VariableRegistry::Get(
	VariableIndex varix
) {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}
	return *(m_vecVariables[varix]);
}

///////////////////////////////////////////////////////////////////////////////

std::string VariableRegistry::GetVariableString(
	VariableIndex varix
) const {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}
	return m_vecVariables[varix]->ToString(*this);
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::UnloadAllGridData() {
	for (int i = 0; i < m_vecVariables.size(); i++) {
		m_vecVariables[i]->UnloadGridData();
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::GetDependentVariableIndicesRecurse(
	VariableIndex varix,
	std::vector<VariableIndex> & vecDependentIxs
) const {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}

	// Recursively add dependent variable indices
	if (m_vecVariables[varix]->IsOp()) {
		const VariableIndexVector & varArg =
			m_vecVariables[varix]->GetArgumentVarIxs();
		for (int i = 0; i < varArg.size(); i++) {
			if (varArg[i] != InvalidVariableIndex) {
				GetDependentVariableIndicesRecurse(varArg[i], vecDependentIxs);
			}
		}

	// Reached a Variable that is not an operator; check if it exists already
	// in the list and if not add it.
	} else {
		for (int i = 0; i < vecDependentIxs.size(); i++) {
			if (vecDependentIxs[i] == varix) {
				return;
			}
		}
		vecDependentIxs.push_back(varix);
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::GetDependentVariableIndices(
	VariableIndex varix,
	std::vector<VariableIndex> & vecDependentIxs
) const {
	vecDependentIxs.clear();
	GetDependentVariableIndicesRecurse(varix, vecDependentIxs);
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::GetDependentVariableNames(
	VariableIndex varix,
	std::vector<std::string> & vecDependentVarNames
) const {
	vecDependentVarNames.clear();
	std::vector<VariableIndex> vecDependentIxs;
	GetDependentVariableIndices(varix, vecDependentIxs);
	for (int i = 0; i < vecDependentIxs.size(); i++) {
		vecDependentVarNames.push_back(
			GetVariableString(vecDependentIxs[i]));
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::GetAuxiliaryDimInfo(
	const NcFileVector & vecncDataFiles,
	const SimpleGrid & grid,
	const std::string & strVarName,
	DimInfoVector & vecAuxDimInfo
) {
	// Find the first occurrence of this variable in all open NcFiles
	NcVar * var = NULL;
	for (int i = 0; i < vecncDataFiles.size(); i++) {
		var = vecncDataFiles[i]->get_var(strVarName.c_str());
		if (var != NULL) {
			break;
		}
	}
	if (var == NULL) {
		_EXCEPTION1("Variable \"%s\" not found in input files",
			strVarName.c_str());
	}

	// First auxiliary dimension
	long lBegin = 0;
	long lEnd = var->num_dims() - grid.DimCount();

	// If the first dimension is time then ignore it.
	if ((var->num_dims() > 0) && (strcmp(var->get_dim(0)->name(), "time") == 0)) {
		lBegin++;
	}

	// Ignore grid dimensions at the end
	if (lEnd - lBegin < 0) {
		_EXCEPTION1("Missing spatial dimensions in variable \"%s\"",
			strVarName.c_str());
	}

	// Story auxiliary sizes
	vecAuxDimInfo.resize(lEnd - lBegin);
	for (long d = lBegin; d < lEnd; d++) {
		vecAuxDimInfo[d-lBegin].name = var->get_dim(d)->name();
		vecAuxDimInfo[d-lBegin].size = var->get_dim(d)->size();
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::GetAuxiliaryDimInfoAndVerifyConsistency(
	const NcFileVector & vecncDataFiles,
	const SimpleGrid & grid,
	const std::vector<std::string> & vecVariables,
	DimInfoVector & vecAuxDimInfo
) {
	if (vecVariables.size() == 0) {
		_EXCEPTIONT("Input vector of variable names must contain at least 1 entry");
	}

	GetAuxiliaryDimInfo(
		vecncDataFiles,
		grid,
		vecVariables[0],
		vecAuxDimInfo);

	DimInfoVector vecAuxDimInfoOther;
	for (int v = 1; v < vecVariables.size(); v++) {
		GetAuxiliaryDimInfo(
			vecncDataFiles,
			grid,
			vecVariables[v],
			vecAuxDimInfoOther);

		if (vecAuxDimInfo.size() != vecAuxDimInfoOther.size()) {
			_EXCEPTION4("Incompatible base variables \"%s\" and \"%s\": Disagreement in number of dimensions (%li vs %li)",
				vecVariables[0].c_str(),
				vecVariables[v].c_str(),
				vecAuxDimInfo.size(),
				vecAuxDimInfoOther.size());
		}

		for (int d = 0; d < vecAuxDimInfo.size(); d++) {
			if ((vecAuxDimInfo[d].name != vecAuxDimInfoOther[d].name) ||
			    (vecAuxDimInfo[d].size != vecAuxDimInfoOther[d].size)
			) {
				_EXCEPTION7("Incompatible base variables \"%s\" and \"%s\": Disagreement in dimension %i (%s:%li) vs (%s:%li)",
					vecVariables[0].c_str(),
					vecVariables[v].c_str(),
					d,
					vecAuxDimInfo[d].name.c_str(),
					vecAuxDimInfo[d].size,
					vecAuxDimInfoOther[d].name.c_str(),
					vecAuxDimInfoOther[d].size);
			}
		}
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

const long Variable::InvalidArgument = (-1);

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
	if (m_strArg.size() != var.m_strArg.size()) {
		return false;
	}
	//_ASSERT(m_strArg.size() == m_lArg.size());
	_ASSERT(m_strArg.size() == m_varArg.size());
	//_ASSERT(m_strArg.size() == var.m_lArg.size());
	_ASSERT(m_strArg.size() == var.m_varArg.size());
	for (int i = 0; i < m_strArg.size(); i++) {
		if (m_strArg[i] != var.m_strArg[i]) {
			return false;
		}
		//if (m_lArg[i] != var.m_lArg[i]) {
		//	return false;
		//}
		if (m_varArg[i] != var.m_varArg[i]) {
			return false;
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////

std::string Variable::ToString(
	const VariableRegistry & varreg
) const {
	char szBuffer[20];
	std::string strOut = m_strName;
	if (m_varArg.size() == 0) {
		return strOut;
	}
	strOut += "(";
	_ASSERT(m_varArg.size() == m_strArg.size());
	//_ASSERT(m_varArg.size() == m_lArg.size());
	for (size_t d = 0; d < m_varArg.size(); d++) {
		if (m_varArg[d] != InvalidVariableIndex) {
			strOut += varreg.GetVariableString(m_varArg[d]);
		} else {
			strOut += m_strArg[d];
		}
		if (d != m_varArg.size()-1) {
			strOut += ",";
		} else {
			strOut += ")";
		}
	}
	return strOut;
}

///////////////////////////////////////////////////////////////////////////////

NcVar * Variable::GetNcVarFromNcFileVector(
	const NcFileVector & ncfilevec,
	const SimpleGrid & grid
) {
	if (m_fOp) {
		_EXCEPTION1("Cannot call GetNcVarFromNetCDF() on operator \"%s\"",
			m_strName.c_str());
	}

	// Find the NcVar in all open NcFiles with this name
	NcVar * var;
	size_t sPos =
		ncfilevec.FindContainingVariable(
			m_strName,
			&var);

	if (sPos == NcFileVector::InvalidIndex) {
		_EXCEPTION1("Variable \"%s\" not found in input files",
			m_strName.c_str());
	}

	// Check number of variable dimensions
	int nVarDims = var->num_dims();
	if (nVarDims == 0) {
		_EXCEPTION2("Variable \"%s\" in file \"%s\" has no dimensions; cannot read as a field",
			m_strName.c_str(),
			ncfilevec.GetFilename(sPos).c_str());
	}

	// Get the current time
	m_timeStored = ncfilevec.GetTime();

	// Get the time index
	long lTime;
	std::string strDim0Name = var->get_dim(0)->name();
	if (strDim0Name != "time") {
		lTime = NcFileVector::NoTimeIndex;
		m_fNoTimeInNcFile = true;
	} else {
		lTime = ncfilevec.GetTimeIx(sPos);
		m_fNoTimeInNcFile = false;
	}

	// Verify correct dimensionality
	int nRequestedVarDims = m_strArg.size() + grid.DimCount();
	if (lTime != NcFileVector::NoTimeIndex) {
		nRequestedVarDims++;
	}
	if (nVarDims != nRequestedVarDims) {
		_EXCEPTION4("Variable \"%s\" in file \"%s\" has fewer auxiliary dimensions (%i) than specified on command line (%i)",
			m_strName.c_str(),
			ncfilevec.GetFilename(sPos).c_str(),
			nVarDims,
			nRequestedVarDims);
	}

	// Set the current index position for this variable
	int nSetDims = 0;
	std::vector<long> lDim;
	lDim.resize(nVarDims, 0);

	if (lTime != NcFileVector::NoTimeIndex) {
		lDim[0] = lTime;
		nSetDims++;
	}

	for (int a = 0; a < m_strArg.size(); a++) {

		if (m_strArg[a].length() == 0) {
			_EXCEPTION2("Empty string argument found in position %i of variable %s",
				a, m_strName.c_str());
		}

		// Index notation
		if (STLStringHelper::IsIntegerIndex(m_strArg[a])) {
			lDim[nSetDims] = std::stol(m_strArg[a]);

		// Value notation (value + units)
		} else {
			_ASSERT(nSetDims < var->num_dims());

			// Get the units from the NetCDF file
			NcDim * dimArg = var->get_dim(nSetDims);
			if (dimArg == NULL) {
				_EXCEPTION3("Failure to read dimension %i of variable \"%s\" in file \"%s\"",
					nSetDims,
					m_strName.c_str(),
					ncfilevec.GetFilename(sPos).c_str());
			}

			std::string strDimName = dimArg->name();

			NcVar * varDim = ncfilevec[sPos]->get_var(strDimName.c_str());
			if (varDim == NULL) {
				_EXCEPTION3("Dimension \"%s\" in file \"%s\" has no corresponding variable (unable to use value notation for variable \"%s\")",
					dimArg->name(),
					ncfilevec.GetFilename(sPos).c_str(),
					m_strName.c_str());
			}
			if (varDim->num_dims() != 1) {
				_EXCEPTION2("Dimension variable \"%s\" in file \"%s\" has more than one dimension",
					dimArg->name(),
					ncfilevec.GetFilename(sPos).c_str());
			}
			if (strDimName != varDim->get_dim(0)->name()) {
				_EXCEPTION3("Dimension variable \"%s\" in file \"%s\" instead refers to dimension \"%s\"",
					dimArg->name(),
					ncfilevec.GetFilename(sPos).c_str(),
					varDim->get_dim(0)->name());
			}

			NcAtt * attUnits = varDim->get_att("units");
			if (attUnits == NULL) {
				_EXCEPTION3("Dimension variable \"%s\" in file \"%s\" has no units attribute (unable to use value notation for variable \"%s\")",
					dimArg->name(),
					ncfilevec.GetFilename(sPos).c_str(),
					m_strName.c_str());
			}

			std::string strNcDimUnits = attUnits->as_string(0);

			// Get the value and units from the argument string
			std::string strValue;
			std::string strUnits;
			SplitIntoValueAndUnits(m_strArg[a], strValue, strUnits);

			// Find the corresponding index in the dimension variable

			// int type dimension
			if (varDim->type() == ncInt) {
				if (!STLStringHelper::IsInteger(strValue)) {
					_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is not of type integer, as expected for dimension \"%s\"",
						m_strArg[a].c_str(),
						m_strName.c_str(),
						dimArg->name());
				}
				long lValue = std::stol(strValue);

				bool fConvert =
					ConvertUnits<long>(
						lValue,
						strUnits,
						strNcDimUnits);

				if (!fConvert) {
					_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is incompatible with dimension units \"%s\"",
						m_strArg[a].c_str(),
						m_strName.c_str(),
						strNcDimUnits.c_str());
				}

				DataArray1D<int> nDimValues(dimArg->size());
				varDim->get(&(nDimValues[0]), dimArg->size());

				NcError err;
				if (err.get_err() != NC_NOERR) {
					_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
				}

				long vd = 0;
				for (; vd < dimArg->size(); vd++) {
					if (nDimValues[vd] == lValue) {
						lDim[nSetDims] = vd;
						break;
					}
				}
				if (vd == dimArg->size()) {
					_EXCEPTION2("Dimension \"%s\" does not contain coordinate value \"%s\"",
						dimArg->name(),
						m_strArg[a].c_str());
				}

			// double type dimension
			} else if ((varDim->type() == ncFloat) || (varDim->type() == ncDouble)) {
				if (!STLStringHelper::IsFloat(strValue)) {
					_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is not of type float/double, as expected for dimension \"%s\"",
						m_strArg[a].c_str(),
						m_strName.c_str(),
						dimArg->name());
				}
				double dValue = std::stof(strValue);

				bool fConvert =
					ConvertUnits<double>(
						dValue,
						strUnits,
						strNcDimUnits);

				if (!fConvert) {
					_EXCEPTION3("Specified dimension \"%s\" for variable \"%s\" is incompatible with dimension units \"%s\"",
						m_strArg[a].c_str(),
						m_strName.c_str(),
						strNcDimUnits.c_str());
				}

				DataArray1D<double> dDimValues(dimArg->size());
				if (varDim->type() == ncFloat) {
					DataArray1D<float> flDimValues(dimArg->size());
					varDim->get(&(flDimValues[0]), dimArg->size());
					for (long vd = 0; vd < dimArg->size(); vd++) {
						dDimValues[vd] = static_cast<double>(flDimValues[vd]);
					}
				} else {
					varDim->get(&(dDimValues[0]), dimArg->size());
				}

				NcError err;
				if (err.get_err() != NC_NOERR) {
					_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
				}

				long vd = 0;
				for (; vd < dimArg->size(); vd++) {
					if (fabs(dDimValues[vd] - dValue) < 1.0e-12 * fabs(dValue)) {
						lDim[nSetDims] = vd;
						break;
					}
				}
				if (vd == dimArg->size()) {
					_EXCEPTION2("Dimension \"%s\" does not contain coordinate value \"%s\"",
						dimArg->name(),
						m_strArg[a].c_str());
				}

				printf("\"%s\" corresponds to dimension index %li", m_strArg[a].c_str(), vd);

			// Unknown dimension type
			} else {
				_EXCEPTION2("Dimension variable \"%s\" in file \"%s\" must be of type int, float, or double to use value indexing",
					varDim->name(),
					ncfilevec.GetFilename(sPos).c_str());
			}

			nSetDims++;
		}
	}
/*
	if (m_lArg.size() != 0) {
		for (int i = 0; i < m_lArg.size(); i++) {
			lDim[nSetDims] = m_lArg[i];
			nSetDims++;
		}
	}
*/
/*
	printf("Loading \"%s\" (%s) [", ncfilevec.GetFilename(sPos).c_str(), m_timeStored.ToString().c_str());
	for (int d = 0; d < lDim.size(); d++) {
		printf("%lu",lDim[d]);
		if (d != lDim.size()-1) {
			printf(",");
		}
	}
	printf("]\n");
*/
	var->set_cur(&(lDim[0]));

	NcError err;
	if (err.get_err() != NC_NOERR) {
		_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
	}

	return var;
}

///////////////////////////////////////////////////////////////////////////////

void Variable::LoadGridData(
	VariableRegistry & varreg,
	const NcFileVector & vecFiles,
	const SimpleGrid & grid
) {

	// Check if data already loaded
	const Time & time = vecFiles.GetTime();
	if (time.GetCalendarType() == Time::CalendarUnknown) {
		_EXCEPTIONT("Invalid time specified");
	}
	if (time == m_timeStored) {
		if (m_data.GetRows() != grid.GetSize()) {
			_EXCEPTIONT("Logic error");
		}
		return;
	}
	if ((m_fNoTimeInNcFile) && (m_timeStored.GetCalendarType() != Time::CalendarUnknown)) {
		if (m_data.GetRows() != grid.GetSize()) {
			_EXCEPTIONT("Logic error");
		}
		return;
	}

	//std::cout << "Loading " << ToString(varreg) << " " << lTime << std::endl;

	// Allocate data
	m_data.Allocate(grid.GetSize());
	//m_lTime = lTime;

	// Get the data directly from a variable
	if (!m_fOp) {

		// Get pointer to variable
		NcVar * var = GetNcVarFromNcFileVector(vecFiles, grid);
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

		std::vector<long> nDataSize;
		nDataSize.resize(nVarDims, 1);
		//long nDataSize[7];
		//for (int i = 0; i < 7; i++) {
		//	nDataSize[i] = 1;
		//}

		// Rectilinear grid
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

	// Evaluate a data operator to get the contents of this variable
	} else {
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
				var.LoadGridData(varreg, vecFiles, grid);

				vecArgData.push_back(&var.GetData());
			} else {
				vecArgData.push_back(NULL);
			}
		}

		// Apply the DataOp
		pop->Apply(grid, m_strArg, vecArgData, m_data);

		// Store the time
		m_timeStored = time;
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

		varField.LoadGridData(varreg, vecFiles, grid, lTime);

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
	m_timeStored = Time(Time::CalendarUnknown);
}

///////////////////////////////////////////////////////////////////////////////

