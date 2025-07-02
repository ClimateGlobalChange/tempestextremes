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
#include "NetCDFUtilities.h"
#include "STLStringHelper.h"
#include "Units.h"

#include <set>
#include <cctype>

///////////////////////////////////////////////////////////////////////////////
// VariableRegistry
///////////////////////////////////////////////////////////////////////////////

VariableRegistry::VariableRegistry() :
	m_sProcessingQueueVarPos(-1)
{
	m_domDataOp.Add("_VECMAG");
	m_domDataOp.Add("_ABS");
	m_domDataOp.Add("_SIGN");
	m_domDataOp.Add("_ALLPOS");
	m_domDataOp.Add("_SUM");
	m_domDataOp.Add("_AVG");
	m_domDataOp.Add("_DIFF");
	m_domDataOp.Add("_PROD");
	m_domDataOp.Add("_DIV");
	m_domDataOp.Add("_MIN");
	m_domDataOp.Add("_MAX");
	m_domDataOp.Add("_COND");
	m_domDataOp.Add("_SQRT");
	m_domDataOp.Add("_POW");
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
				pvar->m_varArg.push_back(InvalidVariableIndex);
				if (strDim == ":") {
					pvar->m_fFreeArg.push_back(true);
				} else {
					pvar->m_fFreeArg.push_back(false);
				}
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
				pvar->m_fFreeArg.push_back(false);

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
				pvar->m_fFreeArg.push_back(false);

			// Check for variable
			} else {
				VariableIndex varix;
				n += FindOrRegisterSubStr(strIn.substr(n), &varix);

				pvar->m_strArg.push_back("");
				pvar->m_varArg.push_back(varix);
				pvar->m_fFreeArg.push_back(false);
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
		_EXCEPTION1("Variable index (%i) out of range", varix);
	}
	return *(m_vecVariables[varix]);
}

///////////////////////////////////////////////////////////////////////////////

std::string VariableRegistry::GetVariableString(
	VariableIndex varix
) const {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTION1("Variable index (%i) out of range", varix);
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
	if (var->num_dims() > 0) {
		if (NcIsTimeDimension(var->get_dim(0))) {
			lBegin++;
		}
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

void VariableRegistry::GetAuxiliaryDimInfo(
	const NcFileVector & vecncDataFiles,
	const SimpleGrid & grid,
	VariableIndex varix,
	DimInfoVector & vecAuxDimInfo
) const {
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}

	// Recursively add dimension information
	if (m_vecVariables[varix]->IsOp()) {
		const VariableIndexVector & varArg =
			m_vecVariables[varix]->GetArgumentVarIxs();
		for (int i = 0; i < varArg.size(); i++) {
			if (varArg[i] != InvalidVariableIndex) {
				DimInfoVector vecAuxDimInfoRecurse;
				GetAuxiliaryDimInfo(
					vecncDataFiles,
					grid,
					varArg[i],
					vecAuxDimInfoRecurse);

				if (vecAuxDimInfo.size() == 0) {
					vecAuxDimInfo = vecAuxDimInfoRecurse;
				} else if (vecAuxDimInfo != vecAuxDimInfoRecurse) {
					_EXCEPTION3("Incompatible variable \"%s\" in operator: expected dimensions %s, found %s",
						GetVariableString(varArg[i]).c_str(),
						vecAuxDimInfo.ToString().c_str(),
						vecAuxDimInfoRecurse.ToString().c_str());
				}
			}
		}

	// Reached a Variable that is not an operator; check if it exists already
	// in the list and if not add it.
	} else {

		_ASSERT((varix >= 0) && (varix < m_vecVariables.size()));
		const Variable * pvar = m_vecVariables[varix];

		DimInfoVector vecAuxDimInfoRecurse;
		GetAuxiliaryDimInfo(
			vecncDataFiles,
			grid,
			pvar->GetName(),
			vecAuxDimInfoRecurse);

		if (vecAuxDimInfoRecurse.size() != pvar->m_fFreeArg.size()) {
			_EXCEPTION3("Incompatible variable \"%s\": File variable contains %lu auxiliary dimensions but %lu specified",
				GetVariableString(varix).c_str(),
				vecAuxDimInfoRecurse.size(),
				pvar->m_fFreeArg.size());
		}

		DimInfoVector vecAuxDimInfoRecurse_FreeDimOnly;
		for (size_t a = 0; a < pvar->m_fFreeArg.size(); a++) {
			if (pvar->m_fFreeArg[a]) {
				vecAuxDimInfoRecurse_FreeDimOnly.push_back(vecAuxDimInfoRecurse[a]);
			}
		}

		if (vecAuxDimInfo.size() == 0) {
			vecAuxDimInfo = vecAuxDimInfoRecurse_FreeDimOnly;
		} else if (vecAuxDimInfo != vecAuxDimInfoRecurse_FreeDimOnly) {
			_EXCEPTION3("Incompatible variable \"%s\": expected %s, found %s",
				GetVariableString(varix).c_str(),
				vecAuxDimInfo.ToString().c_str(),
				vecAuxDimInfoRecurse_FreeDimOnly.ToString().c_str());
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::PopulateFillValues(
	const NcFileVector & vecncDataFiles
) {
	for (Variable * pvar : m_vecVariables) {
		_ASSERT(pvar != NULL);
		if ((!pvar->m_data.HasFillValue()) && (!pvar->IsOp())) {
			NcVar * ncvar = NULL;
			vecncDataFiles.FindContainingVariable(pvar->GetName(), &ncvar);
			if (ncvar != NULL) {
				NcAtt * attFillValue = ncvar->get_att("_FillValue");
				if (attFillValue == NULL) {
					attFillValue = ncvar->get_att("missing_value");
				}
				if (attFillValue != NULL) {
					pvar->m_data.SetFillValue(attFillValue->as_float(0));
				}

				NcAtt * attUnits = ncvar->get_att("units");
				if (attUnits != NULL) {
					pvar->m_data.SetUnits(attUnits->as_string(0));
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::AssignAuxiliaryIndicesRecursive(
	Variable & var,
	const std::vector<std::string> & vecArg
) {
	var.UnloadGridData();

	// Assign free indices to this variable
	if (!var.m_fOp) {
		_ASSERT(var.m_fFreeArg.size() == var.m_strArg.size());
		_ASSERT(var.m_fFreeArg.size() == var.m_varArg.size());
		size_t sFreeArgs = 0;
		for (int d = 0; d < var.m_fFreeArg.size(); d++) {
			if (var.m_fFreeArg[d]) {
				sFreeArgs++;
			}
		}
		if (sFreeArgs != vecArg.size()) {
			_EXCEPTION3("Incompatible auxiliary indices in variable \"%s\": Expected %lu indices, received %lu indcies",
				var.GetName().c_str(), sFreeArgs, vecArg.size());
		}

		// Copy over free arguments
		size_t a = 0;
		for (size_t d = 0; d < var.m_fFreeArg.size(); d++) {
			if (var.m_fFreeArg[d]) {
				var.m_strArg[d] = vecArg[a];
				var.m_varArg[d] = InvalidVariableIndex;
				a++;
			}
		}
	}

	// Recursively assign indices
	for (size_t v = 0; v < var.m_varArg.size(); v++) {
		if (var.m_varArg[v] != InvalidVariableIndex) {
			AssignAuxiliaryIndicesRecursive(
				Get(var.m_varArg[v]),
				vecArg);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::ClearProcessingQueue() {
	m_sProcessingQueueVarPos = (-1);
	m_vecProcessingQueue.clear();
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::AppendVariableToProcessingQueue(
	const NcFileVector & vecncDataFiles,
	const SimpleGrid & grid,
	VariableIndex varix,
	DimInfoVector * pvecAuxDimInfo
) {
	if (m_sProcessingQueueVarPos != (-1)) {
		_EXCEPTIONT("LOGIC ERROR: Processing queue locked");
	}

	DimInfoVector vecAuxDimInfo;
	GetAuxiliaryDimInfo(
		vecncDataFiles,
		grid,
		varix,
		vecAuxDimInfo);

	VariableAuxIndex auxix(vecAuxDimInfo.size());
	for (size_t d = 0; d < vecAuxDimInfo.size(); d++) {
		auxix[d] = vecAuxDimInfo[d].size;
	}
	VariableAuxIndexIterator auxiter;
	auxiter.Initialize(varix, auxix, false);
	m_vecProcessingQueue.push_back(auxiter);

	if (pvecAuxDimInfo != NULL) {
		(*pvecAuxDimInfo) = vecAuxDimInfo;
	}
}

///////////////////////////////////////////////////////////////////////////////

size_t VariableRegistry::GetProcessingQueueVarPos() const {
	return m_sProcessingQueueVarPos;
}

///////////////////////////////////////////////////////////////////////////////

VariableIndex VariableRegistry::GetProcessingQueueVarIx() const {
	_ASSERT(m_sProcessingQueueVarPos < m_vecProcessingQueue.size());

	const VariableAuxIndexIterator & auxit =
		m_vecProcessingQueue[m_sProcessingQueueVarPos];
	
	return auxit.m_varix;
}

///////////////////////////////////////////////////////////////////////////////

Variable & VariableRegistry::GetProcessingQueueVariable() {
	_ASSERT(m_sProcessingQueueVarPos < m_vecProcessingQueue.size());

	const VariableAuxIndexIterator & auxit =
		m_vecProcessingQueue[m_sProcessingQueueVarPos];

	VariableIndex varix = auxit.m_varix;
	if ((varix < 0) || (varix >= m_vecVariables.size())) {
		_EXCEPTIONT("Variable index out of range");
	}

	Variable & var = *(m_vecVariables[varix]);

	// Assign auxiliary indices and clear existing data
	if (auxit.m_vecValue.size() != 0) {
		std::vector<std::string> vecArg(auxit.m_vecValue.size());
		for (size_t d = 0; d < auxit.m_vecValue.size(); d++) {
			vecArg[d] = std::to_string(auxit.m_vecValue[d]);
		}
		AssignAuxiliaryIndicesRecursive(var, vecArg);
		var.UnloadGridData();
	}

	return var;
}

///////////////////////////////////////////////////////////////////////////////

const VariableAuxIndex & VariableRegistry::GetProcessingQueueAuxIx() const {
	_ASSERT(m_sProcessingQueueVarPos < m_vecProcessingQueue.size());

	const VariableAuxIndexIterator & auxit =
		m_vecProcessingQueue[m_sProcessingQueueVarPos];

	return auxit.m_vecValue;
}

///////////////////////////////////////////////////////////////////////////////

const VariableAuxIndex & VariableRegistry::GetProcessingQueueAuxSize() const {
	_ASSERT(m_sProcessingQueueVarPos < m_vecProcessingQueue.size());

	const VariableAuxIndexIterator & auxit =
		m_vecProcessingQueue[m_sProcessingQueueVarPos];

	return auxit.m_vecSize;
}

///////////////////////////////////////////////////////////////////////////////

size_t VariableRegistry::GetProcessingQueueOffset() const {
	_ASSERT(m_sProcessingQueueVarPos < m_vecProcessingQueue.size());

	const VariableAuxIndexIterator & auxit =
		m_vecProcessingQueue[m_sProcessingQueueVarPos];

	return auxit.Offset();
}

///////////////////////////////////////////////////////////////////////////////

bool VariableRegistry::AdvanceProcessingQueue() {
	if (m_sProcessingQueueVarPos == (-1)) {
		m_sProcessingQueueVarPos = 0;
		if (m_vecProcessingQueue.size() == 0) {
			return false;
		}
		return true;
	}
	if (m_sProcessingQueueVarPos >= m_vecProcessingQueue.size()) {
		return false;
	}

	_ASSERT(!m_vecProcessingQueue[m_sProcessingQueueVarPos].at_end());

	m_vecProcessingQueue[m_sProcessingQueueVarPos]++;

	if (m_vecProcessingQueue[m_sProcessingQueueVarPos].at_end()) {
		m_sProcessingQueueVarPos++;
		if (m_sProcessingQueueVarPos >= m_vecProcessingQueue.size()) {
			return false;
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void VariableRegistry::ResetProcessingQueue() {
	m_sProcessingQueueVarPos = (-1);
	for (size_t d = 0; d < m_vecProcessingQueue.size(); d++) {
		m_vecProcessingQueue[d].Reset();
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
) const {
	if (m_fOp != var.m_fOp) {
		return false;
	}
	if (m_strName != var.m_strName) {
		return false;
	}
	if (m_strArg.size() != var.m_strArg.size()) {
		return false;
	}
	_ASSERT(m_strArg.size() == m_varArg.size());
	_ASSERT(m_strArg.size() == var.m_varArg.size());
	for (int i = 0; i < m_strArg.size(); i++) {
		if (m_strArg[i] != var.m_strArg[i]) {
			return false;
		}
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
	_ASSERT(m_varArg.size() == m_strArg.size());
	char szBuffer[20];
	std::string strOut = m_strName;
	if (m_varArg.size() != 0) {
		strOut += "(";
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

	// Get the time index
	long lTime;
	if (!NcIsTimeDimension(var->get_dim(0))) {
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
		_EXCEPTION4("Variable \"%s\" in file \"%s\" has inconsistent auxiliary dimensions (%i) than specified on command line (%i)",
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

		// Get the integer index using either index notation or value-based notation
		lDim[nSetDims] =
			GetIntegerIndexFromValueBasedIndex(
				ncfilevec[sPos],
				ncfilevec.GetFilename(sPos),
				m_strName,
				nSetDims,
				m_strArg[a]);

		nSetDims++;
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

	std::string strSourceFilenamesArg = vecFiles.GetConcatenatedFilenames();
	if (strSourceFilenamesArg == m_strSourceFilenames) {
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
					" \"%s\" on \"ncol\" -- possible mismatch between connectivity file and data",
					m_strName.c_str());
			}

			nDataSize[nVarDims-1] = nSize;
		}

		// Load the data
		var->get(&(m_data[0]), &(nDataSize[0]));

		// Only check attributes when files are changed
		if (m_strSourceFilenames != strSourceFilenamesArg) {

			// Turn off errors as we check attributes
			NcError err(NcError::silent_nonfatal);
			if (err.get_err() != NC_NOERR) {
				_EXCEPTION1("NetCDF Fatal Error (%i)", err.get_err());
			}

			// Check for _FillValue
			NcAtt * attFillValue = var->get_att("_FillValue");
			if (attFillValue == NULL) {
				attFillValue = var->get_att("missing_value");
			}
			if (attFillValue == NULL) {
				m_data.RemoveFillValue();
			} else {
				m_data.SetFillValue(attFillValue->as_float(0));
			}

			// Check for units
			NcAtt * attUnits = var->get_att("units");
			if (attUnits == NULL) {
				m_data.SetUnits("");
			} else {
				m_data.SetUnits(attUnits->as_string(0));
			}

			// Check for scale_factor attribute
			NcAtt * attScaleFactor = var->get_att("scale_factor");
			NcAtt * attAddOffset = var->get_att("add_offset");
			if ((attScaleFactor == NULL) && (attAddOffset == NULL)) {
				m_fHasScaleFactorOrAddOffset = false;
			} else {
				m_fHasScaleFactorOrAddOffset = true;
				if (attScaleFactor == NULL) {
					m_dScaleFactor = 1.0f;
				} else {
					m_dScaleFactor = attScaleFactor->as_float(0);
					if (m_data.HasFillValue()) {
						m_data.SetFillValue(m_data.GetFillValue() * m_dScaleFactor);
					}
				}
				if (attAddOffset == NULL) {
					m_dAddOffset = 0.0f;
				} else {
					m_dAddOffset = attAddOffset->as_float(0);
					if (m_data.HasFillValue()) {
						m_data.SetFillValue(m_data.GetFillValue() + m_dAddOffset);
					}
				}
			}
		}

		// Apply scale_factor and add_offset
		if (m_fHasScaleFactorOrAddOffset) {
			for (int i = 0; i < m_data.GetRows(); i++) {
				m_data[i] = m_data[i] * m_dScaleFactor + m_dAddOffset;
			}
		}

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
	}

	// Store the time
	m_timeStored = time;

	// Store the filenames
	m_strSourceFilenames = strSourceFilenamesArg;
}

///////////////////////////////////////////////////////////////////////////////

void Variable::UnloadGridData() {

	// Force data to be loaded within this structure
	m_timeStored = Time(Time::CalendarUnknown);
	m_strSourceFilenames = "";
}

///////////////////////////////////////////////////////////////////////////////

