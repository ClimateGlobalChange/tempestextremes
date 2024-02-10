///////////////////////////////////////////////////////////////////////////////
///
///	\file    ThresholdOp.cpp
///	\author  Paul Ullrich
///	\version February 9, 2024
///
///	<remarks>
///		Copyright 2024 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "ThresholdOp.h"
#include "SimpleGrid.h"
#include "DataArray1D.h"

#include <queue>
#include <string>
#include <set>

///////////////////////////////////////////////////////////////////////////////
// ThresholdOp
///////////////////////////////////////////////////////////////////////////////

void ThresholdOp::Parse(
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
				m_dThresholdValue = atof(strSubStr.c_str());

				iLast = i + 1;
				eReadMode = ReadMode_Distance;

			// Read in distance to point that satisfies threshold
			} else if (eReadMode == ReadMode_Distance) {
				m_dDistanceDeg = atof(strSubStr.c_str());

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

	if (m_dDistanceDeg < 0.0) {
		_EXCEPTIONT("For threshold op, distance must be nonnegative");
	}
	if (m_dDistanceDeg > 180.0) {
		_EXCEPTIONT("For threshold op, distance must be less than 180 degrees");
	}

	// Output announcement
	std::string strDescription = varreg.GetVariableString(m_varix);
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
	snprintf(szBuffer, 128, "%f within %f degrees",
		m_dThresholdValue, m_dDistanceDeg);
	strDescription += szBuffer;

	Announce("%s", strDescription.c_str());
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
bool ThresholdOp::IsSatisfied(
	const SimpleGrid & grid,
	const DataArray1D<real> & dataState,
	const int ix0
) const {
	// Special case if m_dDistanceDeg is zero
	if (m_dDistanceDeg < 1.0e-12) {
		double dValue = dataState[ix0];

		if (m_eOp == ThresholdOp::GreaterThan) {
			if (dValue > m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::LessThan) {
			if (dValue < m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::LessThanEqualTo) {
			if (dValue <= m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::EqualTo) {
			if (dValue == m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::NotEqualTo) {
			if (dValue != m_dThresholdValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		return false;
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

		if ((ix != ix0) && (dR > m_dDistanceDeg)) {
			continue;
		}

		// Value at this location
		double dValue = dataState[ix];

		// Apply operator
		if (m_eOp == ThresholdOp::GreaterThan) {
			if (dValue > m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::LessThan) {
			if (dValue < m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::GreaterThanEqualTo) {
			if (dValue >= m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::LessThanEqualTo) {
			if (dValue <= m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::EqualTo) {
			if (dValue == m_dThresholdValue) {
				return true;
			}

		} else if (m_eOp == ThresholdOp::NotEqualTo) {
			if (dValue != m_dThresholdValue) {
				return true;
			}

		} else {
			_EXCEPTIONT("Invalid operation");
		}

		// Add all neighbors of this point
		for (int n = 0; n < grid.m_vecConnectivity[ix].size(); n++) {
			queueNodes.push(grid.m_vecConnectivity[ix][n]);
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
// ThresholdOpTreeNode
///////////////////////////////////////////////////////////////////////////////

void ThresholdOpTreeNode::Parse(
	VariableRegistry & varreg,
	std::string strExpression,
	bool fAllowZeroLengthExpression
) {
	if (m_vecChildren.size() != 0) {
		_EXCEPTIONT("ThresholdOpTreeNode already initialized");
	}

	// Set to a leaf node by default
	m_eType = AlwaysTrueNode;

	// Pre-process expression
	for (;;) {
		// Remove whitespace on either side of expression
		STLStringHelper::RemoveWhitespaceInPlace(strExpression);

		// Zero length expressions are always true
		if (strExpression.length() == 0) {
			if (fAllowZeroLengthExpression) {
				return;
			} else {
				_EXCEPTIONT("Zero-length sub-expression in threshold operator");
			}
		}

		// Remove parentheses if they wrap whole expression, then re-process
		if ((strExpression[0] == '(') && (strExpression[strExpression.length()-1] == ')')) {
			int nParentheses = 0;
			for (int i = 0; i < strExpression.length(); i++) {
				if (strExpression[i] == '(') {
					nParentheses++;
				}
				if (strExpression[i] == ')') {
					if (nParentheses == 0) {
						_EXCEPTION1("Malformed sub-expression in threshold operator \"%s\"", strExpression.c_str());
					}
					nParentheses--;
					if ((nParentheses == 0) && (i != strExpression.length()-1)) {
						nParentheses = -1;
						break;
					}
				}
			}
			if (nParentheses == 0) {
				strExpression = strExpression.substr(1,strExpression.length()-2);
				continue;
			}
		}

		break;
	}

	// Check for operators and balanced parentheses
	// If there are any "or" operators in strExpression then vecOpPositionIx only contains "or" positions
	// This is because "and" operations have precedence and so are evaluated in the children
	std::vector<int> vecOpPositionIx;
	int nParentheses = 0;
	for (int i = 0; i < strExpression.length(); i++) {
		if (strExpression[i] == '(') {
			nParentheses++;
		}
		if (strExpression[i] == ')') {
			if (nParentheses == 0) {
				_EXCEPTION1("Malformed sub-expression in threshold operator \"%s\"", strExpression.c_str());
			}
			nParentheses--;
		}
		if (nParentheses == 0) {
			if (strExpression[i] == '|') {
				if (m_eType != OrNode) {
					m_eType = OrNode;
					vecOpPositionIx.clear();
				}
				vecOpPositionIx.push_back(i);
			}
			if ((strExpression[i] == '&') || (strExpression[i] == ';')) {
				if (m_eType != OrNode) {
					m_eType = AndNode;
					vecOpPositionIx.push_back(i);
				}
			}
		}
	}
	if (nParentheses != 0) {
		_EXCEPTION1("Malformed sub-expression in threshold operator \"%s\"", strExpression.c_str());
	}

	// Operator node; split operations
	if (vecOpPositionIx.size() != 0) {
		int iLast = 0;
		for (int i = 0; i < vecOpPositionIx.size(); i++) {
			_ASSERT(vecOpPositionIx[i] - iLast - 1 >= 0);

			ThresholdOpTreeNode * pChild = new ThresholdOpTreeNode();
			pChild->Parse(
				varreg,
				strExpression.substr(
					iLast, vecOpPositionIx[i] - iLast),
				false);
			m_vecChildren.push_back(pChild);

			iLast = vecOpPositionIx[i]+1;
		}

		ThresholdOpTreeNode * pChild = new ThresholdOpTreeNode();
		pChild->Parse(
			varreg,
			strExpression.substr(iLast),
			false);
		m_vecChildren.push_back(pChild);

	// Leaf node; evaluate as threshold operator
	} else {
		m_eType = ThresholdNode;
		m_threshop.Parse(varreg, strExpression);
	}
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
bool ThresholdOpTreeNode::IsSatisfied(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const int ix0,
	std::vector<int> * pvecRejectedCount
) const {
	if (pvecRejectedCount != NULL) {
		_ASSERT(pvecRejectedCount->size() == size());
	}
	if (m_eType == AlwaysTrueNode) {
		return true;
	}
	if (m_eType == ThresholdNode) {

		Variable & var = varreg.Get(m_threshop.m_varix);
		var.LoadGridData(varreg, vecFiles, grid);
		const DataArray1D<real> & dataState = var.GetData();

		bool fCurrentSatisfied = m_threshop.IsSatisfied<real>(grid, dataState, ix0);
		if (pvecRejectedCount != NULL) {
			(*pvecRejectedCount)[0] += (fCurrentSatisfied)?(0):(1);
		}
		return fCurrentSatisfied;
	}
	if (m_eType == OrNode) {
		_ASSERT(m_vecChildren.size() != 0);
		bool fAnySatisfied = false;
		for (int i = 0; i < m_vecChildren.size(); i++) {
			bool fCurrentSatisfied = m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, ix0);
			if (pvecRejectedCount != NULL) {
				(*pvecRejectedCount)[i] += (fCurrentSatisfied)?(0):(1);
			}
			fAnySatisfied = fAnySatisfied || fCurrentSatisfied;
		}
		return fAnySatisfied;
	}
	if (m_eType == AndNode) {
		_ASSERT(m_vecChildren.size() != 0);
		bool fAllSatisfied = true;
		for (int i = 0; i < m_vecChildren.size(); i++) {
			bool fCurrentSatisfied = m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, ix0);
			if (pvecRejectedCount != NULL) {
				(*pvecRejectedCount)[i] += (fCurrentSatisfied)?(0):(1);
			}
			fAllSatisfied = fAllSatisfied && fCurrentSatisfied;
		}
		return fAllSatisfied;
	}
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void ThresholdOpTreeNode::IsSatisfied(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	DataArray1D<int> & bTag
) const {
	_ASSERT(bTag.GetRows() == grid.GetSize());

	if (m_eType == AlwaysTrueNode) {
		return;
	}
	if (m_eType == ThresholdNode) {

		Variable & var = varreg.Get(m_threshop.m_varix);
		var.LoadGridData(varreg, vecFiles, grid);
		const DataArray1D<real> & dataState = var.GetData();

		for (int ix = 0; ix < bTag.GetRows(); ix++) {
			if (bTag[ix] == 0) {
				continue;
			}
			if (!m_threshop.IsSatisfied<real>(grid, dataState, ix)) {
				bTag[ix] = 0;
			}
		}
		return;
	}
	if (m_eType == OrNode) {
		_ASSERT(m_vecChildren.size() != 0);
		DataArray1D<int> bTagCopy(bTag.GetRows());
		for (int i = 0; i < m_vecChildren.size(); i++) {
			std::memcpy(&(bTagCopy[0]), &(bTag[0]), bTag.GetRows() * sizeof(int));
			for (int j = 0; j < bTagCopy.GetRows(); j++) {
				if (bTagCopy[j] == (-1)) {
					bTagCopy[j] = 1;
				}
			}
			m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, bTagCopy);
			for (int j = 0; j < bTag.GetRows(); j++) {
				if (bTagCopy[j] != 0) {
					bTag[j] = (-1);
				}
			}
		}
		for (int j = 0; j < bTag.GetRows(); j++) {
			if (bTag[j] == 1) {
				bTag[j] = 0;
			}
			if (bTag[j] == (-1)) {
				bTag[j] = 1;
			}
		}
		return;
	}
	if (m_eType == AndNode) {
		_ASSERT(m_vecChildren.size() != 0);
		for (int i = 0; i < m_vecChildren.size(); i++) {
			m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, bTag);
		}
		return;
	}
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

template <typename real>
void ThresholdOpTreeNode::IsSatisfied(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const std::set<int> & setCandidates,
	std::set<int> & setNewCandidates
) const {
	if (m_eType == AlwaysTrueNode) {
		setNewCandidates = setCandidates;
		return;
	}
	if (m_eType == ThresholdNode) {
		Variable & var = varreg.Get(m_threshop.m_varix);
		var.LoadGridData(varreg, vecFiles, grid);
		const DataArray1D<real> & dataState = var.GetData();

		for (auto it = setCandidates.begin(); it != setCandidates.end(); it++) {
			if (m_threshop.IsSatisfied<real>(grid, dataState, *it)) {
				setNewCandidates.insert(*it);
			}
		}
		_ASSERT(setNewCandidates.size() <= setCandidates.size());
		return;
	}
	if (m_eType == OrNode) {
		_ASSERT(m_vecChildren.size() != 0);
		for (int i = 0; i < m_vecChildren.size(); i++) {
			std::set<int> setNewCandidatesChild;

			m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, setCandidates, setNewCandidatesChild);

			for (auto it = setNewCandidatesChild.begin(); it != setNewCandidatesChild.end(); it++) {
				setNewCandidates.insert(*it);
			}
		}
		_ASSERT(setNewCandidates.size() <= setCandidates.size());
		return;
	}
	if (m_eType == AndNode) {
		_ASSERT(m_vecChildren.size() != 0);
		m_vecChildren[0]->IsSatisfied<real>(varreg, vecFiles, grid, setCandidates, setNewCandidates);
		for (int i = 1; i < m_vecChildren.size(); i++) {
			std::set<int> setNewCandidatesChild;
			for (auto it = setCandidates.begin(); it != setCandidates.end(); it++) {
				m_vecChildren[i]->IsSatisfied<real>(varreg, vecFiles, grid, setNewCandidates, setNewCandidatesChild);
			}
			setNewCandidates = setNewCandidatesChild;
		}
		_ASSERT(setNewCandidates.size() <= setCandidates.size());
		return;
	}
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////
// Explicit template instantiation

template bool ThresholdOp::IsSatisfied<float>(
	const SimpleGrid & grid,
	const DataArray1D<float> & dataState,
	const int ix0
) const;

template bool ThresholdOpTreeNode::IsSatisfied<float>(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const int ix0,
	std::vector<int> * pvecRejectedCount
) const;

template void ThresholdOpTreeNode::IsSatisfied<float>(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	DataArray1D<int> & bTag
) const;

template void ThresholdOpTreeNode::IsSatisfied<float>(
	VariableRegistry & varreg,
	NcFileVector & vecFiles,
	const SimpleGrid & grid,
	const std::set<int> & setCandidates,
	std::set<int> & setNewCandidates
) const;

///////////////////////////////////////////////////////////////////////////////

