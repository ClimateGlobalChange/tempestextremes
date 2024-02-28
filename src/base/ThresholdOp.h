///////////////////////////////////////////////////////////////////////////////
///
///	\file    ThresholdOp.h
///	\author  Paul Ullrich
///	\version February 6, 2020
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

#ifndef _THRESHOLDOP_H_
#define _THRESHOLDOP_H_

#include "Announce.h"
#include "Variable.h"
#include "NcFileVector.h"
#include "SimpleGrid.h"
#include "DataArray1D.h"

#include <string>
#include <set>

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
		NotEqualTo,
		NoThreshold
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ThresholdOp() :
		m_varix(InvalidVariableIndex),
		m_eOp(GreaterThan),
		m_dThresholdValue(0.0),
		m_dDistanceDeg(0.0)
	{ }

public:
	///	<summary>
	///		Parse a threshold operator string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	);

public:
	///	<summary>
	///		Returns true if this threshold is satisfied, or false otherwise.
	///	</summary>
	template <typename real>
	bool IsSatisfied(
		const SimpleGrid & grid,
		const DataArray1D<real> & dataState,
		const int ix0
	) const;

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
	double m_dThresholdValue;

	///	<summary>
	///		Distance to search for threshold value (in degrees).
	///	</summary>
	double m_dDistanceDeg;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for holding a tree node from a threshold operator evaluation tree.
///	</summary>
class ThresholdOpTreeNode {

public:
	///	<summary>
	///		Types of threshold operators.
	///	</summary>
	enum NodeType {
		AlwaysTrueNode,
		ThresholdNode,
		AndNode,
		OrNode
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ThresholdOpTreeNode(
		NodeType eType = AlwaysTrueNode
	) :
		m_eType(eType)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~ThresholdOpTreeNode() {
		for (auto pchild : m_vecChildren) {
			delete pchild;
		}
	}

	///	<summary>
	///		Parse an expression and populate this node.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		std::string strExpression,
		bool fAllowZeroLengthExpression = true
	);

	///	<summary>
	///		Returns true if this threshold is satisfied.
	///	</summary>
	template <typename real>
	bool IsSatisfied(
		VariableRegistry & varreg,
		NcFileVector & vecFiles,
		const SimpleGrid & grid,
		const int ix0,
		std::vector<int> * pvecRejectedCount = NULL
	) const;

	///	<summary>
	///		Performs IsSatisfied in place over all values in a binary mask.
	///	</summary>
	template <typename real>
	void IsSatisfied(
		VariableRegistry & varreg,
		NcFileVector & vecFiles,
		const SimpleGrid & grid,
		DataArray1D<int> & bTag
	) const;

	///	<summary>
	///		Performs IsSatisfied over a set of candidate nodes, returning
	///		the set of candidates that satisfy the criteria.
	///	</summary>
	template <typename real>
	void IsSatisfied(
		VariableRegistry & varreg,
		NcFileVector & vecFiles,
		const SimpleGrid & grid,
		const std::set<int> & setCandidates,
		std::set<int> & setNewCandidates
	) const;

	///	<summary>
	///		If AlwaysTrueNode return 0.
	///		If ThresholdNode return 1.
	///		If AndNode or OrNode return the number of children.
	///	</summary>
	size_t size() const {
		if (m_eType == AlwaysTrueNode) {
			return 0;
		}
		if (m_eType == ThresholdNode) {
			return 1;
		}
		return m_vecChildren.size();
	}

	///	<summary>
	///		Get the name of this ThresholdOpTreeNode.
	///	</summary>
	std::string ToString(
		VariableRegistry & varreg,
		int ix
	) const {
		_ASSERT((ix >= 0) && (ix < m_vecChildren.size()));
		if (m_eType == AlwaysTrueNode) {
			return 0;
		}
		if (m_eType == ThresholdNode) {
			return varreg.GetVariableString(m_threshop.m_varix);
		}
		if (m_vecChildren[ix]->m_eType == ThresholdNode) {
			return m_vecChildren[ix]->ToString(varreg, 0);
		}
		return std::string("Compound Expr. ") + std::to_string(ix);
	}

public:
	///	<summary>
	///		Node type.
	///	</summary>
	NodeType m_eType;

	///	<summary>
	///		Threshold operator at this node.
	///	</summary>
	ThresholdOp m_threshop;

	///	<summary>
	///		Children nodes.
	///	</summary>
	std::vector<ThresholdOpTreeNode *> m_vecChildren;
};

///////////////////////////////////////////////////////////////////////////////

#endif // _THRESHOLDOP_H_

