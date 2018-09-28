///////////////////////////////////////////////////////////////////////////////
///
///	\file    ArgumentTree.h
///	\author  Paul Ullrich
///	\version September 7, 2018
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

#ifndef _ARGUMENTTREE_H_
#define _ARGUMENTTREE_H_

#include "Exception.h"

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

class ArgumentTree {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ArgumentTree(
		bool fSemiColonDelimited
	) :
		m_fSemiColonDelimited(fSemiColonDelimited)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~ArgumentTree();

private:
	///	<summary>
	///		Parse an argument substring (utility function).
	///	</summary>
	void SubParse(
		const std::string & str,
		char szListType
	);

public:
	///	<summary>
	///		Parse an argument string.
	///	</summary>
	bool Parse(
		const std::string & str
	);

public:
	///	<summary>
	///		Determine if this is semicolon delimited.
	///	</summary>
	bool IsSemiColonDelimited() const {
		return m_fSemiColonDelimited;
	}

	///	<summary>
	///		Get the argument count.
	///	</summary>
	size_t size() const {
		return m_strArgumentList.size();
	}

	///	<summary>
	///		Get the specified argument string.
	///	</summary>
	const std::string & GetArgumentString(size_t s) const {
		if (s >= size()) {
			_EXCEPTION1("Argument %lu out of range",
				static_cast<unsigned long>(s));
		}
		return m_strArgumentList[s];
	}

	///	<summary>
	///		Array operator to access argument string.
	///	</summary>
	const std::string & operator[](size_t s) const {
		return GetArgumentString(s);
	}

	///	<summary>
	///		Get the specified sub-argument tree.
	///	</summary>
	const ArgumentTree * GetSubTree(size_t s) const {
		if (s >= size()) {
			_EXCEPTION1("Argument %lu out of range",
				static_cast<unsigned long>(s));
		}
		return m_vecSubArguments[s];
	}

	///	<summary>
	///		Get the specified sub-argument tree.
	///	</summary>
	ArgumentTree * GetSubTree(size_t s) {
		if (s >= size()) {
			_EXCEPTION1("Argument %lu out of range",
				static_cast<unsigned long>(s));
		}
		return m_vecSubArguments[s];
	}

private:
	///	<summary>
	///		True if this ArgumentTree is semi-colon delimited.
	///	</summary>
	bool m_fSemiColonDelimited;

	///	<summary>
	///		Type of list.
	///	</summary>
	char m_szListType;

	///	<summary>
	///		Strings containing arguments.
	///	</summary>
	std::vector<std::string> m_strArgumentList;

	///	<summary>
	///		Pointers to sub-arguments.
	///	</summary>
	std::vector<ArgumentTree *> m_vecSubArguments;
};

///////////////////////////////////////////////////////////////////////////////

#endif //_ARGUMENTTREE_H_

