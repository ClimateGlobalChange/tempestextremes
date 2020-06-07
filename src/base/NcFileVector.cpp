///////////////////////////////////////////////////////////////////////////////
///
///	\file    NcFileVector.cpp
///	\author  Paul Ullrich
///	\version June 5, 2020
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

#include "NcFileVector.h"

///////////////////////////////////////////////////////////////////////////////

const long NcFileVector::InvalidTimeIndex = (-2);

const long NcFileVector::NoTimeIndex = (-1);

const size_t NcFileVector::InvalidIndex = (-1);

///////////////////////////////////////////////////////////////////////////////

void NcFileVector::clear() {
	for (size_t i = 0; i < size(); i++) {
		m_vecNcFile[i]->close();
		delete m_vecNcFile[i];
	}
	m_vecNcFile.resize(0);
	m_vecFilenames.resize(0);
	m_time = Time(Time::CalendarUnknown);
	m_vecTimeIxs.resize(0);
}

///////////////////////////////////////////////////////////////////////////////

void NcFileVector::InsertFile(
	const std::string & strFile,
	long lTimeIndex
) {
	NcFile * pNewFile = new NcFile(strFile.c_str());
	if (pNewFile == NULL) {
		_EXCEPTIONT("Unable to allocate new NcFile");
	}
	if (!pNewFile->is_valid()) {
		_EXCEPTION1("Cannot open input file \"%s\"", strFile.c_str());
	}
	m_vecNcFile.push_back(pNewFile);
	m_vecFilenames.push_back(strFile);
	m_vecTimeIxs.push_back(lTimeIndex);
}

///////////////////////////////////////////////////////////////////////////////

void NcFileVector::ParseFromString(
	const std::string & strFiles,
	bool fAppend
) {
	if (!fAppend) {
		m_vecNcFile.clear();
		m_vecFilenames.clear();
		m_vecTimeIxs.clear();
	}

	int iLast = 0;
	for (int i = 0; i <= strFiles.length(); i++) {
		if ((i == strFiles.length()) ||
		    (strFiles[i] == ';')
		) {
			std::string strFile =
				strFiles.substr(iLast, i - iLast);

			InsertFile(strFile);

			iLast = i+1;
		}
	}

	if (size() == 0) {
		_EXCEPTION1("No input files found in \"%s\"",
			strFiles.c_str());
	}
	_ASSERT(size() == m_vecFilenames.size());
}

///////////////////////////////////////////////////////////////////////////////

size_t NcFileVector::FindContainingVariable(
	const std::string & strVariable,
	NcVar ** pvar
) const {
	for (size_t s = 0; s < m_vecNcFile.size(); s++) {
		NcVar * var = m_vecNcFile[s]->get_var(strVariable.c_str());
		if (var != NULL) {
			if (pvar != NULL) {
				*pvar = var;
			}
			return s;
		}
	}
	if (pvar != NULL) {
		*pvar = NULL;
	}
	return InvalidIndex;
}

///////////////////////////////////////////////////////////////////////////////

