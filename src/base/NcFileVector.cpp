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

#include "NetCDFUtilities.h"

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
	m_vecFileType.resize(0);
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
	m_vecFileTime.push_back(NcTimeDimension());

	// If a time index is already specified no need to read in the "time" variable
	if (lTimeIndex != InvalidTimeIndex) {
		m_vecFileType.push_back(FileType_Unknown);

	// Identify the FileType by the properties of the "time" variable
	} else {
		NcDim * dimTime = pNewFile->get_dim("time");
		if (dimTime == NULL) {
			NcVar * varTime = pNewFile->get_var("time");
			if (varTime != NULL) {
				ReadCFTimeDataFromNcFile(
					pNewFile,
					strFile,
					m_vecFileTime[m_vecFileTime.size()-1],
					false);
				m_vecFileType.push_back(FileType_Standard);
			} else {
				m_vecFileType.push_back(FileType_NoTimeDim);
			}
		} else {
			NcVar * varTime = pNewFile->get_var("time");
			if (varTime == NULL) {
				m_vecFileType.push_back(FileType_NoTimeVar);
			} else {
				// If the time variable exists read it in
				ReadCFTimeDataFromNcFile(
					pNewFile,
					strFile,
					m_vecFileTime[m_vecFileTime.size()-1],
					false);

				// Check time type
				NcAtt * attType = varTime->get_att("type");
				if (attType == NULL) {
					m_vecFileType.push_back(FileType_Standard);
				} else {
					std::string strType = attType->as_string(0);
					if (strType == "daily mean climatology") {
						m_vecFileType.push_back(FileType_DailyMeanClimo);
					} else if (strType == "annual mean climatology") {
						m_vecFileType.push_back(FileType_AnnualMeanClimo);
					} else {
						m_vecFileType.push_back(FileType_Standard);
						//_EXCEPTION2("Unrecognized time::type (%s) in file \"%s\"",
						//	strType.c_str(), strFile.c_str());
					}
				}
			}
		}
	}
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

long NcFileVector::GetTimeIx(size_t pos) const {
	_ASSERT(pos < m_vecTimeIxs.size());

	if (m_vecTimeIxs[pos] != InvalidTimeIndex) {
		if (m_vecTimeIxs[pos] == NoTimeIndex) {
			return NoTimeIndex;
		}

		NcDim * dimTime = m_vecNcFile[pos]->get_dim("time");
		if (dimTime == NULL) {
			_EXCEPTIONT("Logic error: Dimension \"time\" missing in file");
		}
		if ((m_vecTimeIxs[pos] < 0) || (m_vecTimeIxs[pos] >= dimTime->size())) {
			_EXCEPTIONT("Logic error: Time index out of range");
		}

		return m_vecTimeIxs[pos];

	} else {
		_ASSERT(m_vecNcFile.size() == m_vecFilenames.size());
		_ASSERT(m_vecNcFile.size() == m_vecFileType.size());

		const NcTimeDimension & vecTimes = m_vecFileTime[pos];
		_ASSERT(vecTimes.size() > 0);

		if (m_vecFileType[pos] == NcFileVector::FileType_Standard) {
			for (long t = 0; t < vecTimes.size(); t++) {
				if (vecTimes[t] == m_time) {
					return t;
				}
			}

		} else if (m_vecFileType[pos] == NcFileVector::FileType_DailyMeanClimo) {
			Time timeDailyMean(Time::CalendarNoLeap);
			timeDailyMean.SetYear(1);
			timeDailyMean.SetMonth(m_time.GetMonth());
			timeDailyMean.SetDay(m_time.GetDay());
			if (m_time.IsLeapDay()) {
				timeDailyMean.SetDay(28);
			}
			for (long t = 0; t < vecTimes.size(); t++) {
				if (vecTimes[t] == timeDailyMean) {
					return t;
				}
			}

		} else if (m_vecFileType[pos] == NcFileVector::FileType_AnnualMeanClimo) {
			if (vecTimes.size() != 1) {
				_EXCEPTIONT("Only one time expected in annual mean climatology");
			}
			return 0;
		}
	}
	_EXCEPTION2("Unable to identify time \"%s\" in \"%s\"",
		m_time.ToString().c_str(),
		m_vecFilenames[pos].c_str());
}

///////////////////////////////////////////////////////////////////////////////

