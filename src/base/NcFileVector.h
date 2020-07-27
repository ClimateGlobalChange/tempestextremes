///////////////////////////////////////////////////////////////////////////////
///
///	\file    NcFileVector.h
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

#ifndef _NCFILEVECTOR_H_
#define _NCFILEVECTOR_H_

#include "Exception.h"
#include "netcdfcpp.h"
#include "TimeObj.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of open NetCDF files.
///	</summary>
class NcFileVector {

public:
	///	<summary>
	///		Invalid time index.
	///	</summary>
	static const long InvalidTimeIndex;

	///	<summary>
	///		No time index.
	///	</summary>
	static const long NoTimeIndex;

	///	<summary>
	///		Invalid position in array.
	///	</summary>
	static const size_t InvalidIndex;

	///	<summary>
	///		File type.
	///	</summary>
	enum FileType {
		FileType_Unknown = (-1),
		FileType_Standard = 0,
		FileType_NoTime = 1,
		FileType_DailyMeanClimo = 2,
	};

private:
	///	<summary>
	///		Copy constructor.
	///	</summary>
	NcFileVector(const NcFileVector & vecFiles) {
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	NcFileVector & operator= (const NcFileVector & vecFiles) {
		return (*this);
	}

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	NcFileVector() :
		m_time(Time::CalendarUnknown)
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~NcFileVector() {
		NcFileVector::clear();
	}

	///	<summary>
	///		Size of this NcFileVector.
	///	</summary>
	size_t size() const {
		return m_vecNcFile.size();
	}

	///	<summary>
	///		Clear the contents of this NcFileVector.
	///	</summary>
	void clear();

	///	<summary>
	///		Insert a single filename string into the NcFileVector.
	///	</summary>
	void InsertFile(
		const std::string & strFile,
		long lTimeIndex = InvalidTimeIndex
	);

	///	<summary>
	///		Parse from a semi-colon delineated string of file names.
	///	</summary>
	void ParseFromString(
		const std::string & strFiles,
		bool fAppend = true
	);

	///	<summary>
	///		Get an iterator to the file containing the specified variable.
	///	</summary>
	size_t FindContainingVariable(
		const std::string & strVariable,
		NcVar ** pvar
	) const;

	///	<summary>
	///		Get the filename at the specified position.
	///	</summary>
	const std::string & GetFilename(size_t pos) const {
		_ASSERT(pos < m_vecFilenames.size());
		return m_vecFilenames[pos];
	}

	///	<summary>
	///		Get the FileType at the specified position.
	///	</summary>
	const FileType & GetFileType(size_t pos) const {
		_ASSERT(pos < m_vecFilenames.size());
		return m_vecFileType[pos];
	}

public:
	///	<summary>
	///		Get the Time used to construct this NcFileVector.
	///	</summary>
	const Time & GetTime() const {
		return m_time;
	}

	///	<summary>
	///		Set the Time for this NcFileVector.
	///	</summary>
	void SetTime(const Time & time) {
		m_time = time;
	}

	///	<summary>
	///		Get the time index from the specified file corresponding to m_time.
	///	</summary>
	long GetTimeIx(size_t pos) const;

	///	<summary>
	///		Set the time index across all files.
	///	</summary>
	void SetConstantTimeIx(long lTime) {
		m_time = Time(Time::CalendarNone);
		m_time.SetYear(lTime);
		for (size_t f = 0; f < m_vecTimeIxs.size(); f++) {
			m_vecTimeIxs[f] = lTime;
		}
	}

public:
	///	<summary>
	///		Square bracket accessor.
	///	</summary>
	NcFile * operator[](size_t pos) const {
		return m_vecNcFile.at(pos);
	}

protected:
	///	<summary>
	///		Vector of NcFile* pointers.
	///	</summary>
	std::vector<NcFile *> m_vecNcFile;

	///	<summary>
	///		Vector of file names.
	///	</summary>
	std::vector<std::string> m_vecFilenames;

	///	<summary>
	///		Type of file.
	///	</summary>
	std::vector<FileType> m_vecFileType;

protected:
	///	<summary>
	///		Time associated with the time indices.
	///	</summary>
	Time m_time;

	///	<summary>
	///		Vector of time indices in each file with associated Time m_time.
	///	</summary>
	std::vector<long> m_vecTimeIxs;
};

///////////////////////////////////////////////////////////////////////////////

#endif

