///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileEditor.cpp
///	\author  Paul Ullrich
///	\version August 30, 2018
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

#ifndef _NODEFILEUTILITIES_H_
#define _NODEFILEUTILITIES_H_

#include "Exception.h"
#include "TimeObj.h"
#include "STLStringHelper.h"

#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

class SimpleGrid;
class AutoCurator;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of PathNode indices (path #, pathnode #).
///	</summary>
typedef std::vector< std::pair<int, int> > PathNodeIndexVector;

///	<summary>
///		A map from Times to file lines, used for StitchNodes formatted output.
///	</summary>
typedef std::map<Time, PathNodeIndexVector> TimeToPathNodeMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class containing a list of strings for each column data header.
///	</summary>
class ColumnDataHeader {

public:
	///	<summary>
	///		Invalid index.
	///	</summary>
	static const int InvalidIndex = (-1);

public:
	///	<summary>
	///		Assignment operator overload.
	///	</summary>
	ColumnDataHeader & operator=(const ColumnDataHeader & cdh) {
		m_vecColumnHeader = cdh.m_vecColumnHeader;
		return *this;
	}

	///	<summary>
	///		Parse the format string for the column data header.
	///	</summary>
	void Parse(
		const std::string & strFormat
	) {
		m_vecColumnHeader.clear();
		if (strFormat.length() == 0) {
			return;
		}

		int iLast = 0;
		for (int i = 0; i < strFormat.length(); i++) {
			if (strFormat[i] == ',') {
				m_vecColumnHeader.push_back(
					strFormat.substr(iLast, i-iLast));
				iLast = i+1;
			}
		}
		m_vecColumnHeader.push_back(
			strFormat.substr(iLast, strFormat.length()-iLast));
	}

	///	<summary>
	///		Size of the header.
	///	</summary>
	size_t size() const {
		return m_vecColumnHeader.size();
	}

	///	<summary>
	///		Get the specified index of the column data header.
	///	</summary>
	const std::string & operator[](const int & ix) const {
		return m_vecColumnHeader[ix];
	}

	///	<summary>
	///		Get the specified index of the column data header.
	///	</summary>
	std::string & operator[](const int & ix) {
		return m_vecColumnHeader[ix];
	}

	///	<summary>
	///		Find the specified string index from a string.
	///	</summary>
	int GetIndexFromString(const std::string & str) const {
		for (int i = 0; i < m_vecColumnHeader.size(); i++) {
			if (m_vecColumnHeader[i] == str) {
				return i;
			}
		}
		return InvalidIndex;
	}

	///	<summary>
	///		Populate a vector of indices from a ColumnDataHeader.
	///	</summary>
	void GetIndicesFromColumnDataHeader(
		const ColumnDataHeader & cdh,
		std::vector<int> & vecColumnDataOutIx
	) const {
		for (int i = 0; i < cdh.size(); i++) {
			int ix = GetIndexFromString(cdh[i]);
			if (ix == (-1)) {
				_EXCEPTION1("Unknown column data header \"%s\"",
					cdh[i].c_str());
			} else {
				vecColumnDataOutIx.push_back(ix);
			}
		}
	}

	///	<summary>
	///		Clear the contents of this ColumnDataHeader.
	///	</summary>
	void clear() {
		m_vecColumnHeader.clear();
	}

	///	<summary>
	///		Add a new header string.
	///	</summary>
	void push_back(const std::string & str) {
		m_vecColumnHeader.push_back(str);
	}

protected:
	///	<summary>
	///		Column headers.
	///	</summary>
	std::vector<std::string> m_vecColumnHeader;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnData {

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~ColumnData() {
	}

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		return std::string("null");
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnData();
	}
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataString : public ColumnData {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	ColumnDataString() :
		m_strData()
	{ }

	///	<summary>
	///		Constructor with string initializer.
	///	</summary>
	ColumnDataString(
		const std::string & strData
	) {
		m_strData = strData;
	}

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		return m_strData;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnDataString(*this);
	}

public:
	///	<summary>
	///		The string object in this column data.
	///	</summary>
	std::string m_strData;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataDouble : public ColumnData {

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	ColumnDataDouble()
	{ }

	///	<summary>
	///		Constructor with double initializer.
	///	</summary>
	ColumnDataDouble(
		const double & dData
	) {
		m_dData = dData;
	}

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char szBuffer[32];
		snprintf(szBuffer, 32, "%3.6e", m_dData);
		return std::string(szBuffer);
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnDataDouble(*this);
	}

public:
	///	<summary>
	///		The double in this ColumnData.
	///	</summary>
	double m_dData;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataRLLVelocity : public ColumnData {

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char szVelocity[100];
		snprintf(szVelocity, 100, "\"%3.6e %3.6e\"", m_dU, m_dV);
		return szVelocity;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnDataRLLVelocity(*this);
	}

public:
	///	<summary>
	///		Zonal velocity.
	///	</summary>
	double m_dU;

	///	<summary>
	///		Meridional velocity.
	///	</summary>
	double m_dV;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataDoubleArrayTemplate : public ColumnData {
public:
	///	<summary>
	///		Get the indices of this array.
	///	</summary>
	virtual const std::vector<double> & GetIndices() const = 0;

	///	<summary>
	///		Get the values of this array.
	///	</summary>
	virtual const std::vector<double> & GetValues() const = 0;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnData1DArray : public ColumnData {
public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char buf[100];
		std::string strOut = "\"[";
		for (int i = 0; i < m_dValues.size(); i++) {
			if (i != 0) {
				strOut += ",";
			}
			snprintf(buf, 100, "%3.6e", m_dValues[i]);
			strOut += buf;
		}
		strOut += "]\"";
		return strOut;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnData1DArray(*this);
	}

public:
	///	<summary>
	///		Vector of vectors.
	///	</summary>
	std::vector<double> m_dValues;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnData2DArray : public ColumnData {
public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char buf[100];
		std::string strOut = "\"[";
		for (int i = 0; i < m_dValues.size(); i++) {
			if (i != 0) {
				strOut += ",";
			}
			strOut += "[";
			for (int j = 0; j < m_dValues[i].size(); j++) {
				if (j != 0) {
					strOut += ",";
				}
				snprintf(buf, 100, "%3.6e", m_dValues[i][j]);
				strOut += buf;
			}
			strOut += "]";
		}
		strOut += "]\"";
		return strOut;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnData2DArray(*this);
	}

public:
	///	<summary>
	///		Vector of vectors.
	///	</summary>
	std::vector< std::vector<double> > m_dValues;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataRadialProfile : public ColumnDataDoubleArrayTemplate {

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char buf[100];
		std::string strOut = "\"[";
		for (int i = 0; i < m_dValues.size(); i++) {
			snprintf(buf, 100, "%3.6e", m_dValues[i]);
			strOut += buf;
			if (i == m_dValues.size()-1) {
				strOut += "]\"";
			} else {
				strOut += ",";
			}
		}
		return strOut;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnDataRadialProfile(*this);
	}

	///	<summary>
	///		Get the indices of this array.
	///	</summary>
	virtual const std::vector<double> & GetIndices() const {
		return m_dR;
	}

	///	<summary>
	///		Get the values of this array.
	///	</summary>
	virtual const std::vector<double> & GetValues() const {
		return m_dValues;
	}

public:
	///	<summary>
	///		Vector of radial coordinates.
	///	</summary>
	std::vector<double> m_dR;

	///	<summary>
	///		Vector of values.
	///	</summary>
	std::vector<double> m_dValues;
};

///////////////////////////////////////////////////////////////////////////////

class ColumnDataRadialVelocityProfile : public ColumnDataDoubleArrayTemplate {

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char buf[100];
		std::string strOut = "\"[";
		for (int i = 0; i < m_dUa.size(); i++) {
			snprintf(buf, 100, "%3.6e", m_dUa[i]);
			strOut += buf;
			if (i == m_dUa.size()-1) {
				strOut += "]\"";
			} else {
				strOut += ",";
			}
		}
		return strOut;
	}

	///	<summary>
	///		Return a duplicate of this ColumnData.
	///	</summary>
	virtual ColumnData * Duplicate() const {
		return new ColumnDataRadialVelocityProfile(*this);
	}

	///	<summary>
	///		Get the indices of this array.
	///	</summary>
	virtual const std::vector<double> & GetIndices() const {
		return m_dR;
	}

	///	<summary>
	///		Get the values of this array.
	///	</summary>
	virtual const std::vector<double> & GetValues() const {
		return m_dUa;
	}

public:
	///	<summary>
	///		Vector of radial coordinates.
	///	</summary>
	std::vector<double> m_dR;

	///	<summary>
	///		Vector of azimuthal velocities.
	///	</summary>
	std::vector<double> m_dUa;

	///	<summary>
	///		Vector of radial velocities.
	///	</summary>
	std::vector<double> m_dUr;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A struct for expressing one node along a path.
///	</summary>
class PathNode {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PathNode() :
		m_gridix(0),
		m_fileix(-1),
		m_time()
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	PathNode(const PathNode & node) {
		for (int i = 0; i < node.m_vecColumnData.size(); i++) {
			m_vecColumnData.push_back(node.m_vecColumnData[i]->Duplicate());
		}
		m_gridix = node.m_gridix;
		m_fileix = node.m_fileix;
		m_time = node.m_time;
	}

	///	<summary>
	///		Clear column data.
	///	</summary>
	void clear() {
		for (int i = 0; i < m_vecColumnData.size(); i++) {
			delete m_vecColumnData[i];
		}
		m_vecColumnData.clear();
	}

	///	<summary>
	///		Destructor.
	///	</summary>
	~PathNode() {
		clear();
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	PathNode & operator=(const PathNode & node) {
		clear();
		for (int i = 0; i < node.m_vecColumnData.size(); i++) {
			m_vecColumnData.push_back(node.m_vecColumnData[i]->Duplicate());
		}
		m_gridix = node.m_gridix;
		m_fileix = node.m_fileix;
		m_time = node.m_time;
		return (*this);
	}

public:
	///	<summary>
	///		Push back a ColumnData object into column.
	///	</summary>
	void PushColumnData(
		ColumnData * pcdat
	) {
		m_vecColumnData.push_back(pcdat);
	}

	///	<summary>
	///		Push back a string object into column.
	///	</summary>
	void PushColumnDataString(
		const std::string & strString
	) {
		m_vecColumnData.push_back(
			new ColumnDataString(strString));
	}

	///	<summary>
	///		Duplicate data from given index.
	///	</summary>
	void Duplicate(
		int ix
	) {
		if ((ix < 0) || (ix >= m_vecColumnData.size())) {
			_EXCEPTIONT("Index out of range");
		}
		m_vecColumnData.push_back(
			m_vecColumnData[ix]->Duplicate());
	}

	///	<summary>
	///		Get column data as integer.
	///	</summary>
	int GetColumnDataAsInteger(
		int ix
	) const {
		if ((ix < 0) || (ix >= m_vecColumnData.size())) {
			_EXCEPTION();
		}
		std::string str = m_vecColumnData[ix]->ToString();
		if (!STLStringHelper::IsInteger(str)) {
			_EXCEPTION1("Column data \"%s\" cannot be cast to type integer",
				str.c_str());
		}
		return atoi(str.c_str());
	}

	///	<summary>
	///		Get column data as integer.
	///	</summary>
	int GetColumnDataAsInteger(
		const ColumnDataHeader & cdh,
		const std::string & strColumn
	) const {
		if (STLStringHelper::IsInteger(strColumn)) {
			return atoi(strColumn.c_str());
		}

		int ix = cdh.GetIndexFromString(strColumn);
		if (ix == (-1)) {
			_EXCEPTION1("Invalid column header \"%s\"", strColumn.c_str());
		}
		return GetColumnDataAsInteger(ix);
	}

	///	<summary>
	///		Get column data as double (using a column index).
	///	</summary>
	double GetColumnDataAsDouble(
		int ix
	) const {
		if ((ix < 0) || (ix >= m_vecColumnData.size())) {
			_EXCEPTION();
		}
		std::string str = m_vecColumnData[ix]->ToString();
		if (!STLStringHelper::IsFloat(str)) {
			_EXCEPTION1("Column data \"%s\" cannot be cast to type double",
				str.c_str());
		}
		return atof(str.c_str());
	}

	///	<summary>
	///		Get column data as double (using a column name).
	///	</summary>
	double GetColumnDataAsDouble(
		const ColumnDataHeader & cdh,
		const std::string & strColumn
	) const {
		if (STLStringHelper::IsFloat(strColumn)) {
			return atof(strColumn.c_str());
		}

		int ix = cdh.GetIndexFromString(strColumn);
		if (ix == (-1)) {
			_EXCEPTION1("Invalid column header \"%s\"", strColumn.c_str());
		}
		return GetColumnDataAsDouble(ix);
	}

public:
	///	<summary>
	///		Index on the grid of this point.
	///	</summary>
	long m_gridix;

	///	<summary>
	///		Index of this PathNode in the order it appears in the file.
	///	</summary>
	size_t m_fileix;

	///	<summary>
	///		Time associated with this node.
	///	</summary>
	Time m_time;

	///	<summary>
	///		Array of column data associated with this node and time.
	///	</summary>
	std::vector<ColumnData *> m_vecColumnData;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A struct for expressing a path.
///	</summary>
class Path : public std::vector<PathNode> {
public:
	///	<summary>
	///		Duplicate the given data column across all PathNodes.
	///	</summary>
	void Duplicate(
		int ix
	) {
		for (int j = 0; j < size(); j++) {
			(*this)[j].Duplicate(ix);
		}
	}

public:
	///	<summary>
	///		Start time of the path.
	///	</summary>
	Time m_timeStart;

	///	<summary>
	///		End time of the path.
	///	</summary>
	Time m_timeEnd;
};

///////////////////////////////////////////////////////////////////////////////

class PathVector : public std::vector<Path> {

public:
	///	<summary>
	///		Duplicate the given data column across all Paths
	///	</summary>
	void Duplicate(
		int ix
	) {
		for (int p = 0; p < size(); p++) {
			(*this)[p].Duplicate(ix);
		}
	}

};

///////////////////////////////////////////////////////////////////////////////

class NodeFile {

public:
	///	<summary>
	///		Different types of NodeFiles.
	///	</summary>
	enum PathType {
		PathTypeDN,
		PathTypeSN
	};

	///	<summary>
	///		Format of the NodeFile.
	///	</summary>
	enum FileFormat {
		FileFormatGFDL,
		FileFormatCSV
	};

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	NodeFile() :
		m_ePathType(PathTypeSN),
		m_fContainsNegativeGridIx(false)
	{ }

public:
	///	<summary>
	///		Read in a nodefile as a CSV and parse it into a PathVector.
	///	</summary>
	void ReadCSV(
		const std::string & strNodeFile,
		PathType ePathType,
		ColumnDataHeader & cdh,
		const std::vector<size_t> & nGridDim,
		size_t sGridSize,
		Time::CalendarType caltype
	);

	///	<summary>
	///		Read in a node file and parse it into a PathVector.
	///	</summary>
	void Read(
		const std::string & strNodeFile,
		PathType ePathType,
		ColumnDataHeader & cdh,
		const std::vector<size_t> & nGridDim,
		size_t sGridSize,
		Time::CalendarType caltype
	);

	///	<summary>
	///		Read in a node file and parse it into a PathVector.
	///	</summary>
	void Read(
		const std::string & strNodeFile,
		PathType ePathType,
		ColumnDataHeader & cdh,
		const SimpleGrid & grid,
		Time::CalendarType caltype
	);

	///	<summary>
	///		Write a node file.
	///	</summary>
	void Write(
		const std::string & strNodeFile,
		const SimpleGrid * pgrid = NULL,
		const std::vector<int> * pvecColumnDataOutIx = NULL,
		FileFormat eFileFormat = FileFormatGFDL,
		bool fIncludeHeader = false,
		bool fOutputSeconds = false
	);

	///	<summary>
	///		Generate the TimeToPathNodeMap.
	///	</summary>
	void GenerateTimeToPathNodeMap();

public:
	///	<summary>
	///		Apply a time delta to the node file.
	///	</summary>
	void ApplyTimeDelta(
		const Time & timeDelta,
		bool fAddTimeDelta = true
	);

public:
	///	<summary>
	///		Initialize the interpolation indices at the given time.
	///	</summary>
	void Interpolate(
		const Time & time
	);

	///	<summary>
	///		Get the array of path ids from interpolation.
	///	</summary>
	const std::vector<size_t> & GetInterpolatedPathIds() const {
		return m_vecInterpPathId;
	}

	///	<summary>
	///		Get the interpolated coordinates from nodefile.
	///	</summary>
	void InterpolatedNodeCoordinatesRad(
		const std::string & strLonName,
		const std::string & strLatName,
		std::vector<double> & vecInterpLonRad,
		std::vector<double> & vecInterpLatRad
	) const;

	///	<summary>
	///		Get the interpolated column value.
	///	</summary>
	void InterpolatedColumnDouble(
		const std::string & strHeaderName,
		std::vector<double> & vecInterpDouble
	) const;

public:
	///	<summary>
	///		Get the PathType associated with this NodeFile.
	///	</summary>
	const PathType & GetPathType() const {
		return m_ePathType;
	}

	///	<summary>
	///		Get the ColumnDataHeader.
	///	</summary>
	ColumnDataHeader & GetColumnDataHeader() {
		return m_cdh;
	}

	///	<summary>
	///		Get the PathVector.
	///	</summary>
	PathVector & GetPathVector() {
		return m_pathvec;
	}

	///	<summary>
	///		Get the TimeToPathNodeMap.
	///	</summary>
	TimeToPathNodeMap & GetTimeToPathNodeMap() {
		return m_mapTimeToPathNode;
	}

	///	<summary>
	///		Return true if the NodeFile may contain negative grid indices.
	///	</summary>
	bool ContainsNegativeGridIx() const {
		return m_fContainsNegativeGridIx;
	}

public:
	///	<summary>
	///		The type of path described by this NodeFile.
	///	</summary>
	PathType m_ePathType;

	///	<summary>
	///		The ColumnDataHeaders of this file.
	///	</summary>
	ColumnDataHeader m_cdh;

	///	<summary>
	///		Vector of paths containing the ColumnData from the NodeFile.
	///	</summary>
	PathVector m_pathvec;

	///	<summary>
	///		A map from Times to PathNodes.
	///	</summary>
	///	<remarks>
	///		Because nodes in StitchNodes format output are not ordered in
	///		time, efficient data I/O requires us to reorganize the input
	///		lines by time.
	///	</remarks>
	TimeToPathNodeMap m_mapTimeToPathNode;

protected:
	///	<summary>
	///		Vector of paths from interpolation.
	///	</summary>
	std::vector<size_t> m_vecInterpPathId;

	///	<summary>
	///		Vector of time indices from interpolation.
	///	</summary>
	std::vector<size_t> m_vecInterpTimeId;

	///	<summary>
	///		Vector of alpha indices from interpolation.
	///	</summary>
	std::vector<double> m_vecInterpAlpha;

protected:
	///	<summary>
	///		A flag indicating the NodeFile may contain negative grid indices.
	///	</summary>
	bool m_fContainsNegativeGridIx;
};

///////////////////////////////////////////////////////////////////////////////

#endif //_NODEFILEUTILITIES_H_

