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
		return (-1);
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
		sprintf(szBuffer, "%3.6e", m_dData);
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
		sprintf(szVelocity, "\"%3.6e %3.6e\"", m_dU, m_dV);
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

class ColumnDataRadialProfile : public ColumnDataDoubleArrayTemplate {

public:
	///	<summary>
	///		Express this column data as a string.
	///	</summary>
	virtual std::string ToString() const {
		char buf[100];
		std::string strOut = "\"[";
		for (int i = 0; i < m_dValues.size(); i++) {
			sprintf(buf, "%3.6e", m_dValues[i]);
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
			sprintf(buf, "%3.6e", m_dUa[i]);
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
		m_time()
	{ }

	///	<summary>
	///		Destructor.
	///	</summary>
	~PathNode() {
		for (int i = 0; i < m_vecColumnData.size(); i++) {
			delete m_vecColumnData[i];
		}
	}

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
		const ColumnDataHeader & cdh,
		const std::string & strColumn
	) {
		if (STLStringHelper::IsInteger(strColumn)) {
			return atoi(strColumn.c_str());
		}

		int ix = cdh.GetIndexFromString(strColumn);
		if (ix == (-1)) {
			_EXCEPTION1("Invalid column header \"%s\"", strColumn.c_str());
		}
		std::string str = m_vecColumnData[ix]->ToString();
		if (!STLStringHelper::IsInteger(str)) {
			_EXCEPTIONT("Column header \"%s\" cannot be cast to type integer");
		}
		return atoi(str.c_str());
	}

	///	<summary>
	///		Get column data as double.
	///	</summary>
	double GetColumnDataAsDouble(
		const ColumnDataHeader & cdh,
		const std::string & strColumn
	) {
		if (STLStringHelper::IsFloat(strColumn)) {
			return atof(strColumn.c_str());
		}

		int ix = cdh.GetIndexFromString(strColumn);
		if (ix == (-1)) {
			_EXCEPTION1("Invalid column header \"%s\"", strColumn.c_str());
		}
		if ((ix < 0) || (ix >= m_vecColumnData.size())) {
			_EXCEPTION();
		}
		std::string str = m_vecColumnData[ix]->ToString();
		if (!STLStringHelper::IsFloat(str)) {
			_EXCEPTION1("Column header \"%s\" cannot be cast to type double",
				strColumn.c_str());
		}
		return atof(str.c_str());
	}

public:
	///	<summary>
	///		Index on the grid of this point.
	///	</summary>
	size_t m_gridix;

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
class Path {
public:
	///	<summary>
	///		Duplicate the given data column across all PathNodes.
	///	</summary>
	void Duplicate(
		int ix
	) {
		for (int j = 0; j < m_vecPathNodes.size(); j++) {
			m_vecPathNodes[j].Duplicate(ix);
		}
	}

	///	<summary>
	///		Array overload to access PathNode.
	///	</summary>
	PathNode & operator[](const int & ix) {
		return m_vecPathNodes[ix];
	}

	///	<summary>
	///		Array overload to access PathNode.
	///	</summary>
	const PathNode & operator[](const int & ix) const {
		return m_vecPathNodes[ix];
	}

public:
	///	<summary>
	///		Start time of the path.
	///	</summary>
	Time m_timeStart;

	///	<summary>
	///		Array of PathNodes in the path, in temporal order.
	///	</summary>
	std::vector<PathNode> m_vecPathNodes;
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
		PathTypeDCU,
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
	///		Read in a node file and parse it into a PathVector.
	///	</summary>
	void Read(
		const std::string & strNodeFile,
		PathType ePathType,
		const ColumnDataHeader & cdh,
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
		bool fIncludeHeader = false
	);

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
	TimeToPathNodeMap m_mapTimeToPathNode;
};

///////////////////////////////////////////////////////////////////////////////

#endif //_NODEFILEUTILITIES_H_

