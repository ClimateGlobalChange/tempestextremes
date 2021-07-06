///////////////////////////////////////////////////////////////////////////////
///
///	\file    NodeFileUtilities.cpp
///	\author  Paul Ullrich
///	\version October 4, 2018
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

#include "NodeFileUtilities.h"
#include "Defines.h"

#include <iostream>

#include "SimpleGrid.h"
#include "AutoCurator.h"
#include "CoordTransforms.h"

///////////////////////////////////////////////////////////////////////////////

void NodeFile::ReadCSV(
	const std::string & strNodeFile,
	PathType ePathType,
	ColumnDataHeader & cdh,
	const std::vector<size_t> & nGridDim,
	size_t sGridSize,
	Time::CalendarType caltype
) {
	// Not implemented
	if (ePathType == PathTypeDN) {
		_EXCEPTIONT("Not implemented");
	}

	// Check that this NodeFile is uninitialized
	if ((m_cdh.size() != 0) || (m_pathvec.size() != 0)) {
		_EXCEPTIONT("Attempting to ReadCSV() on an initialized NodeFile");
	}

	// Store the PathType
	m_ePathType = ePathType;

	// String buffer
	std::string strBuffer;

	// Only support grids of dimension 1 or 2
	if ((nGridDim.size() < 1) || (nGridDim.size() > 2)) {
		_EXCEPTIONT("Grid dimension out of range:  Only grids of dimension 1 or 2 supported");
	}

	// Clear the TimeToPathNodeMap
	m_mapTimeToPathNode.clear();

	// Open the file as an input stream
	std::ifstream ifInput(strNodeFile);
	if (!ifInput.is_open()) {
		_EXCEPTION1("Unable to open nodefile \"%s\"", strNodeFile.c_str());
	}

	// Attempt to read column data header from file
	if ((cdh.size() == 1) && (cdh[0] == "(auto)")) {
		cdh.clear();

		getline(ifInput, strBuffer);
		if (ifInput.eof()) {
			_EXCEPTION1("Input nodefile \"%s\" contains no data", strNodeFile.c_str());
		}

		int iLast = 0;
		for (int i = 0; i <= strBuffer.length(); i++) {
			if ((i == strBuffer.length()) || (strBuffer[i] == ',')) {
				std::string strHeader = strBuffer.substr(iLast, i-iLast);
				STLStringHelper::RemoveWhitespaceInPlace(strHeader);
				iLast = i+1;
				cdh.push_back(strHeader);
				std::cout << strHeader << std::endl;
			}
		}
	}

	// Store CDH
	m_cdh = cdh;

	// Check for all necessary headers
	int iColTrackId = (-1);
	int iColYear = (-1);
	int iColMonth = (-1);
	int iColDay = (-1);
	int iColHour = (-1);
	int iColLon = (-1);
	int iColLat = (-1);

	{
		for (int iCol = 0; iCol < m_cdh.size(); iCol++) {
			if (m_cdh[iCol] == "track_id") {
				iColTrackId = iCol;
			}
			if (m_cdh[iCol] == "year") {
				iColYear = iCol;
			}
			if (m_cdh[iCol] == "month") {
				iColMonth = iCol;
			}
			if (m_cdh[iCol] == "day") {
				iColDay = iCol;
			}
			if (m_cdh[iCol] == "hour") {
				iColHour = iCol;
			}
			if ((m_cdh[iCol] == "lon") || (m_cdh[iCol] == "longitude")) {
				iColLon = iCol;
			}
			if ((m_cdh[iCol] == "lat") || (m_cdh[iCol] == "latitude")) {
				iColLat = iCol;
			}
		}

		if (iColTrackId == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"track_id\"", strNodeFile.c_str());
		}
		if (iColYear == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"year\"", strNodeFile.c_str());
		}
		if (iColMonth == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"month\"", strNodeFile.c_str());
		}
		if (iColDay == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"day\"", strNodeFile.c_str());
		}
		if (iColHour == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"hour\"", strNodeFile.c_str());
		}
		if (iColLon == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"lon\"", strNodeFile.c_str());
		}
		if (iColLat == (-1)) {
			_EXCEPTION1("Input CSV-formatted nodefile \"%s\" missing column \"lat\"", strNodeFile.c_str());
		}
	}

	// Current track_id
	int iTrackId = 0;
	std::string strTrackId = "";

	// Loop through all lines
	for (int iLine = 1; ; iLine++) {
		getline(ifInput, strBuffer);
		if (ifInput.eof()) {
			break;
		}

		// Parse the line into a vector
		std::vector<std::string> vecValues;
		int iLast = 0;
		for (int i = 0; i <= strBuffer.length(); i++) {
			if ((i == strBuffer.length()) || (strBuffer[i] == ',')) {
				std::string strValue = strBuffer.substr(iLast, i-iLast);
				STLStringHelper::RemoveWhitespaceInPlace(strValue);
				vecValues.push_back(strValue);
				iLast = i+1;
			}
		}
		if (vecValues.size() != m_cdh.size()) {
			_EXCEPTION2("Mismatch in number of columns of nodefile \"%s\" on line %i", strNodeFile.c_str(), iLine);
		}

		// Check track id; if changed change the path
		if (vecValues[iColTrackId] != strTrackId) {
			iTrackId = m_pathvec.size();
			m_pathvec.resize(m_pathvec.size()+1);
			strTrackId = vecValues[iColTrackId];
		}

		_ASSERT(iTrackId < m_pathvec.size());

		Path & path = m_pathvec[iTrackId];

		path.resize(path.size()+1);

		PathNode & pathnode = path[path.size()-1];

		if (!STLStringHelper::IsInteger(vecValues[iColYear])) {
			_EXCEPTION2("year must be of type integer in nodefile \"%s\" on line %i", strNodeFile.c_str(), iLine);
		}
		if (!STLStringHelper::IsInteger(vecValues[iColMonth])) {
			_EXCEPTION2("month must be of type integer in nodefile \"%s\" on line %i", strNodeFile.c_str(), iLine);
		}
		if (!STLStringHelper::IsInteger(vecValues[iColDay])) {
			_EXCEPTION2("day must be of type integer in nodefile \"%s\" on line %i", strNodeFile.c_str(), iLine);
		}
		if (!STLStringHelper::IsInteger(vecValues[iColHour])) {
			_EXCEPTION2("hour must be of type integer in nodefile \"%s\" on line %i", strNodeFile.c_str(), iLine);
		}

		int iYear = stoi(vecValues[iColYear]);
		int iMonth = stoi(vecValues[iColMonth]);
		int iDay = stoi(vecValues[iColDay]);
		int iHour = stoi(vecValues[iColHour]);

		Time time(iYear, iMonth, iDay, iHour * 3600, caltype);

		pathnode.m_time = time;

		std::cout << strBuffer << std::endl;

		for (int v = 0; v < vecValues.size(); v++) {
			pathnode.m_vecColumnData.push_back(new ColumnDataString(vecValues[v]));
		}
	}

	// Set begin and end of all paths
	for (int p = 0; p < m_pathvec.size(); p++) {
		_ASSERT(m_pathvec[p].size() > 0);
		m_pathvec[p].m_timeStart = m_pathvec[p][0].m_time;
		m_pathvec[p].m_timeEnd = m_pathvec[p][m_pathvec[p].size()-1].m_time;
	}
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::Read(
	const std::string & strNodeFile,
	PathType ePathType,
	ColumnDataHeader & cdh,
	const std::vector<size_t> & nGridDim,
	size_t sGridSize,
	Time::CalendarType caltype
) {
	// Load in file as CSV
	if (strNodeFile.length() >= 4) {
		std::string strExt = strNodeFile.substr(strNodeFile.length()-4,4);
		if (strExt == ".csv") {
			return ReadCSV(
				strNodeFile,
				ePathType,
				cdh,
				nGridDim,
				sGridSize,
				caltype);
		}
	}

	// Cannot automatically load headers from GFDL formatted file
	if ((cdh.size() == 1) && (cdh[0] == "(auto)")) {
		_EXCEPTIONT("--in_fmt \"(auto)\" cannot be used for GFDL formatted nodefiles");
	}

	// Not implemented
	if (ePathType == PathTypeDN) {
		_EXCEPTIONT("Not implemented");
	}

	// Check that this NodeFile is uninitialized
	if ((m_cdh.size() != 0) || (m_pathvec.size() != 0)) {
		_EXCEPTIONT("Attempting to Read() on an initialized NodeFile");
	}

	// Store the PathType
	m_ePathType = ePathType;

	// Store the ColumnDataHeader data
	m_cdh = cdh;

	// String buffer
	std::string strBuffer;

	// Only support grids of dimension 1 or 2
	if ((nGridDim.size() < 1) || (nGridDim.size() > 2)) {
		_EXCEPTIONT("Grid dimension out of range:  Only grids of dimension 1 or 2 supported");
	}

	// Coordinate buffer
	std::vector<int> coord;
	coord.resize(nGridDim.size());

	// Clear the PathVector
	m_pathvec.clear();

	// DetectNodes has a single Path
	if (ePathType == PathTypeDN) {
		m_pathvec.resize(1);
	}

	// Clear the TimeToPathNodeMap
	m_mapTimeToPathNode.clear();

	// Open the file as an input stream
	std::ifstream ifInput(strNodeFile);
	if (!ifInput.is_open()) {
		_EXCEPTION1("Unable to open nodefile \"%s\"", strNodeFile.c_str());
	}

	// Current pathnode index
	size_t ixpathnode = 0;

	// Loop through all lines
	int iLine = 1;
	for (;;) {

		int nCount = 0;
		Time time;

		// Read header lines
		{
			getline(ifInput, strBuffer);
			if (ifInput.eof()) {
				break;
			}

			std::istringstream iss(strBuffer);

			// DetectNodes type
			if (ePathType == PathTypeDN) {
				int iYear;
				int iMonth;
				int iDay;
				int iHour;

				iss >> iYear;
				iss >> iMonth;
				iss >> iDay;
				iss >> nCount;

				if (!iss.good()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				iss >> iHour;

				if (!iss.eof()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				time = Time(
					iYear,
					iMonth-1,
					iDay-1,
					3600 * iHour,
					0,
					caltype);

			// StitchNodes type
			} else if (ePathType == PathTypeSN) {
				std::string strStart;
				int iYear;
				int iMonth;
				int iDay;
				int iHour;

				iss >> strStart;
				iss >> nCount;
				iss >> iYear;
				iss >> iMonth;
				iss >> iDay;

				if (!iss.good()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				iss >> iHour;

				if (strStart != "start") {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}
				if (!iss.eof()) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				m_pathvec.resize(m_pathvec.size() + 1);
				m_pathvec[m_pathvec.size()-1].m_timeStart =
					Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						caltype);

				if (nCount == 0) {
					m_pathvec[m_pathvec.size()-1].m_timeEnd =
						m_pathvec[m_pathvec.size()-1].m_timeStart;
				}
				m_pathvec[m_pathvec.size()-1].resize(nCount);
			}

			iLine++;
		}

		// Read contents under each header line
		for (int i = 0; i < nCount; i++) {

			getline(ifInput, strBuffer);
			if (ifInput.eof()) {
				_EXCEPTION2("Unexpected end-of-file on line %i of \"%s\"",
					iLine, strNodeFile.c_str());
			}

			std::istringstream iss(strBuffer);

			for (int n = 0; n < nGridDim.size(); n++) {
				iss >> coord[n];
			}

			// Note that for 2D grids the coordinate indices are swapped
			if (coord.size() == 1) {
				if (coord[0] < 0) {
					_EXCEPTION4("Negative coordinate index on line %i of \"%s\""
						" (%i/%i) (%i/%i)",
						iLine, strNodeFile.c_str(),
						coord[0], nGridDim[0]);
				}
			} else if (coord.size() == 2) {
				if ((coord[0] < 0) || (coord[1] < 0)) {
					_EXCEPTION6("Negative coordinate index on line %i of \"%s\""
						" (%i/%i) (%i/%i)",
						iLine, strNodeFile.c_str(),
						coord[0], nGridDim[1],
						coord[1], nGridDim[0]);
				}
			} else {
				_EXCEPTION();
			}

			if (iss.eof()) {
				_EXCEPTION2("Format error on line %i of \"%s\"",
					iLine, strNodeFile.c_str());
			}

			std::string strBuf;
			std::vector<std::string> vecDelimitedOutput;
			for (;;) {
				iss >> strBuf;
				vecDelimitedOutput.push_back(strBuf);
				if ((iss.eof()) || (strBuf.length() == 0)) {
					break;
				}
			}

			int nOutputSize = vecDelimitedOutput.size();

			if (cdh.size() != nOutputSize-4) {
				_EXCEPTION3("Mismatch between column header size specified in format (%i)"
					" and node file columns on line %i of \"%s\"",
					static_cast<int>(cdh.size()), iLine, strNodeFile.c_str());
			}

			// DetectNodes format output
			if (ePathType == PathTypeDN) {
				_EXCEPTIONT("Not implemented");

			// StitchNodes format input
			} else if (ePathType == PathTypeSN) {
				Path & path = m_pathvec[m_pathvec.size()-1];
				PathNode & pathnode = path[i];

				if (nOutputSize < 4) {
					_EXCEPTION2("Format error on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				// Store time
				int iYear = std::stoi(vecDelimitedOutput[nOutputSize-4]);
				int iMonth = std::stoi(vecDelimitedOutput[nOutputSize-3]);
				int iDay = std::stoi(vecDelimitedOutput[nOutputSize-2]);
				int iHour = std::stoi(vecDelimitedOutput[nOutputSize-1]);

				time = Time(
					iYear,
					iMonth-1,
					iDay-1,
					3600 * iHour,
					0,
					caltype);

				pathnode.m_time = time;

				// Add end time to Path
				if (i == nCount-1) {
					path.m_timeEnd = time;
				}

				// Store file position
				pathnode.m_fileix = ixpathnode;
				ixpathnode++;

				// Store coordinate
				if (coord.size() == 1) {
					pathnode.m_gridix = coord[0];
				} else if (coord.size() == 2) {
					pathnode.m_gridix = coord[0] + nGridDim[1] * coord[1];
				} else {
					_EXCEPTIONT("Undefined behavior for SimpleGrid dimensionality > 2");
				}

				if (pathnode.m_gridix < 0) {
					_EXCEPTION2("Negative coordinate index on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				// Store all other data as strings
				for (int j = 0; j < nOutputSize-4; j++) {
					pathnode.PushColumnDataString(
						vecDelimitedOutput[j]);
				}
			}

			iLine++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::Read(
	const std::string & strNodeFile,
	NodeFile::PathType ePathType,
	ColumnDataHeader & cdh,
	const SimpleGrid & grid,
	Time::CalendarType caltype
) {
	Read(
		strNodeFile,
		ePathType,
		cdh,
		grid.m_nGridDim,
		grid.GetSize(),
		caltype);
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::Write(
	const std::string & strNodeFile,
	const SimpleGrid * pgrid,
	const std::vector<int> * pvecColumnDataOutIx,
	FileFormat eFileFormat,
	bool fIncludeHeader
) {
	FILE * fpOutput = fopen(strNodeFile.c_str(),"w");
	if (fpOutput == NULL) {
		_EXCEPTION1("Error opening file \"%s\" for writing",
			strNodeFile.c_str());
	}

	// Output StitchNodes format file
	if (m_ePathType == PathTypeSN) {

		// GFDL file format
		if (eFileFormat == FileFormatGFDL) {
			for (int p = 0; p < m_pathvec.size(); p++) {
				Path & path = m_pathvec[p];
				if (m_pathvec[p].size() == 0) {
					_EXCEPTIONT("Zero length Path found");
				}
				fprintf(fpOutput, "start\t%i\t%i\t%i\t%i\t%i\n",
					static_cast<int>(path.size()),
					path.m_timeStart.GetYear(),
					path.m_timeStart.GetMonth(),
					path.m_timeStart.GetDay(),
					path.m_timeStart.GetSecond() / 3600);

				for (int i = 0; i < m_pathvec[p].size(); i++) {
					PathNode & pathnode = path[i];

					if (pgrid != NULL) {
						if (pgrid->m_nGridDim.size() == 1) {
							fprintf(fpOutput, "\t%lu", pathnode.m_gridix);
						} else if (pgrid->m_nGridDim.size() == 2) {
							fprintf(fpOutput, "\t%lu\t%lu",
								pathnode.m_gridix % pgrid->m_nGridDim[1],
								pathnode.m_gridix / pgrid->m_nGridDim[1]);
						}
					}

					if (pvecColumnDataOutIx == NULL) {
						for (int j = 0; j < pathnode.m_vecColumnData.size(); j++) {
							const ColumnData * pcd = pathnode.m_vecColumnData[j];
							fprintf(fpOutput, "\t%s", pcd->ToString().c_str());
						}

					} else {
						for (int j = 0; j < pvecColumnDataOutIx->size(); j++) {
							const ColumnData * pcd =
								pathnode.m_vecColumnData[(*pvecColumnDataOutIx)[j]];
							fprintf(fpOutput, "\t%s", pcd->ToString().c_str());
						}
					}

					fprintf(fpOutput, "\t%i\t%i\t%i\t%i\n",
						pathnode.m_time.GetYear(),
						pathnode.m_time.GetMonth(),
						pathnode.m_time.GetDay(),
						pathnode.m_time.GetSecond() / 3600);
				}
			}

		// CSV file format
		} else if (eFileFormat == FileFormatCSV) {

			if (fIncludeHeader) {
				fprintf(fpOutput, "track_id, year, month, day, hour");
				if (pgrid != NULL) {
					if (pgrid->m_nGridDim.size() == 1) {
						fprintf(fpOutput, ", i");
					} else {
						fprintf(fpOutput, ", i, j");
					}
				}
				if (pvecColumnDataOutIx == NULL) {
					for (int i = 0; i < m_cdh.size(); i++) {
						fprintf(fpOutput, ", %s", m_cdh[i].c_str());
					}

				} else {
					for (int j = 0; j < pvecColumnDataOutIx->size(); j++) {
						fprintf(fpOutput, ", %s", m_cdh[(*pvecColumnDataOutIx)[j]].c_str());
					}
				}
				fprintf(fpOutput,"\n");
			}

			for (int p = 0; p < m_pathvec.size(); p++) {
				Path & path = m_pathvec[p];
				if (m_pathvec[p].size() == 0) {
					_EXCEPTIONT("Zero length Path found");
				}

				for (int i = 0; i < m_pathvec[p].size(); i++) {
					PathNode & pathnode = path[i];

					fprintf(fpOutput, "%i, %i, %i, %i, %i",
						p,
						pathnode.m_time.GetYear(),
						pathnode.m_time.GetMonth(),
						pathnode.m_time.GetDay(),
						pathnode.m_time.GetSecond() / 3600);

					if (pgrid != NULL) {
						if (pgrid->m_nGridDim.size() == 1) {
							fprintf(fpOutput, ", %lu", pathnode.m_gridix);
						} else if (pgrid->m_nGridDim.size() == 2) {
							fprintf(fpOutput, ", %lu, %lu",
								pathnode.m_gridix % pgrid->m_nGridDim[1],
								pathnode.m_gridix / pgrid->m_nGridDim[1]);
						}
					}

					if (pvecColumnDataOutIx == NULL) {
						for (int j = 0; j < pathnode.m_vecColumnData.size(); j++) {
							const ColumnData * pcd = pathnode.m_vecColumnData[j];
							fprintf(fpOutput, ", %s", pcd->ToString().c_str());
						}

					} else {
						for (int j = 0; j < pvecColumnDataOutIx->size(); j++) {
							const ColumnData * pcd =
								pathnode.m_vecColumnData[(*pvecColumnDataOutIx)[j]];
							fprintf(fpOutput, ", %s", pcd->ToString().c_str());
						}
					}
					fprintf(fpOutput,"\n");
				}
			}

		} else {
			_EXCEPTIONT("Sorry, not yet implemented!");
		}
	
	} else {
		_EXCEPTIONT("Sorry, not yet implemented!");
	}

}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::GenerateTimeToPathNodeMap() {

	// Clear the existing map
	m_mapTimeToPathNode.clear();

	// Loop through all paths and path nodes
	for (int p = 0; p < m_pathvec.size(); p++) {
		Path & path = m_pathvec[p];
		for (int n = 0; n < path.size(); n++) {
			PathNode & pathnode = path[n];

			TimeToPathNodeMap::iterator iter =
				m_mapTimeToPathNode.find(pathnode.m_time);
			if (iter == m_mapTimeToPathNode.end()) {
				iter = m_mapTimeToPathNode.insert(
					TimeToPathNodeMap::value_type(
						pathnode.m_time, std::vector< std::pair<int,int> >())).first;
			}
			iter->second.push_back(std::pair<int,int>(p,n));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::ApplyTimeDelta(
	const Time & timeDelta,
	bool fAddTimeDelta
) {
	if (!fAddTimeDelta) {
		for (int p = 0; p < m_pathvec.size(); p++) {
			Path & path = m_pathvec[p];
			path.m_timeStart -= timeDelta;
			path.m_timeEnd -= timeDelta;
			for (int n = 0; n < path.size(); n++) {
				path[n].m_time -= timeDelta;
			}
		}

	} else {
		for (int p = 0; p < m_pathvec.size(); p++) {
			Path & path = m_pathvec[p];
			path.m_timeStart += timeDelta;
			path.m_timeEnd += timeDelta;
			for (int n = 0; n < path.size(); n++) {
				path[n].m_time += timeDelta;
			}
		}
	}

	// Rebuild the time to pathnode map
	if (m_mapTimeToPathNode.size() != 0) {
		GenerateTimeToPathNodeMap();
	}
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::Interpolate(
	const Time & time
) {
	// Clear existing arrays
	m_vecInterpPathId.clear();
	m_vecInterpTimeId.clear();
	m_vecInterpAlpha.clear();

	// Loop through all paths
	for (size_t p = 0; p < m_pathvec.size(); p++) {
		const Path & path = m_pathvec[p];

		if (path.size() == 0) {
			continue;
		}

		if ((path.m_timeStart <= time) && (path.m_timeEnd >= time)) {
			m_vecInterpPathId.push_back(p);

			if (path.size() == 1) {
				m_vecInterpTimeId.push_back(0);
				m_vecInterpAlpha.push_back(0.0);
				continue;
			}

			if (time < path[0].m_time) {
				m_vecInterpTimeId.push_back(0);
				m_vecInterpAlpha.push_back(0.0);
				continue;
			}

			if (time >= path[path.size()-1].m_time) {
				m_vecInterpTimeId.push_back(path.size()-2);
				m_vecInterpAlpha.push_back(1.0);
				continue;
			}

			for (size_t n = 0; n < path.size()-1; n++) {
				const PathNode & pathnodePrev = path[n];
				const PathNode & pathnodeNext = path[n+1];

				if (pathnodeNext.m_time > time) {
					if (pathnodePrev.m_time > time) {
						_EXCEPTIONT("Non-monotone time ordering in path");
					}
					if (pathnodePrev.m_time == time) {
						m_vecInterpTimeId.push_back(n);
						m_vecInterpAlpha.push_back(0.0);
						break;

					} else {
						double dTimeDelta = pathnodeNext.m_time - pathnodePrev.m_time;
						double dAlpha = (time - pathnodePrev.m_time) / dTimeDelta;

						if ((dAlpha < -ReferenceTolerance) || (dAlpha > 1.0 + ReferenceTolerance)) {
							_EXCEPTION1("Interpolation coefficient out of range (%1.5e)", dAlpha);
						}

						m_vecInterpTimeId.push_back(n);
						m_vecInterpAlpha.push_back(dAlpha);
						break;
					}
				}
			}
		}
	}

	if ((m_vecInterpPathId.size() != m_vecInterpTimeId.size()) ||
	    (m_vecInterpPathId.size() != m_vecInterpAlpha.size())
	) {
		_EXCEPTIONT("LOGIC ERROR: Mismatched path size arrays");
	}
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::InterpolatedNodeCoordinatesRad(
	const std::string & strLonName,
	const std::string & strLatName,
	std::vector<double> & vecInterpLonRad,
	std::vector<double> & vecInterpLatRad
) const {
	_ASSERT(m_vecInterpPathId.size() == m_vecInterpTimeId.size());
	_ASSERT(m_vecInterpPathId.size() == m_vecInterpAlpha.size());

	vecInterpLonRad.clear();
	vecInterpLatRad.clear();

	// Get lon and lat column indices
	int iLonIx = m_cdh.GetIndexFromString(strLonName);
	if (iLonIx == (-1)) {
		_EXCEPTION1("Column header \"%s\" not found", strLonName.c_str());
	}

	int iLatIx = m_cdh.GetIndexFromString(strLatName);
	if (iLatIx == (-1)) {
		_EXCEPTION1("Column header \"%s\" not found", strLatName.c_str());
	}

	// Loop through all paths
	for (size_t p = 0; p < m_vecInterpPathId.size(); p++) {
		_ASSERT(m_vecInterpPathId[p] < m_pathvec.size());

		const Path & path = m_pathvec[m_vecInterpPathId[p]];
		_ASSERT(path.size() > 0);

		if (path.size() == 1) {
			const PathNode & pathnode = path[0];

			double dNodeLonDeg = pathnode.GetColumnDataAsDouble(iLonIx);
			double dNodeLatDeg = pathnode.GetColumnDataAsDouble(iLatIx);

			vecInterpLonRad.push_back(DegToRad(dNodeLonDeg));
			vecInterpLatRad.push_back(DegToRad(dNodeLatDeg));

		} else {
			_ASSERT(m_vecInterpTimeId[p] < path.size()-1);
			const PathNode & pathnodePrev = path[ m_vecInterpTimeId[p] ];
			const PathNode & pathnodeNext = path[ m_vecInterpTimeId[p]+1 ];

			double dPrevLonDeg = pathnodePrev.GetColumnDataAsDouble(iLonIx);
			double dPrevLatDeg = pathnodePrev.GetColumnDataAsDouble(iLatIx);
			double dNextLonDeg = pathnodeNext.GetColumnDataAsDouble(iLonIx);
			double dNextLatDeg = pathnodeNext.GetColumnDataAsDouble(iLatIx);

			// TODO: Use great circle arc instead of line in lat-lon space
			double dNodeLonDeg =
				dPrevLonDeg * (1.0 - m_vecInterpAlpha[p])
				+ dNextLonDeg * m_vecInterpAlpha[p];
			double dNodeLatDeg =
				dPrevLatDeg * (1.0 - m_vecInterpAlpha[p])
				+ dNextLatDeg * m_vecInterpAlpha[p];

			vecInterpLonRad.push_back(DegToRad(dNodeLonDeg));
			vecInterpLatRad.push_back(DegToRad(dNodeLatDeg));
		}
	}

	_ASSERT(m_vecInterpPathId.size() == vecInterpLonRad.size());
	_ASSERT(m_vecInterpPathId.size() == vecInterpLatRad.size());
}

///////////////////////////////////////////////////////////////////////////////

void NodeFile::InterpolatedColumnDouble(
	const std::string & strHeaderName,
	std::vector<double> & vecInterpDouble
) const {
	_ASSERT(m_vecInterpPathId.size() == m_vecInterpTimeId.size());
	_ASSERT(m_vecInterpPathId.size() == m_vecInterpAlpha.size());

	vecInterpDouble.clear();

	// Get column index
	int iColIx = m_cdh.GetIndexFromString(strHeaderName);
	if (iColIx == (-1)) {
		_EXCEPTION1("Column header \"%s\" not found", strHeaderName.c_str());
	}

	// Loop through all paths
	for (size_t p = 0; p < m_vecInterpPathId.size(); p++) {
		_ASSERT(m_vecInterpPathId[p] < m_pathvec.size());

		const Path & path = m_pathvec[m_vecInterpPathId[p]];
		_ASSERT(path.size() > 0);

		if (path.size() == 1) {
			const PathNode & pathnode = path[0];
			double dValue = pathnode.GetColumnDataAsDouble(iColIx);
			vecInterpDouble.push_back(dValue);

		} else {
			_ASSERT(m_vecInterpTimeId[p] < path.size()-1);
			const PathNode & pathnodePrev = path[ m_vecInterpTimeId[p] ];
			const PathNode & pathnodeNext = path[ m_vecInterpTimeId[p]+1 ];

			double dPrevValue = pathnodePrev.GetColumnDataAsDouble(iColIx);
			double dNextValue = pathnodeNext.GetColumnDataAsDouble(iColIx);

			vecInterpDouble.push_back(
				dPrevValue * (1.0 - m_vecInterpAlpha[p])
				+ dNextValue * m_vecInterpAlpha[p]);
		}
	}

	_ASSERT(m_vecInterpPathId.size() == vecInterpDouble.size());
}

///////////////////////////////////////////////////////////////////////////////

