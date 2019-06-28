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

#include <iostream>

#include "SimpleGrid.h"
#include "AutoCurator.h"

///////////////////////////////////////////////////////////////////////////////

void ParseNodeFile(
	const std::string & strNodeFile,
	InputFileType iftype,
	const ColumnDataHeader & cdh,
	const SimpleGrid & grid,
	Time::CalendarType caltype,
	PathVector & pathvec,
	TimeToPathNodeMap & mapTimeToPathNode
) {
	// String buffer
	std::string strBuffer;

	// Only support grids of dimension 1 or 2
	if ((grid.m_nGridDim.size() < 1) || (grid.m_nGridDim.size() > 2)) {
		_EXCEPTIONT("Grid dimension out of range:  Only grids of dimension 1 or 2 supported");
	}

	// Coordinate buffer
	std::vector<int> coord;
	coord.resize(grid.m_nGridDim.size());

	// Clear the PathVector
	pathvec.clear();

	// Clear the TimeToPathNodeMap
	mapTimeToPathNode.clear();

	// Open the file as an input stream
	std::ifstream ifInput(strNodeFile);
	if (!ifInput.is_open()) {
		_EXCEPTION1("Unable to open input file \"%s\"", strNodeFile.c_str());
	}

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

			// DetectCyclonesUnstructured output
			if (iftype == InputFileTypeDCU) {
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

			// StitchNodes output
			} else if (iftype == InputFileTypeSN) {
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

				pathvec.resize(pathvec.size() + 1);
				pathvec[pathvec.size()-1].m_timeStart =
					Time(
						iYear,
						iMonth-1,
						iDay-1,
						3600 * iHour,
						0,
						caltype);

				pathvec[pathvec.size()-1].m_vecPathNodes.resize(nCount);
			}

			iLine++;
		}

		// Read contents under each header line
		for (int i = 0; i < nCount; i++) {

			getline(ifInput, strBuffer);
			if (ifInput.eof()) {
				break;
			}

			std::istringstream iss(strBuffer);

			for (int n = 0; n < grid.m_nGridDim.size(); n++) {
				iss >> coord[n];
			}

			// Note that for 2D grids the coordinate indices are swapped
			if (coord.size() == 1) {
				if ((coord[0] < 0) || (coord[0] >= grid.m_nGridDim[0])) {
					_EXCEPTION2("Coordinate index out of range on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}
			} else if (coord.size() == 2) {
				if ((coord[0] < 0) || (coord[0] >= grid.m_nGridDim[1]) ||
				    (coord[1] < 0) || (coord[1] >= grid.m_nGridDim[0])
				) {
					_EXCEPTION2("Coordinate index out of range on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
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

			// StitchNodes format input
			if (iftype == InputFileTypeSN) {
				PathNode & pathnode =
					pathvec[pathvec.size()-1].m_vecPathNodes[i];

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

				// Store coordinate
				if (coord.size() == 1) {
					pathnode.m_gridix = coord[0];
				} else if (coord.size() == 2) {
					pathnode.m_gridix = coord[0] + grid.m_nGridDim[1] * coord[1];
				} else {
					_EXCEPTIONT("Undefined behavior for SimpleGrid dimensionality > 2");
				}

				if ((pathnode.m_gridix < 0) || (pathnode.m_gridix > grid.GetSize())) {
					_EXCEPTION2("Coordinate index out of range on line %i of \"%s\"",
						iLine, strNodeFile.c_str());
				}

				// Store all other data as strings
				for (int i = 0; i < nOutputSize-4; i++) {
					pathnode.PushColumnDataString(
						vecDelimitedOutput[i]);
				}

				// Because nodes in StitchNodes format output are not ordered in
				// time, efficient data I/O requires us to reorganize the input
				// lines by time.
				TimeToPathNodeMap::iterator iter =
					mapTimeToPathNode.find(time);
				if (iter == mapTimeToPathNode.end()) {
					iter = mapTimeToPathNode.insert(
						TimeToPathNodeMap::value_type(
							time, std::vector< std::pair<int,int> >())).first;
				}
				iter->second.push_back(
					std::pair<int,int>(
						static_cast<int>(pathvec.size()-1),i));
			}

			iLine++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

