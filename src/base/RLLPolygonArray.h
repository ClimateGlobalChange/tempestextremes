///////////////////////////////////////////////////////////////////////////////
///
///	\file    RLLPolygonArray.h
///	\author  Paul Ullrich
///	\version July 2, 2019
///
///	<remarks>
///		Copyright 2000-2019 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

struct RLLPoint {
	double lon;
	double lat;
};

typedef std::vector<RLLPoint> RLLPointVector;

///////////////////////////////////////////////////////////////////////////////

class RLLPolygonArray {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RLLPolygonArray()
	{ }

public:
	///	<summary>
	///		Load in the node array from a file, with vertices
	///		specified in degrees.
	///	</summary>
	void FromFile(
		const std::string & strFilename
	);

public:
	///	<summary>
	///		Determine if the given point is within the polygon.
	///	</summary>
	///	<param name="pt_in_degrees">
	///		Regular longitude-latitude coordinates of the testing
	///		point in degrees.
	///	</param>
	const std::string & NameOfRegionContainingPoint(
		const RLLPoint & pt_in_degrees
	);

protected:
	///	<summary>
	///		Default polygon name.
	///	</summary>
	std::string m_strDefault;

	///	<summary>
	///		Array of polygon names.
	///	</summary>
	std::vector<std::string> m_vecNames;

	///	<summary>
	///		Array of polygon nodes.
	///	</summary>
	std::vector<RLLPointVector> m_vecNodes;
};

///////////////////////////////////////////////////////////////////////////////

