///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilities.h
///	\author  Paul Ullrich
///	\version August 7, 2014
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _MESHUTILITIES_H_
#define _MESHUTILITIES_H_

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

class MeshUtilities {

public:
	///	<summary>
	///		Determine if this face contains the specified Node, and whether
	///		the Node is along an edge or at a corner.
	///	</summary>
	virtual void ContainsNode(
		const Face & face,
		const NodeVector & nodevec,
		const Node & node,
		Face::NodeLocation & loc,
		int & ixLocation
	) const = 0;

	///	<summary>
	///		Find all Face indices that contain this Node.
	///	</summary>
	void FindFaceFromNode(
		const Mesh & mesh,
		const Node & node,
		FindFaceStruct & aFindFaceStruct
	);

};

///////////////////////////////////////////////////////////////////////////////

#endif
