///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilities.cpp
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

#include "MeshUtilities.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

void MeshUtilities::FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	FindFaceStruct & aFindFaceStruct
) {
	// Reset the FaceStruct
	aFindFaceStruct.vecFaceIndices.clear();
	aFindFaceStruct.vecFaceLocations.clear();
	aFindFaceStruct.loc = Face::NodeLocation_Undefined;

	// Loop through all faces to find overlaps
	// Note: This algorithm can likely be dramatically improved
	for (int l = 0; l < mesh.faces.size(); l++) {
		Face::NodeLocation loc;
		int ixLocation;

		ContainsNode(
			mesh.faces[l],
			mesh.nodes,
			node,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;
		}

#ifdef VERBOSE
		printf("%i\n", l);
		printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
		printf("n0: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][0]].x,
			mesh.nodes[mesh.faces[l][0]].y,
			mesh.nodes[mesh.faces[l][0]].z);
		printf("n1: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][1]].x,
			mesh.nodes[mesh.faces[l][1]].y,
			mesh.nodes[mesh.faces[l][1]].z);
		printf("n2: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][2]].x,
			mesh.nodes[mesh.faces[l][2]].y,
			mesh.nodes[mesh.faces[l][2]].z);
		printf("n3: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][3]].x,
			mesh.nodes[mesh.faces[l][3]].y,
			mesh.nodes[mesh.faces[l][3]].z);
#endif

		if (aFindFaceStruct.loc == Face::NodeLocation_Undefined) {
			aFindFaceStruct.loc = loc;
		}

		// Node is in the interior of this face
		if (loc == Face::NodeLocation_Interior) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
			break;
		}

		// Node is on the edge of this face
		if (loc == Face::NodeLocation_Edge) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}

		// Node is at the corner of this face
		if (loc == Face::NodeLocation_Corner) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}
	}

	// Edges can only have two adjacent Faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Edge) {
		if (aFindFaceStruct.vecFaceIndices.size() != 2) {
			printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
			_EXCEPTION2("Node found on edge with %i neighboring face(s) (%i)",
				aFindFaceStruct.vecFaceIndices.size(),
				(int)(aFindFaceStruct.vecFaceIndices.size()));
		}
	}

	// Corners must have at least three adjacent Faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Corner) {
		if (aFindFaceStruct.vecFaceIndices.size() < 3) {
			printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
			_EXCEPTION1("Two Faced corner detected (%i)",
				(int)(aFindFaceStruct.vecFaceIndices.size()));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

