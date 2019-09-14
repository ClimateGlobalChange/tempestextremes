///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilitiesFuzzy.cpp
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

#include "MeshUtilitiesFuzzy.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesFuzzy::AreNodesEqual(
	const Node & node0,
	const Node & node1
) {
	static const Real Tolerance = ReferenceTolerance;

	if ((fabs(node0.x - node1.x) < Tolerance) &&
		(fabs(node0.y - node1.y) < Tolerance) &&
		(fabs(node0.z - node1.z) < Tolerance)
	) {
		return true;
	}

/*
	double dDiffX = (node0.x - node1.x) * (node0.x - node1.x);
	double dDiffY = (node0.y - node1.y) * (node0.y - node1.y);
	double dDiffZ = (node0.z - node1.z) * (node0.z - node1.z);

	if (dDiffX + dDiffY + dDiffZ < Tolerance) {
		return true;
	}
*/
	return false;
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesFuzzy::FindNodeEdgeSide(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const Node & nodeTest
) const {
	static const Real Tolerance = ReferenceTolerance;

	if (edgetype == Edge::Type_GreatCircleArc) {
		Real dDotNorm = DotProduct(CrossProduct(nodeBegin, nodeEnd), nodeTest);

		if (dDotNorm <= - Tolerance) {
			return (-1);
		}
		if (dDotNorm < Tolerance) {
			return (0);
		}
		return (+1);

	} else if (edgetype == Edge::Type_ConstantLatitude) {
		Real dAlignment = (nodeBegin.x * nodeEnd.y - nodeEnd.x * nodeBegin.y);
		Real dDotNorm =
			dAlignment / fabs(dAlignment) * (nodeTest.z - nodeBegin.z);

		if (dDotNorm <= - Tolerance) {
			return (-1);
		}
		if (dDotNorm < Tolerance) {
			return (0);
		}
		return (+1);

	} else {
		_EXCEPTION1("Invalid EdgeType (%i)", (int)(edgetype));
	}
}

///////////////////////////////////////////////////////////////////////////////

void MeshUtilitiesFuzzy::ContainsNode(
	const Face & face,
	const NodeVector & nodevec,
	const Node & node,
	Face::NodeLocation & loc,
	int & ixLocation
) const {
	static const Real Tolerance = ReferenceTolerance;

	// Set of edges which "contain" this node
	std::set<int> setContainedEdgeIx;

	// Loop through all edges of this face
	for (int i = 0; i < face.edges.size(); i++) {

		// Paired edges means edge is invalid
		if (face.edges[i][0] == face.edges[i][1]) {
			_EXCEPTIONT("Zero Edge detected");
		}

		// Check which side of the face this edge is on
		const Node & na = nodevec[face.edges[i][0]];
		const Node & nb = nodevec[face.edges[i][1]];

		int iNodeEdgeSide =
			FindNodeEdgeSide(na, nb, face.edges[i].type, node);

		if (iNodeEdgeSide == (-1)) {
			loc = Face::NodeLocation_Exterior;
			ixLocation = 0;
			return;
		}
		if (iNodeEdgeSide == 0) {
			setContainedEdgeIx.insert(i);
		}
	}

	// Check if the node is contained on an edge
	if (setContainedEdgeIx.size() == 1) {
		loc = Face::NodeLocation_Edge;
		ixLocation = *(setContainedEdgeIx.begin());
		return;
	}

	// Node is coincident with a corner of this face
	if (setContainedEdgeIx.size() == 2) {

		std::set<int>::iterator iter;

		iter = setContainedEdgeIx.begin();
		int ix0 = *(iter);

		iter++;
		int ix1 = *(iter);

		if ((ix0 == 0) && (ix1 != 1)) {
			ixLocation = 0;
		} else {
			ixLocation = ix1;
		}

		loc = Face::NodeLocation_Corner;
		return;
	}

	// Node occurs in more than two edges; error.
	if (setContainedEdgeIx.size() > 2) {
		_EXCEPTIONT("Logic error: Node occurs in more than two edges");
	}

	// Default; node occurs in the interior of the face
	loc = Face::NodeLocation_Interior;
	ixLocation = 0;
	return;
}

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesFuzzy::CalculateEdgeIntersectionsSemiClip(
	const Node & nodeFirstBegin,
	const Node & nodeFirstEnd,
	Edge::Type typeFirst,
	const Node & nodeSecondBegin,
	const Node & nodeSecondEnd,
	Edge::Type typeSecond,
	std::vector<Node> & nodeIntersections
) {
	static const Real Tolerance = ReferenceTolerance;

	// Make a locally modifyable version of the Nodes
	Node node11;
	Node node12;
	Node node21;
	Node node22;

	// Second edge is a line of constant latitude; first is a great circle arc
	if ((typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		node11 = nodeSecondBegin;
		node12 = nodeSecondEnd;
		node21 = nodeFirstBegin;
		node22 = nodeFirstEnd;

		typeFirst = Edge::Type_GreatCircleArc;
		typeSecond = Edge::Type_ConstantLatitude;

	} else {
		node11 = nodeFirstBegin;
		node12 = nodeFirstEnd;
		node21 = nodeSecondBegin;
		node22 = nodeSecondEnd;
	}

	// Check for coincident nodes
	if (AreNodesEqual(node11, node12)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}
	if (AreNodesEqual(node21, node22)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}

	// Clear the intersection vector
	nodeIntersections.clear();

	// Both edges are great circle arcs
	if ((typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		// Cross products
		Node nodeN11xN12(CrossProduct(node11, node12));
		Node nodeN21xN22(CrossProduct(node21, node22));

		// Check for coincident lines
		Real dDot11 = DotProduct(nodeN21xN22, node11);
		Real dDot12 = DotProduct(nodeN21xN22, node12);
		Real dDot21 = DotProduct(nodeN11xN12, node21);
		Real dDot22 = DotProduct(nodeN11xN12, node22);
/*
		printf("%1.15e %1.15e %1.15e %1.15e\n",
			dDot11, dDot12, dDot21, dDot22);

		printf("\n");
		node11.Print("n11");
		node12.Print("n12");
		node21.Print("n21");
		node22.Print("n22");
*/
		// A line which is coincident with both planes
		Node nodeLine;

		// Coincident planes
		if ((fabs(dDot11) < Tolerance) &&
			(fabs(dDot21) < Tolerance)
		) {
			return true;

		// node11 is coplanar with the second arc
		} else if (fabs(dDot11) < Tolerance) {
			nodeLine = node11;

		// node12 is coplanar with the second arc
		} else if (fabs(dDot12) < Tolerance) {
			nodeLine = node12;

		// node21 is coplanar with the first arc
		} else if (fabs(dDot21) < Tolerance) {
			nodeLine = node21;

		// node22 is coplanar with the first arc
		} else if (fabs(dDot22) < Tolerance) {
			nodeLine = node22;

		// Line of intersection is the cross product of cross products
		} else {
			nodeLine = CrossProduct(nodeN11xN12, nodeN21xN22);

			// Verify coplanarity
			Real dDotDebug1 = DotProduct(nodeLine, nodeN11xN12);
			Real dDotDebug2 = DotProduct(nodeLine, nodeN21xN22);

			if ((fabs(dDotDebug1) > Tolerance) ||
				(fabs(dDotDebug2) > Tolerance)
			) {
				printf("%1.5e %1.5e\n", dDotDebug1, dDotDebug2);
				_EXCEPTIONT("Logic error");
			}
		}

		// Find the positive intersection node
		Real dMagLine = nodeLine.Magnitude();

		nodeLine.x /= dMagLine;
		nodeLine.y /= dMagLine;
		nodeLine.z /= dMagLine;

		// Check whether each podal point falls within the range
		// of the first edge
		Real dAngle11;
		Real dAngle12;

		Real dAngle1 = 1.0 - DotProduct(node11, node12);
		Real dAngle2 = 1.0 - DotProduct(node21, node22);

		// Check positive node
		dAngle11 = 1.0 - DotProduct(nodeLine, node11);
		dAngle12 = 1.0 - DotProduct(nodeLine, node12);

		if ((dAngle11 < dAngle1 + Tolerance) &&
			(dAngle12 < dAngle1 + Tolerance)
		) {
			nodeIntersections.push_back(nodeLine);
			return false;
		}

		// Check negative node
		nodeLine.x *= (-1.0);
		nodeLine.y *= (-1.0);
		nodeLine.z *= (-1.0);

		dAngle11 = 1.0 - DotProduct(nodeLine, node11);
		dAngle12 = 1.0 - DotProduct(nodeLine, node12);

		if ((dAngle11 < dAngle1 + Tolerance) &&
			(dAngle12 < dAngle1 + Tolerance)
		) {
			nodeIntersections.push_back(nodeLine);
			return false;
		}

		// No intersections
		return false;

	} else {
		_EXCEPTIONT("Not implemented");
	}
}

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesFuzzy::CalculateEdgeIntersections(
	const Node & nodeFirstBegin,
	const Node & nodeFirstEnd,
	Edge::Type typeFirst,
	const Node & nodeSecondBegin,
	const Node & nodeSecondEnd,
	Edge::Type typeSecond,
	std::vector<Node> & nodeIntersections,
	bool fIncludeFirstBeginNode
) {
	static const Real Tolerance = ReferenceTolerance;

	// Make a locally modifyable version of the Nodes
	Node node11;
	Node node12;
	Node node21;
	Node node22;

	// Second edge is a line of constant latitude; first is a great circle arc
	if ((typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		node11 = nodeSecondBegin;
		node12 = nodeSecondEnd;
		node21 = nodeFirstBegin;
		node22 = nodeFirstEnd;

		typeFirst = Edge::Type_GreatCircleArc;
		typeSecond = Edge::Type_ConstantLatitude;

	} else {
		node11 = nodeFirstBegin;
		node12 = nodeFirstEnd;
		node21 = nodeSecondBegin;
		node22 = nodeSecondEnd;
	}

	// Check for coincident nodes
	if (AreNodesEqual(node11, node12)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}
	if (AreNodesEqual(node21, node22)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}

	// Clear the intersection vector
	nodeIntersections.clear();

	// Both edges are great circle arcs
	if ((typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		// Cross products
		Node nodeN11xN12(CrossProduct(node11, node12));
		Node nodeN21xN22(CrossProduct(node21, node22));

		// Check for coincident lines
		Real dDot11 = DotProduct(nodeN21xN22, node11);
		Real dDot12 = DotProduct(nodeN21xN22, node12);
		Real dDot21 = DotProduct(nodeN11xN12, node21);
		Real dDot22 = DotProduct(nodeN11xN12, node22);
/*
		printf("%1.15e %1.15e %1.15e %1.15e\n",
			dDot11, dDot12, dDot21, dDot22);

		printf("\n");
		node11.Print("n11");
		node12.Print("n12");
		node21.Print("n21");
		node22.Print("n22");
*/
		// A line which is coincident with both planes
		Node nodeLine;

		// Coincident planes
		if ((fabs(dDot11) < Tolerance) &&
			(fabs(dDot21) < Tolerance)
		) {
			return true;

		// node11 is coplanar with the second arc
		} else if (fabs(dDot11) < Tolerance) {
			nodeLine = node11;

		// node12 is coplanar with the second arc
		} else if (fabs(dDot12) < Tolerance) {
			nodeLine = node12;

		// node21 is coplanar with the first arc
		} else if (fabs(dDot21) < Tolerance) {
			nodeLine = node21;

		// node22 is coplanar with the first arc
		} else if (fabs(dDot22) < Tolerance) {
			nodeLine = node22;

		// Line of intersection is the cross product of cross products
		} else {
			nodeLine = CrossProduct(nodeN11xN12, nodeN21xN22);

			// Verify coplanarity
			Real dDotDebug1 = DotProduct(nodeLine, nodeN11xN12);
			Real dDotDebug2 = DotProduct(nodeLine, nodeN21xN22);

			if ((fabs(dDotDebug1) > Tolerance) ||
				(fabs(dDotDebug2) > Tolerance)
			) {
				printf("%1.5e %1.5e\n", dDotDebug1, dDotDebug2);
				_EXCEPTIONT("Logic error");
			}
		}
/*
		// True if nodeLine is in the fan of node11 and node12.  False if
		// -nodeLine is in the fan of node11 and node12.
		bool fFirstPosModeInRange = true;
		bool fSecondPosModeInRange = true;

		Real dDenom;
		Real dC;
		Real dD;

		// Determine if nodeLine is in the fan of [node11, node12]
		if (fabs(nodeN11xN12.x) > Tolerance) {
			dDenom = nodeN11xN12.x;

			dC = nodeLine.y * node12.z - nodeLine.z * node12.y;
			dD = nodeLine.z * node11.y - nodeLine.y * node11.z;

		} else if (fabs(nodeN11xN12.y) > Tolerance) {
			dDenom = nodeN11xN12.y;

			dC = nodeLine.z * node12.x - nodeLine.x * node12.z;
			dD = nodeLine.x * node11.z - nodeLine.z * node11.x;

		} else if (fabs(nodeN11xN12.z) > Tolerance) {
			dDenom = nodeN11xN12.z;

			dC = nodeLine.x * node12.y - nodeLine.y * node12.x;
			dD = nodeLine.y * node11.x - nodeLine.x * node11.y;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}

		if ((dC <= -Tolerance) && (dD >= Tolerance)) {
			return false;
		} else if ((dC >= Tolerance) && (dD <= -Tolerance)) {
			return false;
		} else if ((dC >= -Tolerance) && (dD >= -Tolerance)) {
			fFirstPosModeInRange = (dDenom > 0.0);
		} else {
			fFirstPosModeInRange = (dDenom < 0.0);
		}

		// If we have reached this point then nodeLine is in the fan
		// of [node11, node12].  Now determine if nodeLine is in the
		// fan of [node21, node22].
		if (fabs(nodeN21xN22.x) > Tolerance) {
			dDenom = nodeN21xN22.x;

			dC = nodeLine.y * node22.z - nodeLine.z * node22.y;
			dD = nodeLine.z * node21.y - nodeLine.y * node21.z;

		} else if (fabs(nodeN21xN22.y) > Tolerance) {
			dDenom = nodeN21xN22.y;

			dC = nodeLine.z * node22.x - nodeLine.x * node22.z;
			dD = nodeLine.x * node21.z - nodeLine.z * node21.x;

		} else if (fabs(nodeN21xN22.z) > Tolerance) {
			dDenom = nodeN21xN22.z;

			dC = nodeLine.x * node22.y - nodeLine.y * node22.x;
			dD = nodeLine.y * node21.x - nodeLine.x * node21.y;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}

		if ((dC <= -Tolerance) && (dD >= Tolerance)) {
			return false;
		} else if ((dC >= Tolerance) && (dD <= -Tolerance)) {
			return false;
		} else if ((dC >= -Tolerance) && (dD >= -Tolerance)) {
			fSecondPosModeInRange = (dDenom > 0.0);
		} else {
			fSecondPosModeInRange = (dDenom < 0.0);
		}

		// Verify fans are not antipodal
		if (fFirstPosModeInRange != fSecondPosModeInRange) {
			return false;
		}

		// Solution exists
		if (!fFirstPosModeInRange) {
			nodeLine.x = - nodeLine.x;
			nodeLine.y = - nodeLine.y;
			nodeLine.z = - nodeLine.z;
		} 

		nodeIntersections.push_back(nodeLine);
		return false;
*/

		// Find the intersection point
		Real dMagLine = nodeLine.Magnitude();

		nodeLine.x /= dMagLine;
		nodeLine.y /= dMagLine;
		nodeLine.z /= dMagLine;

		// Check whether each podal point falls within the range
		Real dAngle11;
		Real dAngle12;
		Real dAngle21;
		Real dAngle22;

		Real dAngle1 = 1.0 - DotProduct(node11, node12);
		Real dAngle2 = 1.0 - DotProduct(node21, node22);

		// Check positive node
		dAngle11 = 1.0 - DotProduct(nodeLine, node11);
		dAngle12 = 1.0 - DotProduct(nodeLine, node12);
		dAngle21 = 1.0 - DotProduct(nodeLine, node21);
		dAngle22 = 1.0 - DotProduct(nodeLine, node22);

		if ((dAngle11 < dAngle1 + Tolerance) &&
			(dAngle12 < dAngle1 + Tolerance) &&
			(dAngle21 < dAngle2 + Tolerance) &&
			(dAngle22 < dAngle2 + Tolerance)
		) {
			if ((AreNodesEqual(nodeLine, nodeFirstBegin)) &&
				!fIncludeFirstBeginNode
			) {
				return false;
			}

			nodeIntersections.push_back(nodeLine);
			return false;
		}

		// Check negative node
		nodeLine.x *= (-1.0);
		nodeLine.y *= (-1.0);
		nodeLine.z *= (-1.0);

		dAngle11 = 1.0 - DotProduct(nodeLine, node11);
		dAngle12 = 1.0 - DotProduct(nodeLine, node12);
		dAngle21 = 1.0 - DotProduct(nodeLine, node21);
		dAngle22 = 1.0 - DotProduct(nodeLine, node22);

		if ((dAngle11 < dAngle1 + Tolerance) &&
			(dAngle12 < dAngle1 + Tolerance) &&
			(dAngle21 < dAngle2 + Tolerance) &&
			(dAngle22 < dAngle2 + Tolerance)
		) {
			if ((AreNodesEqual(nodeLine, nodeFirstBegin)) &&
				!fIncludeFirstBeginNode
			) {
				return false;
			}

			nodeIntersections.push_back(nodeLine);
			return false;
		}

		// No intersections
		return false;

/*
		// n11 dot n12
		Real dN11oN12 =
			+ node11.x * node12.x
			+ node11.y * node12.y
			+ node11.z * node12.z;

		// Cross product of second vectors
		Node nodeN21xN22(
			+ node21.y * node22.z - node21.z * node22.y,
			- node21.x * node22.z + node21.z * node22.x,
			+ node21.x * node22.y - node21.y * node22.x);

		// Other Cross products
		Node nodeN12xN21(
			+ node12.y * node21.z - node12.z * node21.y,
			- node12.x * node21.z + node12.z * node21.x,
			+ node12.x * node21.y - node12.y * node21.x);

		Node nodeN12xN22(
			+ node12.y * node22.z - node12.z * node22.y,
			- node12.x * node22.z + node12.z * node22.x,
			+ node12.x * node22.y - node12.y * node22.x);

		// n12 dot (n21 cross n22)
		Real dN12oN21xN22 =
			+ node12.x * nodeN21xN22.x
			+ node12.y * nodeN21xN22.y
			+ node12.z * nodeN21xN22.z;

		// n11 dot (n21 cross n22)
		Real dN11oN21xN22 =
			+ node11.x * nodeN21xN22.x
			+ node11.y * nodeN21xN22.y
			+ node11.z * nodeN21xN22.z;

		// Check if all four vectors lay in the same plane (coincident lines)
		if ((fabs(dN11oN21xN22) < Tolerance) &&
		    (fabs(dN12oN21xN22) < Tolerance)
		) {
			// If coincident, check for intersections
			Node nodeN11xN12(
				+ node11.y * node12.z - node11.z * node12.y,
				- node11.x * node12.z + node11.z * node12.x,
				+ node11.x * node12.y - node11.y * node12.x);

			// Use the largest element of the cross product
			Real dA0;
			Real dA1;
			Real dB0;
			Real dB1;

			if ((fabs(nodeN11xN12.x) > fabs(nodeN11xN12.y)) &&
				(fabs(nodeN11xN12.x) > fabs(nodeN11xN12.z))
			) {
				dA0 = - node12.y * node21.z + node12.z * node21.y;
				dA1 = - node12.y * node22.z + node12.z * node22.y;

				dB0 = + node11.y * node21.z - node11.z * node21.y;
				dB1 = + node11.y * node22.z - node11.z * node22.y;

				dA0 /= - nodeN11xN12.x;
				dA1 /= - nodeN11xN12.x;
				dB0 /= - nodeN11xN12.x;
				dB1 /= - nodeN11xN12.x;

			} else if (
				(fabs(nodeN11xN12.y) > fabs(nodeN11xN12.x)) &&
				(fabs(nodeN11xN12.y) > fabs(nodeN11xN12.z))
			) {
				dA0 = - node12.x * node21.z + node12.z * node21.x;
				dA1 = - node12.x * node22.z + node12.z * node22.x;

				dB0 = + node11.x * node21.z - node11.z * node21.x;
				dB1 = + node11.x * node22.z - node11.z * node22.x;

				dA0 /= - nodeN11xN12.y;
				dA1 /= - nodeN11xN12.y;
				dB0 /= - nodeN11xN12.y;
				dB1 /= - nodeN11xN12.y;

			} else {
				dA0 = - node12.x * node21.y + node12.y * node21.x;
				dA1 = - node12.x * node22.y + node12.y * node22.x;

				dB0 = + node11.x * node21.y - node11.y * node21.x;
				dB1 = + node11.x * node22.y - node11.y * node22.x;

				dA0 /= nodeN11xN12.z;
				dA1 /= nodeN11xN12.z;
				dB0 /= nodeN11xN12.z;
				dB1 /= nodeN11xN12.z;
			}

			// Insert in order from FirstBegin to FirstEnd
			if (dA0 < dA1) {
				if ((dA0 > -Tolerance) && (dB0 > -Tolerance)) {
					nodeIntersections.push_back(node21);
				}
				if ((dA1 > -Tolerance) && (dB1 > -Tolerance)) {
					nodeIntersections.push_back(node22);
				}
			} else {
				if ((dA1 > -Tolerance) && (dB1 > -Tolerance)) {
					nodeIntersections.push_back(node22);
				}
				if ((dA0 > -Tolerance) && (dB0 > -Tolerance)) {
					nodeIntersections.push_back(node21);
				}
			}

			// Remove intersections that coincide with the first begin node
			if (!fIncludeFirstBeginNode) {
				for (int i = 0; i < nodeIntersections.size(); i++) {
					if (nodeIntersections[i] == nodeFirstBegin) {
						nodeIntersections.erase(nodeIntersections.begin()+i);
					}
				}
			}

			return true;
		}

		// Solution coefficients
		Real dA0;
		Real dB0;
		Real dC0;
		Real dD0;

		Real dNumerC =
			+ node11.x * nodeN12xN22.x
			+ node11.y * nodeN12xN22.y
			+ node11.z * nodeN12xN22.z;

		Real dNumerD =
			+ node11.x * nodeN12xN21.x
			+ node11.y * nodeN12xN21.y
			+ node11.z * nodeN12xN21.z;

		// node12 is farthest from the plane defining node21 and node22
		if (fabs(dN11oN21xN22) < fabs(dN12oN21xN22)) {

			Real dNumerB =
				+ node11.x * nodeN21xN22.x
				+ node11.y * nodeN21xN22.y
				+ node11.z * nodeN21xN22.z;

			Real dMB = - dNumerB / dN12oN21xN22;
			Real dMC = - dNumerC / dN12oN21xN22;
			Real dMD = + dNumerD / dN12oN21xN22;

			Real dDenom = dMB * dMB + 2.0 * dMB * dN11oN12 + 1.0;

			dA0 = 1.0 / sqrt(dDenom);
			dB0 = dA0 * dMB;
			dC0 = dA0 * dMC;
			dD0 = dA0 * dMD;

		// node11 is farthest from the plane defining node21 and node22
		} else {

			Real dNumerA =
				+ node12.x * nodeN21xN22.x
				+ node12.y * nodeN21xN22.y
				+ node12.z * nodeN21xN22.z;

			Real dMA = - dNumerA / dN11oN21xN22;
			Real dMC = + dNumerC / dN11oN21xN22;
			Real dMD = - dNumerD / dN11oN21xN22;

			Real dDenom = dMA * dMA + 2.0 * dMA * dN11oN12 + 1.0;

			dB0 = 1.0 / sqrt(dDenom);
			dA0 = dB0 * dMA;
			dC0 = dB0 * dMC;
			dD0 = dB0 * dMD;
		}

		// Check if first solution lies within interval
		if ((dA0 > -Tolerance) &&
			(dB0 > -Tolerance) &&
			(dC0 > -Tolerance) &&
			(dD0 > -Tolerance)
		) {
			nodeIntersections.push_back(Node(
				(node11.x * dA0 + node12.x * dB0),
				(node11.y * dA0 + node12.y * dB0),
				(node11.z * dA0 + node12.z * dB0)));

		// Check if second solution lies within interval
		} else if (
			(dA0 < Tolerance) &&
			(dB0 < Tolerance) &&
			(dC0 < Tolerance) &&
			(dD0 < Tolerance)
		) {
			nodeIntersections.push_back(Node(
				- (node11.x * dA0 + node12.x * dB0),
				- (node11.y * dA0 + node12.y * dB0),
				- (node11.z * dA0 + node12.z * dB0)));
		}
*/
	// First edge is a line of constant latitude; second is a great circle arc
	} else if (
		(typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {
		// Check for coincident edges (edges along the equator)
		if ((fabs(node11.z) < Tolerance) &&
			(fabs(node12.z) < Tolerance) &&
			(fabs(node21.z) < Tolerance)
		) {
			return true;
		}

		// Cross product of basis vectors for great circle plane
		Real dCrossX = node11.y * node12.z - node11.z * node12.y;
		Real dCrossY = node11.z * node12.x - node11.x * node12.z;
		Real dCrossZ = node11.x * node12.y - node11.y * node12.x;

		// Maximum Z value reached by great circle arc along sphere
		Real dAbsCross2 =
			dCrossX * dCrossX + dCrossY * dCrossY + dCrossZ * dCrossZ;

		Real dAbsEqCross2 =
			dCrossX * dCrossX + dCrossY * dCrossY;

		Real dApexZ = sqrt(dAbsEqCross2 / dAbsCross2);
/*
		// Check if apex is under this latitude (zero intersections)
		if (fabs(node21.z) - dApexZ > Tolerance) {
			return false;
		}

		// Find apex of great circle arc
		//Real dApexX = - dCrossX / dCrossZ;
		//Real dApexY = - dCrossY / dCrossZ;

		Real dApexExpectedLength = sqrt(1.0 - dApexZ * dApexZ);
		Real dApexX = - dCrossX * dApexExpectedLength / sqrt(dAbsEqCross2);
		Real dApexY = - dCrossY * dApexExpectedLength / sqrt(dAbsEqCross2);

		Real dApexMag = dApexX * dApexX + dApexY * dApexY + dApexZ * dApexZ;
		if (fabs(dApexMag - 1.0) > Tolerance) {
			printf("Magnitude: %1.15e\n", dApexMag);
			_EXCEPTIONT("Logic error");
		}
*/
/*
		// Check if apex is exactly at this latitude (one intersection)
		if ((node21.z > 0.0) && (fabs(node21.z - dApexZ) <= Tolerance)) {
			nodeIntersections.resize(1);
			nodeIntersections[0].x = dApexX;
			nodeIntersections[0].y = dApexY;
			nodeIntersections[0].z = dApexZ;

		} else if ((node21.z < 0.0) && (fabs(node21.z - dApexZ) <= Tolerance) {
			nodeIntersections.resize(1);
			nodeIntersections[0].x = - dApexX;
			nodeIntersections[0].y = - dApexY;
			nodeIntersections[0].z = - dApexZ;
		}
*/
		// node12.z is larger than node11.z
		if (fabs(node11.z) < fabs(node12.z)) {

			// Quadratic coefficients, used to solve for A
			Real dDTermA = (dCrossY * dCrossY + dCrossX * dCrossX)
				/ (node12.z * node12.z);

			Real dDTermB = + 2.0 * node21.z / (node12.z * node12.z) * (
				- node12.x * dCrossY + node12.y * dCrossX);

			Real dDTermC = node21.z * node21.z / (node12.z * node12.z) - 1.0;

			Real dDisc = dDTermB * dDTermB - 4.0 * dDTermA * dDTermC;

			Real dCross2 = node21.x * node22.y - node21.y * node22.x;

			// Only one solution
			//if (fabs(dDisc) < Tolerance) {
			if (fabs(dApexZ - fabs(node21.z)) < Tolerance) {

				// Components of intersection in node1 basis
				Real dA = - dDTermB / (2.0 * dDTermA);

				Real dB = (-dA * node11.z + node21.z) / node12.z;

				// Components of intersection in (1,0)-(0,1) basis
				Real dC = (-dA * dCrossY + node12.x * node21.z) / node12.z;

				Real dD = ( dA * dCrossX + node12.y * node21.z) / node12.z;

				// Components of intersection in node2 basis
				Real dE = ( dC * node22.y - dD * node22.x) / dCross2;

				Real dF = (-dC * node21.y + dD * node21.x) / dCross2;

				if ((dA > -Tolerance) &&
					(dB > -Tolerance) &&
					(dE > -Tolerance) &&
					(dF > -Tolerance)
				) {
					nodeIntersections.resize(1);
					nodeIntersections[0].x = dC;
					nodeIntersections[0].y = dD;
					nodeIntersections[0].z = node21.z;
				}

			// Possibly multiple solutiosn
			} else {
				Real dSqrtDisc =
					sqrt(dDTermB * dDTermB - 4.0 * dDTermA * dDTermC);

				// Components of intersection in node1 basis
				Real dA0 = (- dDTermB + dSqrtDisc) / (2.0 * dDTermA);
				Real dA1 = (- dDTermB - dSqrtDisc) / (2.0 * dDTermA);

				Real dB0 = (-dA0 * node11.z + node21.z) / node12.z;
				Real dB1 = (-dA1 * node11.z + node21.z) / node12.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				Real dC0 = (-dA0 * dCrossY + node12.x * node21.z) / node12.z;
				Real dC1 = (-dA1 * dCrossY + node12.x * node21.z) / node12.z;

				Real dD0 = ( dA0 * dCrossX + node12.y * node21.z) / node12.z;
				Real dD1 = ( dA1 * dCrossX + node12.y * node21.z) / node12.z;

				// Components of intersection in node2 basis
				Real dE0 = ( dC0 * node22.y - dD0 * node22.x) / dCross2;
				Real dE1 = ( dC1 * node22.y - dD1 * node22.x) / dCross2;

				Real dF0 = (-dC0 * node21.y + dD0 * node21.x) / dCross2;
				Real dF1 = (-dC1 * node21.y + dD1 * node21.x) / dCross2;

				if ((dA0 > -Tolerance) &&
					(dB0 > -Tolerance) &&
					(dE0 > -Tolerance) &&
					(dF0 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC0, dD0, node21.z));
				}
				if ((dA1 > -Tolerance) &&
					(dB1 > -Tolerance) &&
					(dE1 > -Tolerance) &&
					(dF1 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC1, dD1, node21.z));
				}
			}

		// node11.z is larger than node12.z
		} else {

			// Quadratic coefficients, used to solve for B
			Real dDTermA = (dCrossY * dCrossY + dCrossX * dCrossX)
				/ (node11.z * node11.z);

			Real dDTermB = - 2.0 * node21.z / (node11.z * node11.z) * (
				- node11.x * dCrossY + node11.y * dCrossX);

			Real dDTermC = node21.z * node21.z / (node11.z * node11.z) - 1.0;

			Real dDisc = dDTermB * dDTermB - 4.0 * dDTermA * dDTermC;

			Real dCross2 = node21.x * node22.y - node21.y * node22.x;

			// Only one solution
			//if (fabs(dDisc) < Tolerance) {
			if (fabs(dApexZ - fabs(node21.z)) < Tolerance) {

				// Components of intersection in node1 basis
				Real dB = - dDTermB / (2.0 * dDTermA);

				Real dA = (-dB * node12.z + node21.z) / node11.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				Real dC = (dB * dCrossY + node11.x * node21.z) / node11.z;

				Real dD = (-dB * dCrossX + node11.y * node21.z) / node11.z;

				// Components of intersection in node2 basis
				Real dE = ( dC * node22.y - dD * node22.x) / dCross2;

				Real dF = (-dC * node21.y + dD * node21.x) / dCross2;

				if ((dA > -Tolerance) &&
					(dB > -Tolerance) &&
					(dE > -Tolerance) &&
					(dF > -Tolerance)
				) {
					nodeIntersections.resize(1);
					nodeIntersections[0].x = dC;
					nodeIntersections[0].y = dD;
					nodeIntersections[0].z = node21.z;
				}

			// Two solutions
			} else {
				Real dSqrtDisc = sqrt(dDisc);

				// Components of intersection in node1 basis
				Real dB0 = (- dDTermB + dSqrtDisc) / (2.0 * dDTermA);
				Real dB1 = (- dDTermB - dSqrtDisc) / (2.0 * dDTermA);

				Real dA0 = (-dB0 * node12.z + node21.z) / node11.z;
				Real dA1 = (-dB1 * node12.z + node21.z) / node11.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				Real dC0 = (dB0 * dCrossY + node11.x * node21.z) / node11.z;
				Real dC1 = (dB1 * dCrossY + node11.x * node21.z) / node11.z;

				Real dD0 = (-dB0 * dCrossX + node11.y * node21.z) / node11.z;
				Real dD1 = (-dB1 * dCrossX + node11.y * node21.z) / node11.z;

				// Components of intersection in node2 basis
				Real dE0 = ( dC0 * node22.y - dD0 * node22.x) / dCross2;
				Real dE1 = ( dC1 * node22.y - dD1 * node22.x) / dCross2;

				Real dF0 = (-dC0 * node21.y + dD0 * node21.x) / dCross2;
				Real dF1 = (-dC1 * node21.y + dD1 * node21.x) / dCross2;

				if ((dA0 > -Tolerance) &&
					(dB0 > -Tolerance) &&
					(dE0 > -Tolerance) &&
					(dF0 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC0, dD0, node21.z));
				}
				if ((dA1 > -Tolerance) &&
					(dB1 > -Tolerance) &&
					(dE1 > -Tolerance) &&
					(dF1 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC1, dD1, node21.z));
				}
			}
		}

	// Both edges are lines of constant latitude
	} else if (
		(typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {
		if (fabs(node11.z - node21.z) < Tolerance) {
			return true;
		} else {
			return false;
		}

	// Unknown
	} else {
		_EXCEPTION2("Invalid Edge::Type (%i, %i)",
			typeFirst, typeSecond);
	}

	// If the begin node is not to be included, erase from intersections
	if (!fIncludeFirstBeginNode) {
		for (int i = 0; i < nodeIntersections.size(); i++) {
			if (AreNodesEqual(nodeIntersections[i], nodeFirstBegin)) {
				nodeIntersections.erase(nodeIntersections.begin()+i);
			}
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesFuzzy::FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {
	static const Real Tolerance = ReferenceTolerance;

	// Get the reference point
	const Node & nodeBegin = mesh.nodes[ixNode];

	// Get the set of faces adjacent this node
	const std::set<int> & setNearbyFaces = mesh.revnodearray[ixNode];

	if (setNearbyFaces.size() < 3) {
		_EXCEPTIONT("Insufficient Faces at Corner; at least three Faces expected");
	}

	// Direction along second curve projected onto
	// surface of the sphere at nodeBegin
	Node nodeDirA;

	//GetLocalDirection(nodeBegin, nodeEnd, nodeRef, edgetype, nodeDirA);
	GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDirA);
/*
	Real dDotNbNb = DotProduct(nodeBegin, nodeBegin);
	Real dDotNeNb = DotProduct(nodeEnd, nodeBegin);

	Node nodeLocalE =
		ScalarProduct(dDotNbNb, nodeEnd)
		- ScalarProduct(dDotNeNb, nodeBegin);
*/
	// Loop through all faces
	std::set<int>::const_iterator iter = setNearbyFaces.begin();
	for (; iter != setNearbyFaces.end(); iter++) {

		const Face & face = mesh.faces[*iter];

		int nEdges = face.edges.size();

		// Find the edge attached to this node
		int ixThisEdge = (-1);
		for (int j = 0; j < nEdges; j++) {
			int jNext = (j + 1) % nEdges;
			if (face[j] == face[jNext]) {
				continue;
			}
			if (ixNode == face[j]) {
				ixThisEdge = j;
				break;
			}
		}
		if (ixThisEdge == (-1)) {
			_EXCEPTIONT("Logic error");
		}

		int ixPrevEdge = (ixThisEdge + nEdges - 1) % nEdges;
		for (int j = ixPrevEdge; j != ixThisEdge; j--) {
			if (face[j] != face[ixThisEdge]) {
				break;
			}
			if (j == (-1)) {
				j = nEdges - 1;
			}
		}

		if (face[ixPrevEdge] == face[ixThisEdge]) {
			_EXCEPTIONT("Invalid Face: Only one distinct Edge found");
		}

		// Bounding edges on Node
		const Edge & edgeThis = face.edges[ixThisEdge];
		const Edge & edgePrev = face.edges[ixPrevEdge];

		if (edgeThis[0] != edgePrev[1]) {
			_EXCEPTIONT("Invalid Node indices on Edges");
		}

		const Node & node0 = mesh.nodes[edgePrev[0]];
		const Node & node1 = mesh.nodes[edgePrev[1]];
		const Node & node2 = mesh.nodes[edgeThis[1]];

		printf("%i %i %i\n", edgePrev[0], edgePrev[1], edgeThis[1]);

		if (!AreNodesEqual(node1, nodeBegin)) {
			_EXCEPTIONT("Logic error");
		}
/*
		// Local left and right vectors (up to a scaling factor)
		Real dDotN0Nb = DotProduct(node0, nodeBegin);
		Real dDotN2Nb = DotProduct(node2, nodeBegin);

		Node nodeLocal0 =
			ScalarProduct(dDotNbNb, node0)
			- ScalarProduct(dDotN0Nb, nodeBegin);

		Node nodeLocal2 =
			ScalarProduct(dDotNbNb, node2)
			- ScalarProduct(dDotN2Nb, nodeBegin);

		Node nodeLocalCross = CrossProduct(nodeLocal0, nodeLocal2);

		// Determine if nodeLocalE is in the fan of [nodeLocal0, nodeLocal2]
		Real dDenom;
		Real d2;
		Real d0;
		if (fabs(nodeLocalCross.x) > Tolerance) {
			dDenom = nodeLocalCross.x;

			d2 = nodeLocalE.y * nodeLocal2.z - nodeLocalE.z * nodeLocal2.y;
			d0 = nodeLocalE.z * nodeLocal0.y - nodeLocalE.y * nodeLocal0.z;

		} else if (fabs(nodeLocalCross.y) > Tolerance) {
			dDenom = nodeLocalCross.y;

			d2 = nodeLocalE.z * nodeLocal2.x - nodeLocalE.x * nodeLocal2.z;
			d0 = nodeLocalE.x * nodeLocal0.z - nodeLocalE.z * nodeLocal0.x;

		} else if (fabs(nodeLocalCross.z) > Tolerance) {
			dDenom = nodeLocalCross.z;

			d2 = nodeLocalE.x * nodeLocal2.y - nodeLocalE.y * nodeLocal2.x;
			d0 = nodeLocalE.y * nodeLocal0.x - nodeLocalE.x * nodeLocal0.y;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}

		if (dDenom > 0.0) {
			if ((d2 > -Tolerance) && (d0 > Tolerance)) {
				return (*iter);
			}

		} else {
			if ((d2 < Tolerance) && (d0 < -Tolerance)) {
				return (*iter);
			}
		}
*/

		// Direction vectors towards each of the nodes
		Node nodeDirL;
		Node nodeDirR;

		GetLocalDirection(node1, node0, edgePrev.type, nodeDirL);
		GetLocalDirection(node1, node2, edgeThis.type, nodeDirR);

		// Verify coplanarity
		Node nodeX(CrossProduct(nodeDirA, nodeDirL));

		Real dMagA = nodeDirA.Magnitude();
		Real dMagL = nodeDirL.Magnitude();
		Real dMagR = nodeDirR.Magnitude();

		if (fabs(dMagA) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of path Edge (possible zero Edge)");
		}
		if (fabs(dMagL) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of incoming Edge");
		}
		if (fabs(dMagR) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of outgoing Edge");
		}

		Real dDotLR = DotProduct(nodeDirL, nodeDirR);
		Real dDotLA = DotProduct(nodeDirL, nodeDirA);
		Real dDotRA = DotProduct(nodeDirR, nodeDirA);

		Real dNormDotLR = dDotLR / dMagL / dMagR;
		Real dNormDotLA = dDotLA / dMagL / dMagA;
		Real dNormDotRA = dDotRA / dMagR / dMagA;

		// These values mimic the monotone structure of acos()
		Real dAngleLR = 1.0 - dNormDotLR;
		Real dAngleLA = 1.0 - dNormDotLA;
		Real dAngleRA = 1.0 - dNormDotRA;

		// Check that this branch does not lie exterior to the corner
		if (dAngleLA + dAngleRA > 1.1 * dAngleLR) {
			continue;
		}

#ifdef VERBOSE
		printf("Face: %i\n", (*iter));
		printf("Type: %i %i %i\n", edgePrev.type, edgeThis.type, edgetype);
		nodeDirA.Print("DirA");
		nodeDirL.Print("DirL");
		nodeDirR.Print("DirR");
		printf("Angles: %1.15e %1.15e %1.15e\n", dAngleLR, dAngleLA, dAngleRA);
		//printf("Mags: %1.15e %1.15e %1.15e\n", dMagL, dMagR, dMagA);

		Real dDotPlanar = DotProduct(nodeX, nodeDirR);

		if (fabs(dDotPlanar) > Tolerance) {
			printf("Planarity: %1.10Le\n", dDotPlanar);
			_EXCEPTIONT("Planarity failure");
		}
#endif

		// Sanity check for angles
		if (dAngleLR < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (dAngleLA < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (dAngleRA < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (fabs(dAngleLR) < Tolerance) {
			_EXCEPTIONT("Zero angle detected, cannot continue");
		}
		if (fabs(dAngleLR - 2.0) < Tolerance) {
			_EXCEPTIONT("180 degree angle detected, cannot continue");
		}

		// This line makes zero angle with the L edge (special case)
		if (dAngleLA <= Tolerance) {

			// Coincident edges
			if (edgetype == edgePrev.type) {
				continue;
			}

			// Great circle arc emerging along a line of constant latitude
			if ((edgetype == Edge::Type_GreatCircleArc) &&
				(edgePrev.type == Edge::Type_ConstantLatitude)
			) {
				if ((node2.z > node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node2.z < node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else {
					continue;
				}
			}

			// Line of constant latitude emerging along a great circle arc
			if ((edgetype == Edge::Type_ConstantLatitude) &&
				(edgePrev.type == Edge::Type_GreatCircleArc)
			) {
				if ((node2.z < node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node2.z > node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else {
					continue;
				}
			}
		}

		// This line makes zero angle with the R edge (special case)
		if (dAngleRA <= Tolerance) {

			// Coincident edges
			if (edgetype == edgeThis.type) {
				return (*iter);
			}

			// Great circle arc emerging along a line of constant latitude
			if ((edgetype == Edge::Type_GreatCircleArc) &&
				(edgeThis.type == Edge::Type_ConstantLatitude)
			) {
				if ((node0.z > node1.z) && (nodeBegin.z < Tolerance)) {
					return (*iter);
				} else if ((node0.z < node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else if (fabs(nodeBegin.z) <= Tolerance) {
					return (*iter);
				} else {
					continue;
				}
			}

			// Line of constant latitude emerging along a great circle arc
			if ((edgetype == Edge::Type_ConstantLatitude) &&
				(edgeThis.type == Edge::Type_GreatCircleArc)
			) {
				if ((node0.z < node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node0.z > node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else if (fabs(nodeBegin.z) <= Tolerance) {
					return (*iter);
				} else {
					continue;
				}
			}
		}

		// Direction vector inside face
		if ((dAngleLA < dAngleLR) &&
			(dAngleRA < dAngleLR)
		) {
			return (*iter);
		}
	}

	_EXCEPTIONT("Logic Error: No exit Face found from Node");
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesFuzzy::FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const FindFaceStruct & aFindFaceStruct
) {
	static const Real Tolerance = ReferenceTolerance;

	const std::vector<int> & vecPossibleFaces =
		aFindFaceStruct.vecFaceIndices;

	// Check that there is actually a decision to be made
	if (vecPossibleFaces.size() < 2) {
		_EXCEPTIONT("vecPossibleFaces must contain at least two Faces.");

	// Node is on an edge
	} else if (aFindFaceStruct.loc == Face::NodeLocation_Edge) {

		// Get Faces and Edges
		if (aFindFaceStruct.vecFaceIndices.size() != 2) {
			_EXCEPTIONT("Logic failure");
		}

		const Face & face0 = mesh.faces[aFindFaceStruct.vecFaceIndices[0]];
		const Face & face1 = mesh.faces[aFindFaceStruct.vecFaceIndices[1]];

		const Edge & edge0 = face0.edges[aFindFaceStruct.vecFaceLocations[0]];
		const Edge & edge1 = face1.edges[aFindFaceStruct.vecFaceLocations[1]];

		const Node & node0 = mesh.nodes[edge0[0]];
		const Node & node1 = mesh.nodes[edge0[1]];

		if (edge0 != edge1) {
			_EXCEPTIONT("Logic failure");
		}

		// Calculate all intersections between edges
		std::vector<Node> nodeIntersections;

		bool fCoincident =
			CalculateEdgeIntersections(
				node0, node1, edge0.type,
				nodeBegin, nodeEnd, edgetype,
				nodeIntersections);
/*
		node0.Print("n0");
		node1.Print("n1");
		nodeBegin.Print("nB");
		nodeEnd.Print("nE");

		std::cout << fCoincident << std::endl;
*/
		// If edges are coincident check orientation to determine face 
		if (fCoincident) {

			// Calculate orientation
			Real dOrientation =
				DotProduct((node1 - node0), (nodeEnd - nodeBegin));

			if (dOrientation == 0.0) {
				_EXCEPTIONT("Logic error");

			// Oriented along with face0
			} else if (dOrientation > 0.0) {
				return aFindFaceStruct.vecFaceIndices[0];

			// Oriented along with face1
			} else {
				return aFindFaceStruct.vecFaceIndices[1];
			}
		}

#ifdef VERBOSE
		printf("Nodes: %i %i\n", edge0[0], edge0[1]);
		printf("Intersects: %lu\n", nodeIntersections.size());
		printf("type: %i %i\n", edge0.type, edgetype);
		node0.Print("node0");
		node1.Print("node1");
		nodeBegin.Print("nodeB");
		nodeEnd.Print("nodeE");
#endif

		if (nodeIntersections.size() == 0) {
			_EXCEPTIONT("Logic failure");
		}

		// Both Edges are lines of constant latitude
		if ((edge0.type == Edge::Type_ConstantLatitude) &&
			(edgetype == Edge::Type_ConstantLatitude)
		) {
			_EXCEPTIONT("Non-coincident lines of constant latitude found.");
		}

		// Great circle arc originating on a great circle arc
		if ((edge0.type == Edge::Type_GreatCircleArc) &&
			(edgetype == Edge::Type_GreatCircleArc)
		) {
			// Outward-pointing normal to first great circle arc
			Node nodeNormal(CrossProduct(node0, node1));

			// Direction along second great circle arc projected onto
			// surface of the sphere at nodeBegin
			Node nodeDir;

			GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDir);

			// Dot product to determine direction of motion
			Real dDot = DotProduct(nodeNormal, nodeDir);

/*
			printf("Faces: %i %i\n", aFindFaceStruct.vecFaceIndices[0], aFindFaceStruct.vecFaceIndices[1]);
			printf("Norm: %1.5e %1.5e %1.5e\n", nodeNormal.x, nodeNormal.y, nodeNormal.z);
			printf("Dir: %1.5e %1.5e %1.5e\n", nodeDir.x, nodeDir.y, nodeDir.z);
			printf("%1.5e\n", dDot);
*/
			if (dDot > Tolerance) {
				return aFindFaceStruct.vecFaceIndices[0];

			} else if (dDot < -Tolerance) {
				return aFindFaceStruct.vecFaceIndices[1];

			} else {
				_EXCEPTIONT("Logic error");
			}
		}

		// Great circle arc originating on a line of constant latitude
		if ((edge0.type == Edge::Type_ConstantLatitude) &&
			(edgetype == Edge::Type_GreatCircleArc)
		) {
			// Check for equatorial failure
			if ((fabs(nodeBegin.z) < Tolerance) &&
				(nodeIntersections.size() > 1)
			) {
				_EXCEPTIONT("Equatorial failure in intersection");
			}

			// Find top Face and bottom Face
			int ixTopFace;
			int ixBottomFace;

			// Orientation of edge indicates which Face is on top
			if (IsPositivelyOrientedEdge(node0, node1)) {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[0];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[1];
			} else {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[1];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[0];
			}
/*
			printf("Pos: %i\n" ,IsPositivelyOrientedEdge(node0, node1));
			printf("Face: %i %i\n", aFindFaceStruct.vecFaceIndices[0], aFindFaceStruct.vecFaceIndices[1]);
			printf("nodeInt: %i\n", nodeIntersections.size());
*/
			// Only begin node intersects
			if (nodeIntersections.size() == 1) {
				if (nodeEnd.z > node0.z) {
					return ixTopFace;
				} else if (nodeEnd.z < node0.z) {
					return ixBottomFace;
				} else if (nodeBegin.z >= Tolerance) {
					return ixTopFace;
				} else if (nodeBegin.z <= Tolerance) {
					return ixBottomFace;
				} else {
					_EXCEPTIONT("Logic error");
				}

			// Great circle arc intersects twice
			} else if (nodeIntersections.size() == 2) {
				if (nodeBegin.z >= Tolerance) {
					return ixTopFace;
				} else if (nodeBegin.z <= Tolerance) {
					return ixBottomFace;
				} else {
					_EXCEPTIONT("Logic error");
				}

			} else {
				_EXCEPTIONT("Too many intersections detected");
			}
		}

		// Line of constant latitude originating on a great circle arc
		if ((edgetype == Edge::Type_ConstantLatitude) &&
			(edge0.type == Edge::Type_GreatCircleArc)
		) {
			// Exterior pointing normal to great circle arc
			Node nodeNormal(CrossProduct(node1, node0));

			// Direction of line of constant latitude
			Node nodeDir;

			GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDir);

			// Dot product
			Real dDot = DotProduct(nodeDir, nodeNormal);
/*
			printf("F: %i %i\n",
				aFindFaceStruct.vecFaceIndices[0],
				aFindFaceStruct.vecFaceIndices[1]);
			printf("E0: %i %i\n", edge0[0], edge0[1]);
			printf("E1: %i %i\n", edge1[0], edge1[1]);
			nodeNormal.Print("Norm");
			nodeDir.Print("Dir ");
			printf("Dot: %1.5e\n", dDot);

			printf("%i %i\n",
				aFindFaceStruct.vecFaceIndices[0],
				aFindFaceStruct.vecFaceIndices[1]);
*/
			if (dDot > Tolerance) {
				return aFindFaceStruct.vecFaceIndices[1];
			} else if (dDot < - Tolerance) {
				return aFindFaceStruct.vecFaceIndices[0];
			}

			// Zero dot product (quadratic intersection)
			int ixTopFace;
			int ixBottomFace;

			// Orientation of edge indicates which Face is on top
			if (IsPositivelyOrientedEdge(node0, node1)) {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[0];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[1];
			} else {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[1];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[0];
			}

			if (nodeBegin.z > Tolerance) {
				return ixTopFace;

			} else if (nodeBegin.z < Tolerance) {
				return ixBottomFace;

			} else {
				_EXCEPTIONT("Logic error");
			}
		}

	// Node is at a corner
	} else if (aFindFaceStruct.loc == Face::NodeLocation_Corner) {

		if (aFindFaceStruct.vecFaceLocations.size() < 3) {
			_EXCEPTIONT("Logic error");
		}

		int ixLocation = aFindFaceStruct.vecFaceLocations[0];

		int ixFace = aFindFaceStruct.vecFaceIndices[0];

		int ixNode = mesh.faces[ixFace][ixLocation];

		return FindFaceNearNode(
			mesh,
			ixNode,
			nodeEnd,
			edgetype);

	}

	// This function does not handle Interior or Exterior nodes
	_EXCEPTIONT("Invalid Node location");
}

///////////////////////////////////////////////////////////////////////////////

