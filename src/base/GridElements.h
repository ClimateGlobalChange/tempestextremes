///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.h
///	\author  Paul Ullrich
///	\version March 7, 2014
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

#ifndef _GRIDELEMENTS_H_
#define _GRIDELEMENTS_H_

///////////////////////////////////////////////////////////////////////////////

#include "Defines.h"

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <cassert>

#if defined(OVERLAPMESH_USE_UNSORTED_MAP)
#include <unordered_map>
#endif
#if defined(OVERLAPMESH_USE_NODE_MULTIMAP)
#include "node_multimap_3d.h"
#endif

#include "Exception.h"
#include "DataArray1D.h"
#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single point in 3D Cartesian geometry.
///	</summary>
class Node {

public:
	///	<summary>
	///		Cartesian coordinates (x,y,z) of this Node.
	///	</summary>
	Real x;
	Real y;
	Real z;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Node() :
		x(0.0),
		y(0.0),
		z(0.0)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	Node(
		Real _x,
		Real _y,
		Real _z
	) :
		x(_x),
		y(_y),
		z(_z)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Node(const Node & node) {
		x = node.x;
		y = node.y;
		z = node.z;
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	const Node & operator=(const Node & node) {
		x = node.x;
		y = node.y;
		z = node.z;

		return (*this);
	}

	///	<summary>
	///		Comparator operator using floating point tolerance.
	///	</summary>
	bool operator< (const Node & node) const {
		static const Real Tolerance = ReferenceTolerance;

		if (x - node.x <= -Tolerance) {
			return true;
		} else if (x - node.x >= Tolerance) {
			return false;
		}

		if (y - node.y <= -Tolerance) {
			return true;
		} else if (y - node.y >= Tolerance) {
			return false;
		}

		if (z - node.z <= -Tolerance) {
			return true;
		} else if (z - node.z >= Tolerance) {
			return false;
		}

		return false;
	}

	///	<summary>
	///		Comparator operator using floating point tolerance.
	///	</summary>
	bool operator== (const Node & node) const {
		static const Real Tolerance = ReferenceTolerance;

		if (fabs(x - node.x) >= Tolerance) {
			return false;
		}
		if (fabs(y - node.y) >= Tolerance) {
			return false;
		}
		if (fabs(z - node.z) >= Tolerance) {
			return false;
		}

		return true;
	}

	///	<summary>
	///		Comparator operator using floating point tolerance.
	///	</summary>
	bool operator!= (const Node & node) const {
		return (!((*this) == node));
	}

	///	<summary>
	///		Difference between two nodes.
	///	</summary>
	Node operator-(const Node & node) const {
		Node nodeDiff;
		nodeDiff.x = x - node.x;
		nodeDiff.y = y - node.y;
		nodeDiff.z = z - node.z;

		return nodeDiff;
	}

	///	<summary>
	///		Sum of two nodes.
	///	</summary>
	Node operator+(const Node & node) const {
		Node nodeSum;
		nodeSum.x = x + node.x;
		nodeSum.y = y + node.y;
		nodeSum.z = z + node.z;

		return nodeSum;
	}

	///	<summary>
	///		muiltiple node position-vector by a constant
	///	</summary>
	Node operator*(const double c) const {
		Node result;
		result.x = x*c;
		result.y = y*c;
		result.z = z*c;
		return result;
	}

	///	<summary>
	///		divide node position-vector by a constant
	///	</summary>
	Node operator/(const double c) const {
		Node result;
		result.x = x/c;
		result.y = y/c;
		result.z = z/c;
		return result;
	}

	///	<summary>
	///		Project node onto the unit sphere
	///	</summary>
	Node Normalized() const {
		return (*this)/Magnitude();;
	}

	///	<summary>
	///		Magnitude of this node.
	///	</summary>
	Real Magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}

	///	<summary>
	///		Output node to stdout.
	///	</summary>
	void Print(const char * szName) const {
		printf("%s: %1.15e %1.15e %1.15e\n", szName, x, y, z);
	}
};


///	<summary>
///		A vector for the storage of Nodes.
///	</summary>
typedef std::vector<Node> NodeVector;

///	<summary>
///		A map between Nodes and indices.
///	</summary>
#if defined(OVERLAPMESH_RETAIN_REPEATED_NODES)
typedef std::map<Node, int> NodeMap;
#endif

#if defined(OVERLAPMESH_USE_NODE_MULTIMAP)
typedef node_multimap_3d<Node, int> NodeMap;
#endif

#if defined(OVERLAPMESH_USE_UNSORTED_MAP)
///	<summary>
///		Hasher for a Node.
///	</summary>
struct NodeHash {
	std::size_t operator()(const Node& node) const {
		static const double bin_width = OVERLAPMESH_BIN_WIDTH;

		int i = static_cast<int>((node.x + 2.123456789101112) / bin_width);
		int j = static_cast<int>((node.y + 2.123456789101112) / bin_width);
		int k = static_cast<int>((node.z + 2.123456789101112) / bin_width);

		return (std::size_t)(i * 18397 + j * 20483 + k * 29303);
	}
};

typedef std::unordered_map<Node, int, NodeHash> NodeMap;
#endif

///	<summary>
///		Value type for NodeMap.
///	</summary>
typedef NodeMap::value_type NodeMapPair;

///	<summary>
///		Constant iterator for NodeMap.
///	</summary>
typedef NodeMap::const_iterator NodeMapConstIterator;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A node index.
///	</summary>
typedef int NodeIndex;

///	<summary>
///		A vector for the storage of Node indices.
///	</summary>
typedef std::vector<NodeIndex> NodeIndexVector;

///	<summary>
///		An index indicating this Node is invalid.
///	</summary>
static const NodeIndex InvalidNode = (-1);

///	<summary>
///		An index indicating this Face is invalid.
///	</summary>
static const NodeIndex InvalidFace = (-1);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An edge connects two nodes.
///	</summary>
class Edge {

public:
	///	<summary>
	///		Type of edge.
	///	</summary>
	enum Type {
		Type_GreatCircleArc = 0,
		Type_Default = Type_GreatCircleArc,
		Type_ConstantLatitude = 1
	};

public:
	///	<summary>
	///		Node indices representing the endpoints of this edge.
	///	</summary>
	int node[2];

	///	<summary>
	///		The type of this edge.
	///	</summary>
	Type type;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Edge(
		int node0 = InvalidNode,
		int node1 = InvalidNode,
		Type _type = Type_Default
	) {
		node[0] = node0;
		node[1] = node1;
		type = _type;
	}

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Edge()
	{ }

	///	<summary>
	///		Flip the order of the nodes stored in the segment.  Note that this
	///		does not affect the comparator properties of the segment, and so
	///		this method can be treated as const.
	///	</summary>
	void Flip() const {
		int ixTemp = node[0];
		const_cast<int&>(node[0]) = node[1];
		const_cast<int&>(node[1]) = ixTemp;
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return node[i];
	}

	int & operator[](int i) {
		return node[i];
	}

	///	<summary>
	///		Get the nodes as an ordered pair.
	///	</summary>
	void GetOrderedNodes(
		int & ixNodeSmall,
		int & ixNodeBig
	) const {
		if (node[0] < node[1]) {
			ixNodeSmall = node[0];
			ixNodeBig   = node[1];
		} else {
			ixNodeSmall = node[1];
			ixNodeBig   = node[0];
		}
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Edge & edge) const {

		// Order the local nodes
		int ixNodeSmall;
		int ixNodeBig;
		GetOrderedNodes(ixNodeSmall, ixNodeBig);

		// Order the nodes in edge
		int ixEdgeNodeSmall;
		int ixEdgeNodeBig;
		edge.GetOrderedNodes(ixEdgeNodeSmall, ixEdgeNodeBig);

		// Compare
		if (ixNodeSmall < ixEdgeNodeSmall) {
			return true;
		} else if (ixNodeSmall > ixEdgeNodeSmall) {
			return false;
		} else if (ixNodeBig < ixEdgeNodeBig) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Edge & edge) const {
		if (edge.type != type) {
			return false;
		}

		if ((node[0] == edge.node[0]) &&
			(node[1] == edge.node[1])
		) {
			return true;

		} else if (
			(node[0] == edge.node[1]) &&
			(node[1] == edge.node[0])
		) {
			return true;
		}

		return false;
	}

	///	<summary>
	///		Inequality operator.
	///	</summary>
	bool operator!=(const Edge & edge) const {
		return !((*this) == edge);
	}

	///	<summary>
	///		Return the node that is shared between segments.
	///	</summary>
	int CommonNode(
		const Edge & edge
	) const {
		if (edge[0] == node[0]) {
			return node[0];
		} else if (edge[0] == node[1]) {
			return node[1];
		} else if (edge[1] == node[0]) {
			return node[0];
		} else if (edge[1] == node[1]) {
			return node[1];
		} else {
			return InvalidNode;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An edge connects two nodes with a sub-array of interior nodes.
///	</summary>
class MultiEdge : public std::vector<int> {

public:
	///	<summary>
	///		Flip the edge.
	///	</summary>
	MultiEdge Flip() const {
		MultiEdge edgeFlip;
		for (int i = size()-1; i >= 0; i--) {
			edgeFlip.push_back((*this)[i]);
		}
		return edgeFlip;
	}

};

typedef std::vector<MultiEdge> MultiEdgeVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A pair of face indices, typically on opposite sides of an Edge.
///	</summary>
class FacePair {

public:
	///	<summary>
	///		Indices of the Faces in this pair.
	///	</summary>
	int face[2];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FacePair() {
		face[0] = InvalidFace;
		face[1] = InvalidFace;
	}

	///	<summary>
	///		Add a face to this FacePair.
	///	</summary>
	void AddFace(int ixFace) {
		if (face[0] == InvalidFace) {
			face[0] = ixFace;

		} else if (face[1] == InvalidFace) {
			face[1] = ixFace;

		} else {
			_EXCEPTIONT("FacePair already has a full set of Faces.");
		}
	}

	///	<summary>
	///		Does this FacePair have a complete set of Faces?
	///	</summary>
	bool IsComplete() const {
		return ((face[0] != InvalidFace) && (face[1] != InvalidFace));
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return face[i];
	}
};

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<Edge> EdgeVector;

typedef std::map<Edge, FacePair> EdgeMap;

typedef EdgeMap::value_type EdgeMapPair;

typedef EdgeMap::iterator EdgeMapIterator;

typedef EdgeMap::const_iterator EdgeMapConstIterator;

typedef std::vector<EdgeMap::iterator> EdgeMapIteratorVector;

typedef std::set<Edge> EdgeSet;

typedef std::pair<Edge, FacePair> EdgePair;

typedef std::vector<EdgePair> EdgeMapVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A face.
///	</summary>
class Face {

public:
	///	<summary>
	///		Vector of node indices bounding this face, stored in
	///		counter-clockwise order.
	///	</summary>
	EdgeVector edges;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Face(
		int edge_count = 0
	) {
		edges.resize(edge_count);
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	inline int operator[](int ix) const {
		return edges[ix][0];
	}

	///	<summary>
	///		Set a node.
	///	</summary>
	void SetNode(int ixLocal, int ixNode) {
		int nEdges = static_cast<int>(edges.size());
		edges[ixLocal][0] = ixNode;

		int ixPrev = (ixLocal + nEdges - 1) % nEdges;
		edges[ixPrev][1] = ixNode;
	}

public:
	///	<summary>
	///		Possible locations of nodes.
	///	</summary>
	enum NodeLocation {
		NodeLocation_Undefined = (-1),
		NodeLocation_Exterior = 0,
		NodeLocation_Default = NodeLocation_Exterior,
		NodeLocation_Interior = 1,
		NodeLocation_Edge = 2,
		NodeLocation_Corner = 3
	};

	///	<summary>
	///		Determine the Edge index corresponding to the given Edge.  If the
	///		Edge is not found an Exception is thrown.
	///	</summary>
	int GetEdgeIndex(const Edge & edge) const;

	///	<summary>
	///		Remove zero Edges (Edges with repeated Node indices)
	///	</summary>
	void RemoveZeroEdges();
};

///	<summary>
///		A vector of Faces.
///	</summary>
typedef std::vector<Face> FaceVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A reverse node array stores all faces associated with a given node.
///	</summary>
typedef std::vector< std::set<int> > ReverseNodeArray;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A mesh.
///	</summary>
class Mesh {

public:
    ///	<summary>
    ///		Type of standard mesh
    ///	</summary>
    enum MeshType {
        MeshType_Unknown = -1,
        MeshType_CubedSphere = 0,
        MeshType_RLL = 1,
        MeshType_IcosaHedral = 2,
        MeshType_IcosaHedralDual = 3,
        MeshType_Overlap = 4,
		MeshType_UTM = 5
    };

    MeshType type;

public:
	///	<summary>
	///		Filename for this mesh.
	///	</summary>
	std::string strFileName;

	///	<summary>
	///		Vector of Nodes for this mesh.
	///	</summary>
	NodeVector nodes;

	///	<summary>
	///		Vector of Faces for this mesh.
	///	<summary>
	FaceVector faces;

	///	<summary>
	///		Vector of first mesh Face indices.
	///	</summary>
	std::vector<int> vecSourceFaceIx;

	///	<summary>
	///		Vector of second mesh Face indices.
	///	</summary>
	std::vector<int> vecTargetFaceIx;

	///	<summary>
	///		Vector of Face areas.
	///	</summary>
	DataArray1D<double> vecFaceArea;

	///	<summary>
	///		Vector storing mask variable for this mesh.
	///	</summary>
	DataArray1D<int> vecMask;

	///	<summary>
	///		EdgeMap for this mesh.
	///	</summary>
	EdgeMap edgemap;

	///	<summary>
	///		ReverseNodeArray for this mesh.
	///	</summary>
	ReverseNodeArray revnodearray;

	///	<summary>
	///		Indices of the original Faces for this mesh (for use when
	///		the original mesh has been subdivided).
	///	</summary>
	std::vector<int> vecMultiFaceMap;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
    Mesh() : type(MeshType_Unknown) {
	}

	///	<summary>
	///		Constructor with input mesh parameter.
	///	</summary>
	Mesh(const std::string & strFile) {
		Read(strFile);
	}

public:
	///	<summary>
	///		Clear the contents of the mesh.
	///	</summary>
	void Clear();

	///	<summary>
	///		Construct the EdgeMap from the NodeVector and FaceVector.
	///	</summary>
	void ConstructEdgeMap();

	///	<summary>
	///		Construct the ReverseNodeArray from the NodeVector and FaceVector.
	///	</summary>
	void ConstructReverseNodeArray();

	///	<summary>
	///		Calculate Face areas.
	///	</summary>
	Real CalculateFaceAreas(
		bool fContainsConcaveFaces
	);

	///	<summary>
	///		Calculate Face areas from an Overlap mesh.
	///	</summary>
	Real CalculateFaceAreasFromOverlap(
		const Mesh & meshOverlap
	);

	///	<summary>
	///		Sort Faces by the opposite source mesh.
	///	</summary>
	void ExchangeFirstAndSecondMesh();

	///	<summary>
	///		Remove coincident nodes from the Mesh and adjust indices in faces.
	///	</summary>
	void RemoveCoincidentNodes();

	///	<summary>
	///		Write the mesh to a NetCDF file in Exodus format.
	///	</summary>
	void Write(
		const std::string & strFile,
		NcFile::FileFormat eFileFormat = NcFile::Classic
	) const;

	///	<summary>
	///		Write the mesh to a NetCDF file in SCRIP format.
	///	</summary>
	void WriteScrip(
		const std::string & strFile,
		NcFile::FileFormat eFileFormat = NcFile::Classic
	) const;

	///	<summary>
	///		Read the mesh from a NetCDF file.
	///	</summary>
	void Read(const std::string & strFile);

	///	<summary>
	///		Remove zero edges from all Faces.
	///	</summary>
	void RemoveZeroEdges();

	///	<summary>
	///		Validate the Mesh.
	///	</summary>
	void Validate() const;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Location data returned from FindFaceFromNode()
///		Generate a PathSegmentVector describing the path around the face
///		ixCurrentFirstFace.
///	</summary>
struct FindFaceStruct {
	
	///	<summary>
	///		A vector of face indices indicating possible Faces.
	///	</summary>
	std::vector<int> vecFaceIndices;

	///	<summary>
	///		A vector of locations on each Face.  If loc is NodeLocation_Corner,
	///		this corresponds to the associated corner of the Face.  If loc
	///		is NodeLocation_Edge, this corresponds to the associated Edge of
	///		the Face.  If loc is NodeLocation_Interior, this value is
	///		undefined.
	///	</summary>
	std::vector<int> vecFaceLocations;

	///	<summary>
	///		The NodeLocation where this Node lies.
	///	</summary>
	Face::NodeLocation loc;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate latitude and longitude from normalized 3D Cartesian
///		coordinates, in degrees.
///	</summary>
inline void XYZtoRLL_Deg(
	const double & dX,
	const double & dY,
	const double & dZ,
	double & dLon,
	double & dLat
) {
	assert(fabs(dX * dX + dY * dY + dZ * dZ - 1.0) < HighTolerance);

	if (fabs(dZ) < 1.0 - ReferenceTolerance) {
		dLon = atan2(dY, dX);
		dLat = asin(dZ);

		if (dLon < 0.0) {
			dLon += 2.0 * M_PI;
		}

		dLon = dLon / M_PI * 180.0;
		dLat = dLat / M_PI * 180.0;

	} else if (dZ > 0.0) {
		dLon = 0.0;
		dLat = 90.0;

	} else {
		dLon = 0.0;
		dLat = -90.0;
	}
}

///	<summary>
///		Calculate an average longitude from two given longitudes (in radians).
///	</summary>
inline double AverageLongitude_Rad(
	double dLon1,
	double dLon2
) {
	double dDeltaLon;
	if (dLon2 > dLon1) {
		dDeltaLon = fmod(dLon2 - dLon1, 2.0 * M_PI);
		if (dDeltaLon > M_PI) {
			dDeltaLon = dDeltaLon - 2.0 * M_PI;
		}
	} else {
		dDeltaLon = - fmod(dLon1 - dLon2, 2.0 * M_PI);
		if (dDeltaLon < -M_PI) {
			dDeltaLon = dDeltaLon + 2.0 * M_PI;
		}
	}

	double dLonAvg = dLon1 + 0.5 * dDeltaLon;

	if ((dLonAvg < 0.0) && (dLon1 >= 0.0) && (dLon2 >= 0.0)) {
		dLonAvg += 2.0 * M_PI;
	}
	if ((dLonAvg > 2.0 * M_PI) && (dLon1 <= 2.0 * M_PI) && (dLon2 <= 2.0 * M_PI)) {
		dLonAvg -= 2.0 * M_PI;
	}

	return dLonAvg;
}

///	<summary>
///		Calculate the great circle distance between two RLL points.
///	</summary>
inline double GreatCircleDistance_Deg(
	double dLon1,
	double dLat1,
	double dLon2,
	double dLat2
) {
	double dR =
		sin(dLat1) * sin(dLat2)
		+ cos(dLat1) * cos(dLat2) * cos(dLon2 - dLon1);

	if (dR >= 1.0) {
		dR = 0.0;
	} else if (dR <= -1.0) {
		dR = 180.0;
	} else {
		dR = 180.0 / M_PI * acos(dR);
	}
	if (dR != dR) {
		_EXCEPTIONT("NaN value detected");
	}

	return dR;
}

///	<summary>
///		Calculate the dot product between two Nodes.
///	</summary>
inline Real DotProduct(
	const Node & node1,
	const Node & node2
) {
	return (node1.x * node2.x + node1.y * node2.y + node1.z * node2.z);
}

///	<summary>
///		Calculate the cross product between two Nodes.
///	</summary>
inline Node CrossProduct(
	const Node & node1,
	const Node & node2
) {
	Node nodeCross;
	nodeCross.x = node1.y * node2.z - node1.z * node2.y;
	nodeCross.y = node1.z * node2.x - node1.x * node2.z;
	nodeCross.z = node1.x * node2.y - node1.y * node2.x;

	return nodeCross;
}

///	<summary>
///		Calculate the product of a Node with a scalar.
///	</summary>
inline Node ScalarProduct(
	const Real & d,
	const Node & node
) {
	Node nodeProduct(node);
	nodeProduct.x *= d;
	nodeProduct.y *= d;
	nodeProduct.z *= d;
	return nodeProduct;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if an edge is positively oriented
///		(aligned with increasing longitude).
///	</summary>
bool IsPositivelyOrientedEdge(
	const Node & nodeBegin,
	const Node & nodeEnd
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the local direction vector along the surface of the sphere
///		for the given edge.
///	</summary>
void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Node & nodeRef,
	const Edge::Type edgetype,
	Node & nodeDir
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the local direction vector along the surface of the sphere
///		for the given edge.
///	</summary>
void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	Node & nodeDir
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		For all Nodes on meshSecond that are "nearby" a Node on meshFirst,
///		set the Node equal to the meshFirst Node.
///	</summary>
void EqualizeCoincidentNodes(
	const Mesh & meshFirst,
	Mesh & meshSecond
);

///////////////////////////////////////////////////////////////////////////////
/*
///	<summary>
///		Equate coincident nodes on mesh.
///	</summary>
void EqualizeCoincidentNodes(
	Mesh & mesh
);
*/
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build the mapping function for nodes on meshSecond which are
///		coincident with nodes on meshFirst.
///	</summary>
///	<returns>
///		The number of coincident nodes on meshSecond.
///	</returns>
int BuildCoincidentNodeVector(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	std::vector<int> & vecSecondToFirstCoincident
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Check if the specified Face is concave.
///	</summary>
bool IsFaceConcave(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a single Face, detecting concave Faces.
///	</summary>
Real CalculateFaceArea_Concave(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a single Face.
///	</summary>
Real CalculateFaceArea(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a concave Face into a Convex face.  Based on routine by
///		Mark Bayazit (https://mpen.ca/406/bayazit).
///	</summary>
///	<returns>
///		true if the Face is convex and has been removed from the mesh.faces
///		vector.
///	</returns>
bool ConvexifyFace(
	Mesh & meshin,
	Mesh & meshout,
	int iFace,
	bool fRemoveConcaveFaces,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert all concave Faces into the Mesh into Concave faces via
///		subdivision.
///	</summary>
void ConvexifyMesh(
	Mesh & mesh,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert concave Mesh meshin to convex Mesh meshout by dividing
///		Faces and populating the MultiFaceMap.
///	</summary>
void ConvexifyMesh(
	Mesh & meshin,
	Mesh & meshout,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find a node within the specified quadrilateral.
///	</summary>
inline Node InterpolateQuadrilateralNode(
	const Node & node0,
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB
) {
	Node nodeRef;

	nodeRef.x =
		  (1.0 - dA) * (1.0 - dB) * node0.x
		+        dA  * (1.0 - dB) * node1.x
		+        dA  *        dB  * node2.x
		+ (1.0 - dA) *        dB  * node3.x;

	nodeRef.y =
		  (1.0 - dA) * (1.0 - dB) * node0.y
		+        dA  * (1.0 - dB) * node1.y
		+        dA  *        dB  * node2.y
		+ (1.0 - dA) *        dB  * node3.y;

	nodeRef.z =
		  (1.0 - dA) * (1.0 - dB) * node0.z
		+        dA  * (1.0 - dB) * node1.z
		+        dA  *        dB  * node2.z
		+ (1.0 - dA) *        dB  * node3.z;

	double dMag = sqrt(
		  nodeRef.x * nodeRef.x
		+ nodeRef.y * nodeRef.y
		+ nodeRef.z * nodeRef.z);

	nodeRef.x /= dMag;
	nodeRef.y /= dMag;
	nodeRef.z /= dMag;

	return nodeRef;
}

///////////////////////////////////////////////////////////////////////////////

#endif

