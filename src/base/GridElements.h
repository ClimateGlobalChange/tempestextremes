///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.h
///	\author  Paul Ullrich
///	\version May 12, 2023
///
///	<remarks>
///		Copyright 2023 Paul Ullrich
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

#include "Exception.h"
#include "DataArray1D.h"
#include "CoordTransforms.h"
#include "netcdfcpp.h"
#include "kdtree.h"

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
	///		Set to a given value.
	///	</summary>
	inline const Node & Set(
		Real _x,
		Real _y,
		Real _z
	) {
		x = _x;
		y = _y;
		z = _z;

		return (*this);
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

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// !!!! NOTE: Removed Tolerance from the comparator because there are
		// !!!!       some cases where nodes with exactly equal coordinates
		// !!!!       are not sorted properly in a std::map.
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (x < node.x - ReferenceTolerance) {
			return true;
		} else if (x > node.x + ReferenceTolerance) {
			return false;
		}

		if (y < node.y - ReferenceTolerance) {
			return true;
		} else if (y > node.y + ReferenceTolerance) {
			return false;
		}

		if (z < node.z - ReferenceTolerance) {
			return true;
		} else if (z > node.z + ReferenceTolerance) {
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
	///		Muiltiply by a constant.
	///	</summary>
	Node operator*(const Real c) const {
		Node result;
		result.x = x * c;
		result.y = y * c;
		result.z = z * c;
		return result;
	}

	///	<summary>
	///		Divide by a constant.
	///	</summary>
	Node operator/(const Real c) const {
		Node result;
		result.x = x / c;
		result.y = y / c;
		result.z = z / c;
		return result;
	}

	///	<summary>
	///		Sum of two nodes (in place).
	///	</summary>
	Node & operator+=(const Node & node) {
		x += node.x;
		y += node.y;
		z += node.z;
		return (*this);
	}

	///	<summary>
	///		Difference between two nodes (in place).
	///	</summary>
	Node & operator-=(const Node & node) {
		x -= node.x;
		y -= node.y;
		z -= node.z;
		return (*this);
	}

	///	<summary>
	///		Multiply by a constant (in place).
	///	</summary>
	Node & operator*=(const Real c) {
		x *= c;
		y *= c;
		z *= c;
		return (*this);
	}

	///	<summary>
	///		Divide by a constant (in place).
	///	</summary>
	Node & operator/=(const Real c) {
		x /= c;
		y /= c;
		z /= c;
		return (*this);
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
/*
	///	<summary>
	///		int64 hash of this node (assumes all coordinates in [-1,1] interval).
	///	</summary>
	inline int64_t Hash64() const {
		int64_t i64x = (int64_t)(std::round(0.5 * (1.0 + nodes.x) * (double)(1 << 20)));
		int64_t i64y = (int64_t)(std::round(0.5 * (1.0 + nodes.y) * (double)(1 << 20)));

		return ((i64x << 40) + (i64y << 20) + (i64z));
	}
*/
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
typedef std::map<Node, int> NodeMap;

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
///		A structure for storing Nodes that ensures that minimum Node
///		spacing is met.
///	</summary>
class NodeTree {
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	NodeTree(
		double minimum_spacing = ReferenceTolerance
	);

	///	<summary>
	///		Destructor.
	///	</summary>
	~NodeTree();

	///	<summary>
	///		Find a node in the NodeTree that is within minimum_spacing of
	///		the provided node.  If multiple nodes exist satisfying this criteria
	///		return the node with lowest index.
	///	</summary>
	size_t find(
		const Node & node
	);

	///	<summary>
	///		Similar to find() except insert the node into the NodeTree if no
	///		nearby node is found.  The return value is then equal to index.
	///	</summary>
	size_t find_or_insert(
		const Node & node,
		size_t index
	);

	///	<summary>
	///		Number of nodes in this NodeTree.
	///	</summary>
	size_t size() const {
		return m_size;
	}

private:
	///	<summary>
	///		Minimum node spacing (Cartesian distance).
	///	</summary>
	double m_minimum_spacing;

	///	<summary>
	///		kd-tree used for quick lookup of grid points (optionally initialized).
	///	</summary>
	kdtree * m_kdtree;

	///	<summary>
	///		Number of nodes in tree.
	///	</summary>
	size_t m_size;
};

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
		MeshType_UTM = 5,
		MeshType_Transect = 6
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
	///		Tolerance for removing coincident nodes.
	///	</summary>
	double coincident_node_tolerance;

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
    Mesh(
		double _coincident_node_tolerance = ReferenceTolerance
	) :
		type(MeshType_Unknown),
		coincident_node_tolerance(_coincident_node_tolerance) 
	{ }

	///	<summary>
	///		Constructor with input mesh parameter.
	///	</summary>
	Mesh(
		const std::string & strFile,
		double _coincident_node_tolerance = ReferenceTolerance
	) :
		type(MeshType_Unknown),
		coincident_node_tolerance(_coincident_node_tolerance)
	{
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
	///		Remove coincident nodes from the Mesh and adjust indices in faces
	///		using a kd-tree for nearest neighbor search.
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

