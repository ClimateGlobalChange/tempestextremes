///////////////////////////////////////////////////////////////////////////////
///
///	\file    FiniteElementTools.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

static const double InverseMapTolerance = 1.0e-13;

///////////////////////////////////////////////////////////////////////////////

class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the array of GLL nodal locations within an element.
///	</summary>
void GetDefaultNodalLocations(
	int nP,
	DataArray1D<double> & dG
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Apply the local map.
///	</summary>
void ApplyLocalMap(
	const Face & face,
	const NodeVector & nodes,
	double dAlpha,
	double dBeta,
	Node & node
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Apply the local map.
///	</summary>
void ApplyLocalMap(
	const Face & face,
	const NodeVector & nodes,
	double dAlpha,
	double dBeta,
	Node & node,
	Node & dDx1G,
	Node & dDx2G
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Apply inverse map using Newton's method.
///	</summary>
void ApplyInverseMap(
	const Face & face,
	const NodeVector & nodes,
	const Node & node,
	double & dAlpha,
	double & dBeta
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate Mesh meta data for a spectral element grid.
///	</summary>
double GenerateMetaData(
	const Mesh & mesh,
	int nP,
	bool fBubble,
	DataArray3D<int> & dataGLLnodes,
	DataArray3D<double> & dataGLLJacobian
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate unique Jacobian values from non-unique Jacobians.
///	</summary>
void GenerateUniqueJacobian(
	const DataArray3D<int> & dataGLLnodes,
	const DataArray3D<double> & dataGLLJacobian,
	DataArray1D<double> & dataUniqueJacobian
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate Jacobian vector from Jacobians on GLL nodes.
///	</summary>
void GenerateDiscontinuousJacobian(
	const DataArray3D<double> & dataGLLJacobian,
	DataArray1D<double> & dataUniqueJacobian
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the coefficients for sampling a 2D finite element at the
///		specified point.
///	</summary>
void SampleGLLFiniteElement(
	int nMonotoneType,
	int nP,
	double dAlpha,
	double dBeta,
	DataArray2D<double> & dCoeff
);

///////////////////////////////////////////////////////////////////////////////

