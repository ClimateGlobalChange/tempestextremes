///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussLobattoQuadrature.h
///	\author  Paul Ullrich
///	\version July 9, 2012
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _GAUSSLOBATTOQUADRATURE_H_
#define _GAUSSLOBATTOQUADRATURE_H_

#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Quadrature nodes and weights for Gauss-Lobatto quadrature.
///	</summary>
class GaussLobattoQuadrature {

public:
	///	<summary>
	///		Return the Gauss-Lobatto quadrature points and their corresponding
	///		weights for the given number of points.
	///	</summary>
	static void GetPoints(
		int nCount,
		DataArray1D<double> & dG,
		DataArray1D<double> & dW
	);

	///	<summary>
	///		Retrun the Gauss-Lobatto quadrature points and their corresponding
	///		weights for the given number of points and reference element.
	///	</summary>
	static void GetPoints(
		int nCount,
		double dXi0,
		double dXi1,
		DataArray1D<double> & dG,
		DataArray1D<double> & dW
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

