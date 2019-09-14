///////////////////////////////////////////////////////////////////////////////
///
///	\file    PolynomialInterp.h
///	\author  Paul Ullrich
///	\version December 19, 2011
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

#ifndef _POLYNOMIALINTERP_H_
#define _POLYNOMIALINTERP_H_

#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for performing polynomial interpolation.
///	</summary>
class PolynomialInterp {

private:
	///	<summary>
	///		Private constructor.
	///	<summary>
	PolynomialInterp()
	{ }

public:
	///	<summary>
	///		Determine the coefficients of a Lagrangian polynomial through the
	///		specified points dX which is sampled at point dXsample.
	///	</summary>
	static void LagrangianPolynomialCoeffs(
		int nPoints,
		const double * dX,
		double * dCoeffs,
		double dXsample
	);

	///	<summary>
	///		Determine the coefficients of the first derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void DiffLagrangianPolynomialCoeffs(
		int nPoints,
		const double * dX,
		double * dCoeffs,
		double dXsample
	);

	///	<summary>
	///		Determine the coefficients of the second derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void DiffDiffLagrangianPolynomialCoeffs(
		int nPoints,
		const double * dX,
		double * dCoeffs,
		double dXsample
	);

	///	<summary>
	///		Determine the coefficients of the third derivative of a Lagrangian
	///		polynomial through the specified points dX which is sampled at
	///		point dXsample.
	///	</summary>
	static void DiffDiffDiffLagrangianPolynomialCoeffs(
		int nPoints,
		const double * dX,
		double * dCoeffs,
		double dXsample
	);

	///	<summary>
	///		Interpolate a polynomial through the given (X,Y) points and sample
	///		at point dXsample.  This method is faster and more computationally
	///		stable than InterpolateCoeffs.
	///	</summary>
	static double Interpolate(
		int nPoints,
		const double * dX,
		const double * dY,
		double dXsample
	);

	///	<summary>
	///		Obtain the coefficients a_i of a polynomial interpolated through the
	///		points (X,Y), where
	///		    p(x) = a_0 + a_1 (x - dXmid) + ... + a_(n-1) (x - dXmid)^(n-1)
	///	</summary>
	///	<parameter name = "dWorkingSpace">
	///		An external workspace of at least of size nPoints^2.
	///	</parameter>
	///	<parameter name = "iPivot">
	///		An external workspace of at least size nPoints.
	///	</parameter>
	static void InterpolateCoeffs(
		int nPoints,
		const double * dX,
		const double * dY,
		double * dA,
		double dXmid = 0.0,
		double * dWorkspace = NULL,
		int * iPivot = NULL
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

