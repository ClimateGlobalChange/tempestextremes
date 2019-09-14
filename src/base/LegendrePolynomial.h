///////////////////////////////////////////////////////////////////////////////
///
///	\file    LegendrePolynomial.h
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#ifndef _LEGENDREPOLYNOMIAL_H_
#define _LEGENDREPOLYNOMIAL_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class for handling operations involving Legendre polynomials.
///	</summary>
class LegendrePolynomial {

private:
	///	<summary>
	///		Private constructor.
	///	<summary>
	LegendrePolynomial()
	{ }

public:
	///	<summary>
	///		Evaluate the Legendre polynomial and its derivative with given
	///		degree at point dX.
	///	</summary>
	static void EvaluateValueAndDerivative(
		int nDegree,
		double dX,
		double & dValue,
		double & dDerivative
	);

	///	<summary>
	///		Evaluate the Legendre polynomial of the given degree at point dX.
	///	</summary>
	static double Evaluate(
		int nDegree,
		double dX
	);

	///	<summary>
	///		Evaluate the derivative of the Legendre polynomial of the given
	///		degree at point dX.
	///	</summary>
	static double EvaluateDerivative(
		int nDegree,
		double dX
	);

	///	<summary>
	///		Determine the number of real roots of the derivative of the
	///		Legendre polynomial of given degree.
	///	</summary>
	static int DerivativeRootCount(
		int nDegree
	) {
		return (nDegree - 1);
	}

	///	<summary>
	///		Retrieve the specified root of the derivative of the Legendre
	///		polynomial of the given degree.
	///	</summary>
	static double DerivativeRoot(
		int nDegree,
		int nRoot
	);

	///	<summary>
	///		Retrieve the specified root of the extended derivative of the 
	///		Legendre polynomial of the given degree.  The extended derivative
	///		is defined as (x^2 - 1) * P'(x), where P'(x) is the usual
	///		derivative.
	///	</summary>
	static double DerivativeExtendedRoot(
		int nDegree,
		int nRoot
	);

public:
	///	<summary>
	///		Evaluate the characteristic function at the given point.
	///	</summary>
	static double EvaluateCharacteristic(
		int nDegree,
		int nRoot,
		double dX
	);

	///	<summary>
	///		Determine the number of real roots of the Legendre polynomial of
	///		the given degree.
	///	</summary>
	static int RootCount(
		int nDegree
	) {
		return nDegree;
	}

	///	<summary>
	///		Return all roots to the Legendre polynomial of the
	///		given degree.
	///	</summary>
	static void AllRoots(
		int nDegree,
		double * dRoots
	);

	///	<summary>
	///		Return all roots to the derivative of the Legendre polynomial of the
	///		given degree.
	///	</summary>
	static void AllDerivativeRoots(
		int nDegree,
		double * dRoots
	);

	///	<summary>
	///		Return the given root to the Legendre polynomial of the
	///		given degree.
	///	</summary>
	static double Root(
		int nDegree,
		int nRoot
	);

};

///////////////////////////////////////////////////////////////////////////////

#endif

