///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussLobattoQuadrature.cpp
///	\author  Paul Ullrich
///	\version January 1, 2015
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

#include "GaussLobattoQuadrature.h"

#include "LegendrePolynomial.h"

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

void GaussLobattoQuadrature::GetPoints(
	int nCount,
	DataArray1D<double> & dG,
	DataArray1D<double> & dW
) {
	// Check for valid range
	if (nCount < 2) {
		_EXCEPTION1("Invalid count (%i): Minimum count 2", nCount);
	}

	// Initialize the arrays
	dG.Allocate(nCount);
	dW.Allocate(nCount);

	// Degree 2
	if (nCount == 2) {
		dG[0] = -1.0;
		dG[1] = +1.0;

		dW[0] = +1.0;
		dW[1] = +1.0;

	// Degree 3
	} else if (nCount == 3) {
		dG[0] = -1.0;
		dG[1] =  0.0;
		dG[2] = +1.0;

		dW[0] = +0.333333333333333;
		dW[1] = +1.333333333333334;
		dW[2] = +0.333333333333333;

	// Degree 4
	} else if (nCount == 4) {
		dG[0] = -1.0;
		dG[1] = -0.447213595499958;
		dG[2] = +0.447213595499958;
		dG[3] = +1.0;

		dW[0] = 0.166666666666667;
		dW[1] = 0.833333333333333;
		dW[2] = 0.833333333333333;
		dW[3] = 0.166666666666667;

	// Degree 5
	} else if (nCount == 5) {
		dG[0] = -1.0;
		dG[1] = -0.654653670707977;
		dG[2] =  0.0;
		dG[3] = +0.654653670707977;
		dG[4] = +1.0;

		dW[0] = 0.100000000000000;
		dW[1] = 0.544444444444445;
		dW[2] = 0.711111111111110;
		dW[3] = 0.544444444444445;
		dW[4] = 0.100000000000000;

	// Degree 6
	} else if (nCount == 6) {
		dG[0] = -1.0;
		dG[1] = -0.765055323929465;
		dG[2] = -0.285231516480645;
		dG[3] = +0.285231516480645;
		dG[4] = +0.765055323929465;
		dG[5] = +1.0;

		dW[0] = 0.066666666666667;
		dW[1] = 0.378474956297847;
		dW[2] = 0.554858377035486;
		dW[3] = 0.554858377035486;
		dW[4] = 0.378474956297847;
		dW[5] = 0.066666666666667;

	// Degree 7
	} else if (nCount == 7) {
		dG[0] = -1.0;
		dG[1] = -0.830223896278567;
		dG[2] = -0.468848793470714;
		dG[3] =  0.0;
		dG[4] = +0.468848793470714;
		dG[5] = +0.830223896278567;
		dG[6] = +1.0;

		dW[0] = 0.047619047619048;
		dW[1] = 0.276826047361566;
		dW[2] = 0.431745381209862;
		dW[3] = 0.487619047619048;
		dW[4] = 0.431745381209862;
		dW[5] = 0.276826047361566;
		dW[6] = 0.047619047619048;

	// Degree 8
	} else if (nCount == 8) {
		dG[0] = -1.0;
		dG[1] = -0.871740148509607;
		dG[2] = -0.591700181433142;
		dG[3] = -0.209299217902479;
		dG[4] = +0.209299217902479;
		dG[5] = +0.591700181433142;
		dG[6] = +0.871740148509607;
		dG[7] = +1.0;

		dW[0] = 0.035714285714286;
		dW[1] = 0.210704227143506;
		dW[2] = 0.341122692483505;
		dW[3] = 0.412458794658703;
		dW[4] = 0.412458794658703;
		dW[5] = 0.341122692483505;
		dW[6] = 0.210704227143506;
		dW[7] = 0.035714285714286;

	// Degree 9
	} else if (nCount == 9) {
		dG[0] = -1.0;
		dG[1] = -0.899757995411460;
		dG[2] = -0.677186279510738;
		dG[3] = -0.363117463826178;
		dG[4] =  0.0;
		dG[5] = +0.363117463826178;
		dG[6] = +0.677186279510738;
		dG[7] = +0.899757995411460;
		dG[8] = +1.0;

		dW[0] = 0.027777777777778;
		dW[1] = 0.165495361560806;
		dW[2] = 0.274538712500162;
		dW[3] = 0.346428510973046;
		dW[4] = 0.371519274376417;
		dW[5] = 0.346428510973046;
		dW[6] = 0.274538712500162;
		dW[7] = 0.165495361560806;
		dW[8] = 0.027777777777778;

	// Degree 10
	} else if (nCount == 10) {
		dG[0] = -1.0;
		dG[1] = -0.919533908166459;
		dG[2] = -0.738773865105505;
		dG[3] = -0.477924949810444;
		dG[4] = -0.165278957666387;
		dG[5] = +0.165278957666387;
		dG[6] = +0.477924949810444;
		dG[7] = +0.738773865105505;
		dG[8] = +0.919533908166459;
		dG[9] = +1.0;

		dW[0] = 0.022222222222222;
		dW[1] = 0.133305990851070;
		dW[2] = 0.224889342063126;
		dW[3] = 0.292042683679684;
		dW[4] = 0.327539761183897;
		dW[5] = 0.327539761183897;
		dW[6] = 0.292042683679684;
		dW[7] = 0.224889342063126;
		dW[8] = 0.133305990851070;
		dW[9] = 0.022222222222222;

	// Higher degrees
	} else {
		DataArray1D<double> dRoots(nCount);
		
		LegendrePolynomial::AllDerivativeRoots(nCount-1, dRoots);

		dG[0] = -1.0;
		memcpy(dG+1, dRoots, (nCount-2) * sizeof(double));
		dG[nCount-1] = +1.0;

		for (int k = 0; k < nCount; k++) {
			double dValue =
				LegendrePolynomial::Evaluate(nCount-1, dG[k]);

			double dDegree = static_cast<double>(nCount);

			dW[k] = 2.0 / (dDegree * (dDegree - 1.0) * dValue * dValue);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GaussLobattoQuadrature::GetPoints(
	int nCount,
	double dXi0,
	double dXi1,
	DataArray1D<double> & dG,
	DataArray1D<double> & dW
) {
	// Get quadrature points in the [-1, 1] reference element
	GetPoints(nCount, dG, dW);

	// Scale quadrature points
	for (int i = 0; i < nCount; i++) {
		dG[i] = dXi0 + 0.5 * (dXi1 - dXi0) * (dG[i] + 1.0);
		dW[i] = 0.5 * (dXi1 - dXi0) * dW[i];
	}
}

///////////////////////////////////////////////////////////////////////////////

