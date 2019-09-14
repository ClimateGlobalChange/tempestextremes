///////////////////////////////////////////////////////////////////////////////
///
///	\file    PolynomialInterp.cpp
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

#include "PolynomialInterp.h"

#include "Exception.h"

#include <cmath>
#include <cstring>

///////////////////////////////////////////////////////////////////////////////

void PolynomialInterp::LagrangianPolynomialCoeffs(
	int nPoints,
	const double * dX,
	double * dCoeffs,
	double dXsample
) {
	int i;
	int j;

	// Compute Lagrangian polynomial
	for (i = 0; i < nPoints; i++) {
		dCoeffs[i] = 1.0;

		for (j = 0; j < nPoints; j++) {
			if (i == j) {
				continue;
			}

			dCoeffs[i] *= (dXsample - dX[j]) / (dX[i] - dX[j]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void PolynomialInterp::DiffLagrangianPolynomialCoeffs(
	int nPoints,
	const double * dX,
	double * dCoeffs,
	double dXsample
) {
	int i;
	int j;

	int iMatch = (-1);

	// Check if dXsample is equivalent to one of the values of dX
	for (i = 0; i < nPoints; i++) {
		if (fabs(dXsample - dX[i]) < 1.0e-14) {
			iMatch = i;
			break;
		}
	}

	// Match detected
	if (iMatch != (-1)) {

		// Calculate product
		for (i = 0; i < nPoints; i++) {
			dCoeffs[i] = 1.0;

			double dSum = 0.0;

			for (j = 0; j < nPoints; j++) {
				if ((j == i) || (j == iMatch)) {
					continue;
				}

				dCoeffs[i] *= (dXsample - dX[j]) / (dX[i] - dX[j]);

				dSum += 1.0 / (dXsample - dX[j]);
			}

			// Multiply by differential
			if (i != iMatch) {
				dCoeffs[i] *= (1.0 + (dXsample - dX[iMatch]) * dSum)
					/ (dX[i] - dX[iMatch]);
			} else {
				dCoeffs[i] *= dSum;
			}
		}

		return;
	}

	// Compute Lagrangian polynomial coeffs
	LagrangianPolynomialCoeffs(nPoints, dX, dCoeffs, dXsample);

	// Multiply by first moment of differential
	for (i = 0; i < nPoints; i++) {

		double dDifferential = 0.0;
		for (j = 0; j < nPoints; j++) {
			if (i == j) {
				continue;
			}

			dDifferential += 1.0 / (dXsample - dX[j]);
		}

		dCoeffs[i] *= dDifferential;
	}
}

///////////////////////////////////////////////////////////////////////////////

void PolynomialInterp::DiffDiffLagrangianPolynomialCoeffs(
	int nPoints,
	const double * dX,
	double * dCoeffs,
	double dXsample
) {
	int i;
	int j;

	int iMatch = (-1);

	// Check if dXsample is equivalent to one of the values of dX
	for (i = 0; i < nPoints; i++) {
		if (fabs(dXsample - dX[i]) < 1.0e-14) {
			iMatch = i;
			break;
		}
	}

	// Match detected
	if (iMatch != (-1)) {

		// Calculate product
		for (i = 0; i < nPoints; i++) {
			dCoeffs[i] = 1.0;

			double dSum = 0.0;

			double dSum2 = 0.0;

			for (j = 0; j < nPoints; j++) {
				if ((j == i) || (j == iMatch)) {
					continue;
				}

				dCoeffs[i] *= (dXsample - dX[j]) / (dX[i] - dX[j]);

				double dDiff = 1.0 / (dXsample - dX[j]);

				dSum  += dDiff;
				dSum2 += dDiff * dDiff;
			}

			// Multiply by differential
			if (i != iMatch) {
				dCoeffs[i] *= (2.0 * dSum
					+ (dXsample - dX[iMatch]) * (dSum * dSum - dSum2))
					/ (dX[i] - dX[iMatch]);
			} else {
				dCoeffs[i] *= (dSum * dSum - dSum2);
			}
		}

		return;
	}

	// Compute Lagrangian polynomial coeffs
	LagrangianPolynomialCoeffs(nPoints, dX, dCoeffs, dXsample);

	// Multiply by second moment of differential
	for (i = 0; i < nPoints; i++) {

		double dDifferential1 = 0.0;
		double dDifferential2 = 0.0;
		for (j = 0; j < nPoints; j++) {
			if (i == j) {
				continue;
			}

			double dDiff = 1.0 / (dXsample - dX[j]);
			dDifferential1 += dDiff;
			dDifferential2 += dDiff * dDiff;
		}

		dCoeffs[i] *= (dDifferential1 * dDifferential1 - dDifferential2);
	}
}

///////////////////////////////////////////////////////////////////////////////

void PolynomialInterp::DiffDiffDiffLagrangianPolynomialCoeffs(
	int nPoints,
	const double * dX,
	double * dCoeffs,
	double dXsample
) {
	int i;
	int j;

	int iMatch = (-1);

	// Check if dXsample is equivalent to one of the values of dX
	for (i = 0; i < nPoints; i++) {
		if (fabs(dXsample - dX[i]) < 1.0e-14) {
			iMatch = i;
			break;
		}
	}

	// Match detected
	if (iMatch != (-1)) {

		// Calculate product
		for (i = 0; i < nPoints; i++) {
			dCoeffs[i] = 1.0;

			double dSum = 0.0;
			double dSum2 = 0.0;
			double dSum3 = 0.0;

			for (j = 0; j < nPoints; j++) {
				if ((j == i) || (j == iMatch)) {
					continue;
				}

				dCoeffs[i] *= (dXsample - dX[j]) / (dX[i] - dX[j]);

				double dDiff = 1.0 / (dXsample - dX[j]);

				dSum  += dDiff;
				dSum2 += dDiff * dDiff;
				dSum3 += 2.0 * dDiff * dDiff * dDiff;
			}

			// Multiply by differential
			if (i != iMatch) {
				dCoeffs[i] *= (3.0 * (dSum * dSum - dSum2)
					+ (dXsample - dX[iMatch])
						* (dSum * dSum * dSum - 3.0 * dSum * dSum2 + dSum3))
					/ (dX[i] - dX[iMatch]);
			} else {
				dCoeffs[i] *= (dSum * dSum * dSum - 3.0 * dSum * dSum2 + dSum3);
			}
		}

		return;
	}

	// Compute Lagrangian polynomial coeffs
	LagrangianPolynomialCoeffs(nPoints, dX, dCoeffs, dXsample);

	// Multiply by second moment of differential
	for (i = 0; i < nPoints; i++) {

		double dSum  = 0.0;
		double dSum2 = 0.0;
		double dSum3 = 0.0;
		for (j = 0; j < nPoints; j++) {
			if (i == j) {
				continue;
			}

			double dDiff = 1.0 / (dXsample - dX[j]);
			dSum  += dDiff;
			dSum2 += dDiff * dDiff;
			dSum3 += 2.0 * dDiff * dDiff * dDiff;
		}

		dCoeffs[i] *= (dSum * dSum * dSum - 3.0 * dSum * dSum2 + dSum3);
	}
}

///////////////////////////////////////////////////////////////////////////////

double PolynomialInterp::Interpolate(
	int nPoints,
	const double * dX,
	const double * dY,
	double dXsample
) {
	int i;
	int j;

	// Compute Lagrangian polynomial
	double dValue = 0.0;

	for (i = 0; i < nPoints; i++) {
		double dCoeff = 1.0;

		for (j = 0; j < nPoints; j++) {
			if (i == j) {
				continue;
			}

			dCoeff *= (dXsample - dX[j]) / (dX[i] - dX[j]);
		}

		dValue += dCoeff * dY[i];
	}

	return dValue;
}

///////////////////////////////////////////////////////////////////////////////

void PolynomialInterp::InterpolateCoeffs(
	int nPoints,
	const double * dX,
	const double * dY,
	double * dA,
	double dXmid,
	double * dWorkspace,
	int * iPivot
) {
	// Check if an external workspace is specified
	bool fExternalWorkspace = (dWorkspace == NULL)?(false):(true);
	if (!fExternalWorkspace) {
		dWorkspace = new double[nPoints*nPoints];
	}

	// Check if an external pivot vector is specified
	bool fExternalPivot = (iPivot == NULL)?(false):(true);
	if (!fExternalPivot) {
		iPivot = new int[nPoints];
	}

	// Construct the Vandermonde matrix
	int i;
	int j;

	int k = 0;
	for (j = 0; j < nPoints; j++) {
		dWorkspace[k] = 1.0;
		k++;
	}

	for (i = 1; i < nPoints; i++) {
	for (j = 0; j < nPoints; j++) {
		dWorkspace[k] = (dX[j] - dXmid) * dWorkspace[k-nPoints];

		k++;
	}
	}

	// Initialize A
	memcpy(dA, dY, nPoints * sizeof(double));

	// Store CLAPACK parameters
	int n     = nPoints;
	int nRHS  = 1;
	int nLDA  = nPoints;
	int nLDB  = nPoints;

	int nInfo;

/*
#ifdef USEACML
	// Call the matrix solve
	dgesv(n, nRHS, dWorkspace, nLDA, iPivot, dA, nLDB, &nInfo);
#endif
#ifdef USEESSL
	// Call the matrix solve
	dgesv(n, nRHS, dWorkspace, nLDA, iPivot, dA, nLDB, nInfo);
#endif
#if defined USEVECLIB || defined USEMKL
	// Call the matrix solve
	dgesv_(&n, &nRHS, dWorkspace, &nLDA, iPivot, dA, &nLDB, &nInfo);
#endif
*/
	_EXCEPTION();

	// Delete workspace
	if (!fExternalWorkspace) {
		delete[] dWorkspace;
	}
	if (!fExternalPivot) {
		delete[] iPivot;
	}
}

///////////////////////////////////////////////////////////////////////////////

