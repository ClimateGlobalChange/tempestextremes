///////////////////////////////////////////////////////////////////////////////
///
///	\file    LegendrePolynomial.cpp
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

#include "LegendrePolynomial.h"
#include "Exception.h"

#include <cmath>
#include <cstring>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////

void LegendrePolynomial::EvaluateValueAndDerivative(
	int nDegree,
	double dX,
	double & dValue,
	double & dDerivative
) {

	// Low order cases
	if (nDegree < 9) {
		if (nDegree == 0) {
			dValue = 1.0;
			dDerivative = 0.0;

		} else if (nDegree == 1) {
			dValue = dX;
			dDerivative = 1.0;

		} else if (nDegree == 2) {
			dValue = (1.5 * dX * dX - 0.5);
			dDerivative = 3.0 * dX;

		} else if (nDegree == 3) {
			dValue = ((2.5 * dX * dX - 1.5) * dX);
			dDerivative = (7.5 * dX * dX - 1.5);

		} else if (nDegree == 4) {
			dValue = ((4.375 * dX * dX - 3.75) * dX * dX + 0.375);
			dDerivative = (2.5 * (7.0 * dX * dX - 3.0) * dX);

		} else if (nDegree == 5) {
			dValue = (((7.875 * dX * dX - 8.75) * dX * dX + 1.875) * dX);
			dDerivative = ((39.375 * dX * dX - 26.25) * dX * dX + 1.875);

		} else if (nDegree == 6) {
			dValue = (((14.4375 * dX * dX - 19.6875) * dX * dX + 6.5625) * dX * dX - 0.3125);
			dDerivative = (((86.625 * dX * dX - 78.75) * dX * dX + 13.125) * dX);

		} else if (nDegree == 7) {
			dValue = ((((26.8125 * dX * dX - 43.3125) * dX * dX + 19.6875) * dX * dX - 2.1875) * dX);
			dDerivative = (((187.6875 * dX * dX - 216.5625) * dX * dX + 59.0625) * dX * dX - 2.1875);

		} else if (nDegree == 8) {
			dValue = ((((50.2734375 * dX * dX - 93.84375) * dX * dX + 54.140625) * dX * dX - 9.84375) * dX * dX + 0.2734375);
			dDerivative = ((((402.1875 * dX * dX - 563.0625) * dX * dX + 216.5625) * dX * dX - 19.6875) * dX);
		}

		return;
	}

	// Arbitrary order cases
	// Algorithm adapted from Shanjie Zhang and Jianming Jin (1996)
	// "Computation of Special Functions" ISBN 0-471-11963-6
	{
		double * dPn = new double[nDegree + 1];
		double * dPd = new double[nDegree + 1];

		dPn[0] = 1.0;
		dPn[1] = dX;
		dPd[0] = 0.0;
		dPd[1] = 1.0;

		double dP0 = 1.0;
		double dP1 = dX;

		double dXk = dX * dX;

		for (int k = 2; k <= nDegree; k++) {
			double dK = static_cast<double>(k);
			double dPf = (2.0 * dK - 1.0) / dK * dX * dP1
				- (dK - 1.0) / dK * dP0;

			dXk = dXk * dX;

			dPn[k] = dPf;
			if (fabs(dX) == 1.0) {
				dPd[k] = 0.5 * dXk * dK * (dK + 1.0);
			} else {
				dPd[k] = dK * (dP1 - dX * dPf) / (1.0 - dX * dX);
			}
			dP0 = dP1;
			dP1 = dPf;
		}

		dValue = dPn[nDegree];
		dDerivative = dPd[nDegree];

		delete[] dPn;
		delete[] dPd;
	}
}

///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::Evaluate(
	int nDegree,
	double dX
) {
	double dValue;
	double dDerivative;

	EvaluateValueAndDerivative(nDegree, dX, dValue, dDerivative);

	return dValue;
}

///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::EvaluateDerivative(
	int nDegree,
	double dX
) {
	double dValue;
	double dDerivative;

	EvaluateValueAndDerivative(nDegree, dX, dValue, dDerivative);

	return dDerivative;
}

///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::DerivativeRoot(
	int nDegree,
	int nRoot
) {
	if ((nRoot >= (nDegree - 1)) || (nRoot < 0)) {
		_EXCEPTIONT("Invalid root.");
	}

	// Second degree polynomial
	if (nDegree == 2) {
		return 0.0;
	}

	// Third degree polynomial
	if (nDegree == 3) {
		if (nRoot == 0) {
			return (- 0.44721359549995793928);
		} else {
			return (+ 0.44721359549995793928);
		}
	}

	// Fourth degree polynomial
	if (nDegree == 4) {
		if (nRoot == 0) {
			return (- 0.65465367070797714380);
		} else if (nRoot == 1) {
			return (0.0);
		} else {
			return (+ 0.65465367070797714380);
		}
	}

	// Fifth degree polynomial
	if (nDegree == 5) {
		if (nRoot == 0) {
			return (- 0.76505532392946469285);
		} else if (nRoot == 1) {
			return (- 0.28523151648064509632);
		} else if (nRoot == 2) {
			return (+ 0.28523151648064509632);
		} else {
			return (+ 0.76505532392946469285);
		}
	}

	// Sixth degree polynomial
	if (nDegree == 6) {
		if (nRoot == 0) {
			return (- 0.83022389627856692987);
		} else if (nRoot == 1) {
			return (- 0.46884879347071421380);
		} else if (nRoot == 2) {
			return (0.0);
		} else if (nRoot == 3) {
			return (+ 0.46884879347071421380);
		} else {
			return (+ 0.83022389627856692987);
		}
	}

	// Seventh degree polynomial
	if (nDegree == 7) {
		if (nRoot == 0) {
			return (- 0.87174014850960661534);
		} else if (nRoot == 1) {
			return (- 0.59170018143314230214);
		} else if (nRoot == 2) {
			return (- 0.20929921790247886877);
		} else if (nRoot == 3) {
			return (+ 0.20929921790247886877);
		} else if (nRoot == 4) {
			return (+ 0.59170018143314230214);
		} else {
			return (+ 0.87174014850960661534);
		}
	}

	// Eighth degree polynomial
	if (nDegree == 8) {
		if (nRoot == 0) {
			return (- 0.89975799541146015731);
		} else if (nRoot == 1) {
			return (- 0.67718627951073775344);
		} else if (nRoot == 2) {
			return (- 0.36311746382617815871);
		} else if (nRoot == 3) {
			return (0.0);
		} else if (nRoot == 4) {
			return (+ 0.36311746382617815871);
		} else if (nRoot == 5) {
			return (+ 0.67718627951073775344);
		} else {
			return (+ 0.89975799541146015731);
		}
	}

	_EXCEPTIONT("Unimplemented - maximum degree is 8.");
}

///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::DerivativeExtendedRoot(
	int nDegree,
	int nRoot
) {
	if ((nRoot >= (nDegree + 1)) || (nRoot < 0)) {
		_EXCEPTIONT("Invalid root.");
	}

	// First root is -1
	if (nRoot == 0) {
		return -1.0;

	// Last root is 1
	} else if (nRoot == nDegree) {
		return 1.0;

	// Otherwise handle interior roots
	} else {
		return DerivativeRoot(nDegree, nRoot - 1);
	}
}


///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::EvaluateCharacteristic(
	int nDegree,
	int nRoot,
	double dX
) {
	if (nDegree == 0) {
		_EXCEPTIONT("Polynomial of degree zero has no characteristic.");
	}

	if ((nRoot < 0) || (nRoot >= nDegree)) {
		_EXCEPTIONT("Invalid root.");
	}

	// First degree polynomial
	if (nDegree == 1) {
		return 1.0;
	}

	// Second degree polynomial
	if (nDegree == 2) {
		if (nRoot == 0) {
			return (- 0.86602540378443864677 * dX + 0.5);
		} else {
			return (+ 0.86602540378443864677 * dX + 0.5);
		}
	}

	// Third degree polynomial
	if (nDegree == 3) {
		if (nRoot == 0) {
			return (dX / 6.0 * (5.0 * dX - 3.8729833462074168852));
		} else if (nRoot == 1) {
			return (-5.0/3.0 * dX * dX + 1.0);
		} else {
			return (dX / 6.0 * (5.0 * dX + 3.8729833462074168852));
		}
	}

	// Fourth degree polynomial
	if (nDegree == 4) {
		if (nRoot == 0) {
			return (- 0.00075719796612642307197
				* (35.0 * dX * dX - 4.045548849896677731)
				* (35.0 * dX - 30.139770905791840133));

		} else if (nRoot == 1) {
			return (+ 0.0019179030007709245726
				* (35.0 * dX * dX - 25.954451150103322269)
				* (35.0 * dX - 11.899336525469969268));

		} else if (nRoot == 2) {
			return (- 0.0019179030007709245726
				* (35.0 * dX * dX - 25.954451150103322269)
				* (35.0 * dX + 11.899336525469969268));

		} else {
			return (+ 0.00075719796612642307197
				* (35.0 * dX * dX - 4.045548849896677731)
				* (35.0 * dX + 30.139770905791840133));
		}
	}

	// Fifth degree polynomial
	if (nDegree == 5) {
		if (nRoot == 0) {
			return (+ 0.00086638894153586418251 * dX
				* (63.0 * dX * dX - 18.266799469318489040)
				* (21.0 * dX - 19.029776764711943849));

		} else if (nRoot == 1) {
			return (- 0.0024536905288374514841 * dX
				* (63.0 * dX * dX - 51.733200530681510960)
				* (21.0 * dX - 11.307855512219344912));

		} else if (nRoot == 2) {
			return ((21.0 / 5.0 * dX * dX - 14.0 / 3.0) * dX * dX + 1.0);

		} else if (nRoot == 3) {
			return (- 0.0024536905288374514841 * dX
				* (63.0 * dX * dX - 51.733200530681510960)
				* (21.0 * dX + 11.307855512219344912));

		} else {
			return (+ 0.00086638894153586418251 * dX
				* (63.0 * dX * dX - 18.266799469318489040)
				* (21.0 * dX + 19.029776764711943849));
		}
	}

	// Sixth degree polynomial
	if (nDegree == 6) {
		if (nRoot == 0) {
			return (((((
				(-1.5264865741891266785 * dX)
				+ 1.4234021942717687413) * dX
				+ 0.7542934485586941013) * dX
				- 0.7033556455441457329) * dX
				- 0.0379998366097363458) * dX
				+ 0.0354336891832800019);

		} else if (nRoot == 1) {
			return (((((
				( 4.6000760138804121676 * dX)
				- 3.0416134388360470124) * dX
				- 4.2616875722584162543) * dX
				+ 2.8178678249638917282) * dX
				+ 0.2277429869855125288) * dX
				- 0.1505858006966852032);

		} else if (nRoot == 2) {
			return (((((
				(-6.7815638429010193828 * dX)
				+ 1.6182112445642782709) * dX
				+ 8.8614508084127848892) * dX
				- 2.1145121794197459957) * dX
				- 2.5779658442843167374) * dX
				+ 0.6151521115134052013);

		} else if (nRoot == 3) {
			return (((((
				(+6.7815638429010193828 * dX)
				+ 1.6182112445642782709) * dX
				- 8.8614508084127848892) * dX
				- 2.1145121794197459957) * dX
				+ 2.5779658442843167374) * dX
				+ 0.6151521115134052013);

		} else if (nRoot == 4) {
			return (((((
				(-4.6000760138804121676 * dX)
				- 3.0416134388360470124) * dX
				+ 4.2616875722584162543) * dX
				+ 2.8178678249638917282) * dX
				- 0.2277429869855125288) * dX
				- 0.1505858006966852032);

		} else {
			return (((((
				( 1.526486574189126678 * dX)
				+ 1.4234021942717687413) * dX
				- 0.75429344855869410144) * dX
				- 0.70335564554414573290) * dX
				+ 0.037999836609736345824) * dX
				+ 0.035433689183280001923);
		}
	}

	// Seventh degree polynomial
	if (nDegree == 7) {
		if (nRoot == 0) {
			return ((((((
				( 2.1486964239590934589 * dX)
				- 2.0393447772021659821) * dX
				- 1.5354128823580029005) * dX
				+ 1.4572725153589816235) * dX
				+ 0.1946052826856376549) * dX
				- 0.1847014135806379266) * dX);

		} else if (nRoot == 1) {
			return ((((((
				(-6.7273253820885753268 * dX)
				+ 4.9885215664930404665) * dX
				+ 7.1680736133225813075) * dX
				- 5.3153501249508289727) * dX
				- 0.9981466312931015002) * dX
				+ 0.7401568549048151787) * dX);

		} else if (nRoot == 2) {
			return ((((((
				( 10.707200386700910439 * dX)
				-  4.345465361768756480) * dX
				- 15.532660730964578407) * dX
				+  6.303855045652071854) * dX
				+  5.303541348607463845) * dX
				-  2.152416541461881283) * dX);

		} else if (nRoot == 3) {
			return (((
				(-12.257142857142857143 * dX * dX)
				+ 19.800000000000000000) * dX * dX
				-  9.000000000000000000) * dX * dX
				+  1.000000000000000000);

		} else if (nRoot == 4) {
			return ((((((
				( 10.707200386700910439 * dX)
				+  4.345465361768756480) * dX
				- 15.532660730964578407) * dX
				-  6.303855045652071854) * dX
				+  5.303541348607463845) * dX
				+  2.152416541461881283) * dX);

		} else if (nRoot == 5) {
			return ((((((
				(-6.7273253820885753268 * dX)
				- 4.9885215664930404665) * dX
				+ 7.1680736133225813075) * dX
				+ 5.3153501249508289727) * dX
				- 0.9981466312931015002) * dX
				- 0.7401568549048151787) * dX);

		} else {
			return ((((((
				( 2.1486964239590934589 * dX)
				+ 2.0393447772021659821) * dX
				- 1.5354128823580029005) * dX
				- 1.4572725153589816235) * dX
				+ 0.1946052826856376549) * dX
				+ 0.1847014135806379266) * dX);
		}
	}

	// Eighth degree polynomial
	if (nDegree == 8) {
		if (nRoot == 0) {
			return (((((((
				(-3.1556289885043198449 * dX)
				+ 3.0303185085302787149) * dX
				+ 2.9805233195096943401) * dX
				- 2.8621663107795247093) * dX
				- 0.6498604040771765922) * dX
				+ 0.6240543541747028191) * dX
				+ 0.0186123301334081728) * dX
				- 0.0178732318328953036);

		} else if (nRoot == 1) {
			return (((((((
				( 10.132361625580626384 * dX)
				-  8.072112844132326405) * dX
				- 12.482959996263643601) * dX
				+  9.944755767918576155) * dX
				+  2.989120511951059053) * dX
				-  2.381332108820866783) * dX
				-  0.086831457007825462) * dX
				+  0.069175710983117085);

		} else if (nRoot == 2) {
			return (((((((
				(-16.939455205503316327 * dX)
				+  8.902232716817861695) * dX
				+ 26.941904570300910930) * dX
				- 14.158844036565995668) * dX
				- 10.801558793145292750) * dX
				+  5.676569223414559800) * dX
				+  0.333595299836857488) * dX
				-  0.175315141860024065);

		} else if (nRoot == 3) {
			return (((((((
				( 21.045307084278615457 * dX)
				-  3.860438381215814005) * dX
				- 38.576435089651940956) * dX
				+  7.076254579426944221) * dX
				+ 21.366146631006968828) * dX
				-  3.919291468768395837) * dX
				-  3.401825599679738188) * dX
				+  0.624012662709802284);

		} else if (nRoot == 4) {
			return (((((((
				(-21.045307084278615457 * dX)
				-  3.860438381215814005) * dX
				+ 38.576435089651940956) * dX
				+  7.076254579426944221) * dX
				- 21.366146631006968828) * dX
				-  3.919291468768395837) * dX
				+  3.401825599679738188) * dX
				+  0.624012662709802284);

		} else if (nRoot == 5) {
			return (((((((
				( 16.939455205503316327 * dX)
				+  8.902232716817861695) * dX
				- 26.941904570300910930) * dX
				- 14.158844036565995668) * dX
				+ 10.801558793145292750) * dX
				+  5.676569223414559800) * dX
				-  0.333595299836857488) * dX
				-  0.175315141860024065);

		} else if (nRoot == 6) {
			return (((((((
				(-10.132361625580626384 * dX)
				-  8.072112844132326405) * dX
				+ 12.482959996263643601) * dX
				+  9.944755767918576155) * dX
				-  2.989120511951059053) * dX
				-  2.381332108820866783) * dX
				+  0.086831457007825462) * dX
				+  0.069175710983117085);

		} else {
			return (((((((
				( 3.1556289885043198449 * dX)
				+ 3.0303185085302787149) * dX
				- 2.9805233195096943401) * dX
				- 2.8621663107795247093) * dX
				+ 0.6498604040771765922) * dX
				+ 0.6240543541747028191) * dX
				- 0.0186123301334081728) * dX
				- 0.0178732318328953036);
		}
	}

	_EXCEPTIONT("Unimplemented - maximum degree is 8.");
}

///////////////////////////////////////////////////////////////////////////////

void LegendrePolynomial::AllRoots(
	int nDegree,
	double * dRoots
) {
	// Iteration tolerance
	const double IterationTolerance = 1.0e-14;

	// Check for degree 0
	if (nDegree == 0) {
		return;
	}

	// Check for negative degree
	if (nDegree < 0) {
		_EXCEPTION1("Invalid degree (%i)", nDegree);
	}

	// Check for NULL dRoots pointer
	if (dRoots == NULL) {
		_EXCEPTIONT("NULL pointer passed into AllRoots argument dRoots");
	}

	// Double degree
	double dDegree = static_cast<double>(nDegree);

	// Additional storage for roots
	double * dBuffer = new double[nDegree];

	// Set initial estimate
	for (int k = 0; k < nDegree; k++) {
		dRoots[k] = (static_cast<double>(k) + 0.5) * 2.0 / dDegree  - 1.0;
	}

	// Refine estimate using Aberth-Ehrlich method
	for (int iter = 0; iter < nDegree + 10; iter++) {

		int nFoundRoots = 0;

		for (int k = 0; k < nDegree; k++) {
			double dValue;
			double dDerivative;

			EvaluateValueAndDerivative(nDegree, dRoots[k], dValue, dDerivative);

			if (fabs(dValue) < IterationTolerance) {
				nFoundRoots++;
				dBuffer[k] = dRoots[k];
				continue;
			}

			double dMod = 0.0;
			for (int l = 0; l < nDegree; l++) {
				if (l == k) {
					continue;
				}

				dMod += 1.0 / (dRoots[k] - dRoots[l]);
			}

			dBuffer[k] = dRoots[k] - 1.0 / (dDerivative / dValue - dMod);
		}

		// Copy buffer to roots
		memcpy(dRoots, dBuffer, nDegree * sizeof(double));

		// Check if all roots have been found
		if (nFoundRoots == nDegree) {
			break;
		}
	}

	// Sort roots
	std::sort(dRoots, dRoots + nDegree);

	// Delete buffer
	delete[] dBuffer;
}

///////////////////////////////////////////////////////////////////////////////

void LegendrePolynomial::AllDerivativeRoots(
	int nDegree,
	double * dRoots
) {
	// Iteration tolerance
	const double IterationTolerance = 1.0e-14;

	// Check for degree 0 or 1
	if ((nDegree == 0) || (nDegree == 1)) {
		return;
	}

	// Check for negative degree
	if (nDegree < 0) {
		_EXCEPTION1("Invalid degree (%i)", nDegree);
	}

	// Check for NULL dRoots pointer
	if (dRoots == NULL) {
		_EXCEPTIONT("NULL pointer passed into AllRoots argument dRoots");
	}

	// Double degree
	double dDegree = static_cast<double>(nDegree);

	// Additional storage for roots
	double * dBuffer = new double[nDegree-1];

	// Set initial estimate
	for (int k = 0; k < nDegree-1; k++) {
		dRoots[k] = (static_cast<double>(k) + 0.5) * 2.0 / (dDegree - 1.0) - 1.0;
	}

	// Refine estimate using Aberth-Ehrlich method
	for (int iter = 0; iter < nDegree + 10; iter++) {
/*
		for (int k = 0; k < nDegree - 1; k++) {
			printf("%i %i: %1.15e\n", iter, k, dRoots[k]);
		}
*/
		int nFoundRoots = 0;

		for (int k = 0; k < nDegree - 1; k++) {
			double dValue;
			double dDerivative;

			EvaluateValueAndDerivative(nDegree, dRoots[k], dValue, dDerivative);

			if (fabs(dDerivative) < IterationTolerance) {
				nFoundRoots++;
				dBuffer[k] = dRoots[k];
				continue;
			}

			double dMod = 0.0;
			for (int l = 0; l < nDegree - 1; l++) {
				if (l == k) {
					continue;
				}

				dMod += 1.0 / (dRoots[k] - dRoots[l]);
			}

			// Second derivative obtained from 
			double dSecondDerivative =
				1.0 / (1.0 - dRoots[k] * dRoots[k])
				* (2.0 * dRoots[k] * dDerivative
					- dDegree * (dDegree + 1.0) * dValue);

			dBuffer[k] = dRoots[k]
				- 1.0 / (dSecondDerivative / dDerivative - dMod);
		}

		// Copy buffer to roots
		memcpy(dRoots, dBuffer, (nDegree - 1) * sizeof(double));

		// Check if all roots have been found
		if (nFoundRoots == nDegree - 1) {
			break;
		}
	}

	// Sort roots
	std::sort(dRoots, dRoots + nDegree - 1);

	// Delete buffer
	delete[] dBuffer;
}

///////////////////////////////////////////////////////////////////////////////

double LegendrePolynomial::Root(
	int nDegree,
	int nRoot
) {
	if ((nRoot >= nDegree) || (nRoot < 0)) {
		_EXCEPTIONT("Invalid root.");
	}

	// First degree polynomial
	if (nDegree == 1) {
		return 0.0;
	}

	// Second degree polynomial
	if (nDegree == 2) {
		if (nRoot == 0) {
			return (- 0.57735026918962576451);
		} else {
			return (+ 0.57735026918962576451);
		}
	}

	// Third degree polynomial
	if (nDegree == 3) {
		if (nRoot == 0) {
			return (- 0.77459666924148337704);
		} else if (nRoot == 1) {
			return (0.0);
		} else {
			return (+ 0.77459666924148337704);
		}
	}

	// Fourth degree polynomial
	if (nDegree == 4) {
		if (nRoot == 0) {
			return (- 0.86113631159405257524);
		} else if (nRoot == 1) {
			return (- 0.33998104358485626481);
		} else if (nRoot == 2) {
			return (+ 0.33998104358485626481);
		} else {
			return (+ 0.86113631159405257524);
		}
	}

	// Fifth degree polynomial
	if (nDegree == 5) {
		if (nRoot == 0) {
			return (- 0.90617984593866399282);
		} else if (nRoot == 1) {
			return (- 0.53846931010568309105);
		} else if (nRoot == 2) {
			return (0.0);
		} else if (nRoot == 3) {
			return (+ 0.53846931010568309105);
		} else {
			return (+ 0.90617984593866399282);
		}
	}

	// Sixth degree polynomial
	if (nDegree == 6) {
		if (nRoot == 0) {
			return (- 0.93246951420315202781);
		} else if (nRoot == 1) {
			return (- 0.66120938646626451366);
		} else if (nRoot == 2) {
			return (- 0.23861918608319690863);
		} else if (nRoot == 3) {
			return (+ 0.23861918608319690863);
		} else if (nRoot == 4) {
			return (+ 0.66120938646626451366);
		} else {
			return (+ 0.93246951420315202781);
		}
	}

	// Seventh degree polynomial
	if (nDegree == 7) {
		if (nRoot == 0) {
			return (- 0.94910791234275852453);
		} else if (nRoot == 1) {
			return (- 0.74153118559939443986);
		} else if (nRoot == 2) {
			return (- 0.40584515137739716691);
		} else if (nRoot == 3) {
			return (0.0);
		} else if (nRoot == 4) {
			return (+ 0.40584515137739716691);
		} else if (nRoot == 5) {
			return (+ 0.74153118559939443986);
		} else {
			return (+ 0.94910791234275852453);
		}
	}

	// Eighth degree polynomial
	if (nDegree == 8) {
		if (nRoot == 0) {
			return (- 0.96028985649753623168);
		} else if (nRoot == 1) {
			return (- 0.79666647741362673961);
		} else if (nRoot == 2) {
			return (- 0.52553240991632898582);
		} else if (nRoot == 3) {
			return (- 0.18343464249564980494);
		} else if (nRoot == 4) {
			return (+ 0.18343464249564980494);
		} else if (nRoot == 5) {
			return (+ 0.52553240991632898582);
		} else if (nRoot == 6) {
			return (+ 0.79666647741362673961);
		} else {
			return (+ 0.96028985649753623168);
		}
	}

	_EXCEPTIONT("Unimplemented - maximum degree is 8.");
}

///////////////////////////////////////////////////////////////////////////////

