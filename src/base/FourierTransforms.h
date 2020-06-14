///////////////////////////////////////////////////////////////////////////////
///
///	\file    FourierTransforms.h
///	\author  Paul Ullrich
///	\version June 13, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _FOURIERTRANSFORMS_H_
#define _FOURIERTRANSFORMS_H_

#include "DataArray1D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Apply a Fourier filter to a given data sequence.
///	</summary>

// TODO: Test this function more thoroughly
template <typename T>
void fourier_filter(
	T * const data,
	size_t sCount,
	size_t sStride,
	size_t sModes,
	DataArray1D<T> & an,
	DataArray1D<T> & bn
) {
	_ASSERT(an.GetRows() >= sModes);
	_ASSERT(bn.GetRows() >= sModes);

	an.Zero();
	bn.Zero();

	{
		for (int n = 0; n < sCount; n++) {
			an[0] += data[n*sStride];
		}

		double dProd0 = 2.0 * M_PI / static_cast<double>(sCount);
		for (int k = 1; k < sModes; k++) {
			double dProd1 = dProd0 * static_cast<double>(k);
			for (int n = 0; n < sCount; n++) {
				double dProd2 = dProd1 * static_cast<double>(n);
				an[k] += data[n*sStride] * cos(dProd2);
				bn[k] -= data[n*sStride] * sin(dProd2);
			}

			//an[sCount-k] = an[k];
			//bn[sCount-k] = -bn[k];
		}
	}
	{
		for (int n = 0; n < sCount; n++) {
			data[n*sStride] = 0.0;
		}

		double dProd0 = 2.0 * M_PI / static_cast<double>(sCount);
		for (int n = 0; n < sCount; n++) {
			data[n*sStride] += an[0];

			double dProd1 = dProd0 * static_cast<double>(n);
			for (int k = 1; k < sModes; k++) {
				double dProd2 = dProd1 * static_cast<double>(k);
				double dProd3 = dProd1 * static_cast<double>(sCount-k);
				data[n*sStride] += an[k] * cos(dProd2) - bn[k] * sin(dProd2);
				data[n*sStride] += an[k] * cos(dProd3) + bn[k] * sin(dProd3);

				//printf("%1.15e %1.15e : ", cos(dProd2), cos(dProd3));
				//printf("%1.15e %1.15e\n", sin(dProd2), sin(dProd3));
			}
			//for (int k = sCount-sModes+1; k < sCount; k++) {
			//	double dProd2 = dProd1 * static_cast<double>(k);
			//	data[n*sStride] += an[k] * cos(dProd2) - bn[k] * sin(dProd2);
			//}
		}

		for (int n = 0; n < sCount; n++) {
			data[n*sStride] /= static_cast<double>(sCount);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

#endif // _FOURIERTRANSFORMS_H_

