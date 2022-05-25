///////////////////////////////////////////////////////////////////////////////
///
///	\file    ClosedContourOp.h
///	\author  Paul Ullrich
///	\version February 6, 2020
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

#ifndef _CLOSEDCONTOUROP_H_
#define _CLOSEDCONTOUROP_H_

#include "Variable.h"

#include <string>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class describing a general closed contour operation.
///	</summary>
class ClosedContourOp {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	ClosedContourOp() :
		m_varix(InvalidVariableIndex),
		m_dDeltaAmount(0.0),
		m_dDistance(0.0),
		m_dMinMaxDist(0.0)
	{ }

public:
	///	<summary>
	///		Parse a closed contour operation string.
	///	</summary>
	void Parse(
		VariableRegistry & varreg,
		const std::string & strOp
	) {
		// Read mode
		enum {
			ReadMode_Amount,
			ReadMode_Distance,
			ReadMode_MinMaxDist,
			ReadMode_Invalid
		} eReadMode = ReadMode_Amount;

		// Get variable information
		int iLast = varreg.FindOrRegisterSubStr(strOp, &m_varix) + 1;

		// Loop through string
		for (int i = iLast; i <= strOp.length(); i++) {

			// Comma-delineated
			if ((i == strOp.length()) || (strOp[i] == ',')) {

				std::string strSubStr =
					strOp.substr(iLast, i - iLast);

				// Read in amount
				if (eReadMode == ReadMode_Amount) {
					m_dDeltaAmount = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Distance;

				// Read in distance
				} else if (eReadMode == ReadMode_Distance) {
					m_dDistance = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_MinMaxDist;

				// Read in min/max distance
				} else if (eReadMode == ReadMode_MinMaxDist) {
					m_dMinMaxDist = atof(strSubStr.c_str());

					iLast = i + 1;
					eReadMode = ReadMode_Invalid;

				// Invalid
				} else if (eReadMode == ReadMode_Invalid) {
					_EXCEPTION1("\nToo many entries in closed contour op \"%s\""
						"\nRequired: \"<name>,<amount>,<distance>,"
						"<minmaxdist>\"",
						strOp.c_str());
				}
			}
		}

		if (eReadMode != ReadMode_Invalid) {
			_EXCEPTION1("\nInsufficient entries in closed contour op \"%s\""
					"\nRequired: \"<name>,<amount>,<distance>,"
					"<minmaxdist>\"",
					strOp.c_str());
		}

		// Output announcement
		if (m_dDeltaAmount == 0.0) {
			_EXCEPTIONT("For closed contour op, delta amount must be non-zero");
		}
		if (m_dDistance <= 0.0) {
			_EXCEPTIONT("For closed contour op, distance must be positive");
		}
		if (m_dMinMaxDist < 0.0) {
			_EXCEPTIONT("For closed contour op, min/max dist must be nonnegative");
		}

		if (m_dDeltaAmount < 0.0) {
			Announce("%s decreases by %f over %f degrees"
				   " (max search %f deg)",
				varreg.GetVariableString(m_varix).c_str(),
				-m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);

		} else {
			Announce("%s increases by %f over %f degrees"
					" (min search %f deg)",
				varreg.GetVariableString(m_varix).c_str(),
				m_dDeltaAmount,
				m_dDistance,
				m_dMinMaxDist);
		}
	}

public:
	///	<summary>
	///		Variable to use for closed contour op.
	///	</summary>
	VariableIndex m_varix;

	///	<summary>
	///		Threshold amount.  If positive this represents a minimum
	///		increase.  If negative this represents a minimum decrease.
	///	</summary>
	double m_dDeltaAmount;

	///	<summary>
	///		Threshold distance (in degrees).
	///	</summary>
	double m_dDistance;

	///	<summary>
	///		Distance to search for min or max (in degrees).
	///	</summary>
	double m_dMinMaxDist;
};

///////////////////////////////////////////////////////////////////////////////

#endif // _CLOSEDCONTOUROP_H_

