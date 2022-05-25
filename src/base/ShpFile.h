///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShpFile.h
///	\author  Paul Ullrich
///	\version May 18, 2022
///
///	<remarks>
///		Copyright 2022 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _SHPFILE_H_
#define _SHPFILE_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "Announce.h"
#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

void ReadShpFileAsMesh(
	const std::string & strInputFile,
	Mesh & mesh,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

#endif // _SHPFILE_H_

