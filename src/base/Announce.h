///////////////////////////////////////////////////////////////////////////////
///
///	\file    Announce.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		Functions for making announcements to standard output safely when
///		using MPI.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _ANNOUNCE_H_
#define _ANNOUNCE_H_

#include <cstdio>

///////////////////////////////////////////////////////////////////////////////

extern int g_iVerbosityLevel;

extern FILE * g_fpOutputBuffer;

extern bool g_fOnlyOutputOnRankZero;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the output buffer.
///	</summary>
FILE * AnnounceGetOutputBuffer();

///	<summary>
///		Set the output buffer.
///	</summary>
void AnnounceSetOutputBuffer(FILE * fpOutputBuffer);

///	<summary>
///		Set the verbosity level.
///	</summary>
void AnnounceSetVerbosityLevel(int iVerbosityLevel);

///	<summary>
///		Only output on rank zero.
///	</summary>
void AnnounceOnlyOutputOnRankZero();

///	<summary>
///		Allow output on all ranks.
///	</summary>
void AnnounceOutputOnAllRanks();

///	<summary>
///		Begin a new announcement block.
///	</summary>
void AnnounceStartBlock(const char * szText);

///	<summary>
///		Begin a new announcement block.
///	</summary>
void AnnounceStartBlock(int iVerbosity, const char * szText);

///	<summary>
///		End an announcement block.
///	</summary>
void AnnounceEndBlock(const char * szText);

///	<summary>
///		End an announcement block.
///	</summary>
void AnnounceEndBlock(int iVerbosity, const char * szText);

///	<summary>
///		Make an announcement.
///	</summary>
void Announce(const char * szText, ...);

///	<summary>
///		Make an announcement.
///	</summary>
void Announce(int iVerbosity, const char * szText, ...);

///	<summary>
///		Create a banner / separator containing the specified text.
///	</summary>
void AnnounceBanner(const char * szText = NULL);

///////////////////////////////////////////////////////////////////////////////

#endif

