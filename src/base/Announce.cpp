///////////////////////////////////////////////////////////////////////////////
///
///	\file    Announce.cpp
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

#ifdef TEMPEST_MPIOMP
#include <mpi.h>
#endif

#include "Announce.h"

#include <cstdio>
#include <cstring>
#include <cstdarg>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Verbosity level.
///	</summary>
int g_iVerbosityLevel = 0;

///	<summary>
///		Output buffer.
///	</summary>
FILE * g_fpAnnounceOutput = stdout;

///	<summary>
///		Only output on rank 0.
///	</summary>
bool g_fOnlyOutputOnRankZero = false;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Maximum announcement buffer size.
///	</summary>
static const int AnnouncementBufferSize = 1024;

///	<summary>
///		Maximum indentation level.
///	</summary>
static const int MaximumIndentationLevel = 16;

///	<summary>
///		Banner size.
///	</summary>
static const int BannerSize = 60;

///	<summary>
///		Current indentation level.
///	</summary>
static int s_nIndentationLevel = 0;

///	<summary>
///		Flag indicating whether a start block is still dangling.
///	</summary>
static bool s_fBlockFlag = false;

///////////////////////////////////////////////////////////////////////////////

FILE * AnnounceGetOutputBuffer() {
	return g_fpAnnounceOutput;
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceSetOutputBuffer(FILE * fpAnnounceOutput) {
	g_fpAnnounceOutput = fpAnnounceOutput;
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceSetVerbosityLevel(int iVerbosityLevel) {
	g_iVerbosityLevel = iVerbosityLevel;
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceOnlyOutputOnRankZero() {
	g_fOnlyOutputOnRankZero = true;
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceOutputOnAllRanks() {
	g_fOnlyOutputOnRankZero = false;
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceStartBlock(
	const char * szText,
	...
) {

	// Do not start a block at maximum indentation level
	if (s_nIndentationLevel == MaximumIndentationLevel) {
		return;
	}
	if (szText == NULL) {
		return;
	}

#ifdef TEMPEST_MPIOMP
	// Only output on rank zero
	if (g_fOnlyOutputOnRankZero) {
		int nRank;

		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank > 0) {
			return;
		}
	}
#endif

	// Check the block flag
	if (s_fBlockFlag) {
		fprintf(g_fpAnnounceOutput, "\n");
	}

	// Build output string from variable argument list
	char szBuffer[AnnouncementBufferSize];
	va_list arguments;
	va_start(arguments, szText);
	int nc = vsnprintf(szBuffer, AnnouncementBufferSize, szText, arguments);
	if (nc > AnnouncementBufferSize-2) {
		szBuffer[AnnouncementBufferSize-4] = '.';
		szBuffer[AnnouncementBufferSize-3] = '.';
		szBuffer[AnnouncementBufferSize-2] = '.';
		szBuffer[AnnouncementBufferSize-1] = '\0';
	}
	va_end(arguments);

	// Output with proper indentation
	for (int i = 0; i < s_nIndentationLevel; i++) {
		fprintf(g_fpAnnounceOutput, "..");
	}
	fprintf(g_fpAnnounceOutput, "%s", szBuffer);

	s_fBlockFlag = true;
	s_nIndentationLevel++;

	fflush(g_fpAnnounceOutput);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceStartBlock(
	int iVerbosity,
	const char * szText
) {
	if (iVerbosity > g_iVerbosityLevel) {
		return;
	}

	AnnounceStartBlock(szText);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceEndBlock(
	const char * szText,
	...
) {
	// Do not remove a block at minimum indentation level
	if (s_nIndentationLevel == 0) {
		return;
	}

#ifdef TEMPEST_MPIOMP
	// Only output on rank zero
	if (g_fOnlyOutputOnRankZero) {
		int nRank;

		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank > 0) {
			return;
		}
	}
#endif

	// Check block flag
	if (szText != NULL) {

		// Build output string from variable argument list
		char szBuffer[AnnouncementBufferSize];

		va_list arguments;
		va_start(arguments, szText);
		int nc = vsnprintf(szBuffer, AnnouncementBufferSize, szText, arguments);
		if (nc > AnnouncementBufferSize-2) {
			szBuffer[AnnouncementBufferSize-4] = '.';
			szBuffer[AnnouncementBufferSize-3] = '.';
			szBuffer[AnnouncementBufferSize-2] = '.';
			szBuffer[AnnouncementBufferSize-1] = '\0';
		}

		va_end(arguments);

		if (s_fBlockFlag) {
			s_fBlockFlag = false;

			fprintf(g_fpAnnounceOutput, ".. ");
			fprintf(g_fpAnnounceOutput, "%s", szBuffer);
			fprintf(g_fpAnnounceOutput, "\n");

		} else {
			Announce(szBuffer);
		}

	} else {
		if (s_fBlockFlag) {
			s_fBlockFlag = false;
			fprintf(g_fpAnnounceOutput, "\n");
		}
	}

	s_nIndentationLevel--;

	fflush(g_fpAnnounceOutput);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceEndBlock(
	int iVerbosity,
	const char * szText
) {
	if (iVerbosity > g_iVerbosityLevel) {
		return;
	}

	AnnounceEndBlock(szText);
}

///////////////////////////////////////////////////////////////////////////////

void Announce(const char * szText, ...) {

#ifdef TEMPEST_MPIOMP
	// Only output on rank zero
	if (g_fOnlyOutputOnRankZero) {
		int nRank;

		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank > 0) {
			return;
		}
	}
#endif

	// Turn off the block flag
	if (s_fBlockFlag) {
		fprintf(g_fpAnnounceOutput, "\n");
		s_fBlockFlag = false;
	}

	// If no text, return
	if (szText == NULL) {
		return;
	}

	// Output buffer
	char szBuffer[AnnouncementBufferSize];

	va_list arguments;

	// Initialize the argument list
	va_start(arguments, szText);

	// Write to string
	int nc = vsnprintf(szBuffer, AnnouncementBufferSize, szText, arguments);
	if (nc > AnnouncementBufferSize-2) {
		szBuffer[AnnouncementBufferSize-4] = '.';
		szBuffer[AnnouncementBufferSize-3] = '.';
		szBuffer[AnnouncementBufferSize-2] = '.';
		szBuffer[AnnouncementBufferSize-1] = '\0';
	}

	// Cleans up the argument list
	va_end(arguments);

	// Output with proper indentation
	int i;
	for (i = 0; i < s_nIndentationLevel; i++) {
		fprintf(g_fpAnnounceOutput, "..");
	}
	fprintf(g_fpAnnounceOutput, "%s", szBuffer);
	fprintf(g_fpAnnounceOutput, "\n");

	fflush(g_fpAnnounceOutput);
}

///////////////////////////////////////////////////////////////////////////////

void Announce(
	int iVerbosity,
	const char * szText,
	...
) {

#ifdef TEMPEST_MPIOMP
	// Only output on rank zero
	if (g_fOnlyOutputOnRankZero) {
		int nRank;

		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank > 0) {
			return;
		}
	}
#endif

	// Check verbosity
	if (iVerbosity > g_iVerbosityLevel) {
		return;
	}

	// Turn off the block flag
	if (s_fBlockFlag) {
		fprintf(g_fpAnnounceOutput, "\n");
		s_fBlockFlag = false;
	}

	// If no text, return
	if (szText == NULL) {
		return;
	}

	// Output buffer
	char szBuffer[AnnouncementBufferSize];

	va_list arguments;

	// Initialize the argument list
	va_start(arguments, szText);

	// Write to string
	int nc = vsnprintf(szBuffer, AnnouncementBufferSize, szText, arguments);
	if (nc > AnnouncementBufferSize-2) {
		szBuffer[AnnouncementBufferSize-4] = '.';
		szBuffer[AnnouncementBufferSize-3] = '.';
		szBuffer[AnnouncementBufferSize-2] = '.';
		szBuffer[AnnouncementBufferSize-1] = '\0';
	}

	// Cleans up the argument list
	va_end(arguments);

	// Output with proper indentation
	int i;
	for (i = 0; i < s_nIndentationLevel; i++) {
		fprintf(g_fpAnnounceOutput, "..");
	}
	fprintf(g_fpAnnounceOutput, "%s", szBuffer);
	fprintf(g_fpAnnounceOutput, "\n");

	fflush(g_fpAnnounceOutput);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceBanner(const char * szText) {

#ifdef TEMPEST_MPIOMP
	// Only output on rank zero
	if (g_fOnlyOutputOnRankZero) {
		int nRank;

		MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

		if (nRank > 0) {
			return;
		}
	}
#endif

	// Turn off the block flag
	if (s_fBlockFlag) {
		fprintf(g_fpAnnounceOutput, "\n");
		s_fBlockFlag = false;
	}

	// No text in banner
	int i;
	if (szText == NULL) {
		for (i = 0; i < BannerSize; i++) {
			printf("-");
		}
		fprintf(g_fpAnnounceOutput, "\n");
		fflush(g_fpAnnounceOutput);
		return;
	}

	// Text in banner
	int nLen = strlen(szText) + 2;
	printf("--");
	if (nLen > BannerSize - 2) {
		fprintf(g_fpAnnounceOutput, "%s", szText);
		fprintf(g_fpAnnounceOutput, "--");
	} else {
		fprintf(g_fpAnnounceOutput, " %s ", szText);
		for (i = 0; i < BannerSize - nLen - 2; i++) {
			fprintf(g_fpAnnounceOutput, "-");
		}
	}
	fprintf(g_fpAnnounceOutput, "\n");
	fflush(g_fpAnnounceOutput);
}

///////////////////////////////////////////////////////////////////////////////

