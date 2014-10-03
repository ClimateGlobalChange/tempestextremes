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

#include "Announce.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <cstring>
#include <cstdarg>

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

void AnnounceStartBlock(const char * szText) {

	// Do not start a block at maximum indentation level
	if (s_nIndentationLevel == MaximumIndentationLevel) {
		return;
	}

#ifdef USE_MPI
	// Retrieve the rank of this processor
	int nRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank > 0) {
		return;
	}
#endif

	// Check the block flag
	if (s_fBlockFlag) {
		printf("\n");
	}

	// Add indentation
	int i;
	for (i = 0; i < s_nIndentationLevel; i++) {
		printf("..");
	}

	// Output the text
	if (szText != NULL) {
		printf("%s", szText);
		s_fBlockFlag = true;
	}
	s_nIndentationLevel++;

	fflush(NULL);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceEndBlock(const char * szText) {
	// Do not remove a block at minimum indentation level
	if (s_nIndentationLevel == 0) {
		return;
	}

#ifdef USE_MPI
	// Retrieve the rank of this processor
	int nRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank > 0) {
		return;
	}
#endif

	// Check block flag
	if (szText != NULL) {
		if (s_fBlockFlag) {
			s_fBlockFlag = false;

			printf(".. ");
			printf("%s", szText);
			printf("\n");

		} else {
			Announce(szText);
		}
	}

	s_nIndentationLevel--;

	fflush(NULL);
}

///////////////////////////////////////////////////////////////////////////////

void Announce(const char * szText, ...) {

#ifdef USE_MPI
	// Retrieve the rank of this processor
	int nRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank > 0) {
		return;
	}
#endif

	// Turn off the block flag
	if (s_fBlockFlag) {
		printf("\n");
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
	vsprintf(szBuffer, szText, arguments);

	// Cleans up the argument list
	va_end(arguments);

	// Output with proper indentation
	int i;
	for (i = 0; i < s_nIndentationLevel; i++) {
		printf("..");
	}
	printf("%s", szBuffer);
	printf("\n");

	fflush(NULL);
}

///////////////////////////////////////////////////////////////////////////////

void AnnounceBanner(const char * szText) {

#ifdef USE_MPI
	// Retrieve the rank of this processor
	int nRank;

	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if (nRank > 0) {
		return;
	}
#endif

	// Turn off the block flag
	if (s_fBlockFlag) {
		printf("\n");
		s_fBlockFlag = false;
	}

	// No text in banner
	int i;
	if (szText == NULL) {
		for (i = 0; i < BannerSize; i++) {
			printf("-");
		}
		printf("\n");
		fflush(NULL);
		return;
	}

	// Text in banner
	int nLen = strlen(szText) + 2;
	printf("--");
	if (nLen > BannerSize - 2) {
		printf("%s", szText);
		printf("--");
	} else {
		printf(" %s ", szText);
		for (i = 0; i < BannerSize - nLen - 2; i++) {
			printf("-");
		}
	}
	printf("\n");
	fflush(NULL);
}

///////////////////////////////////////////////////////////////////////////////

