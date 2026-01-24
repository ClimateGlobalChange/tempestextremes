#ifndef STITCHBLOBS_MPI_DEBUG_H
#define STITCHBLOBS_MPI_DEBUG_H

///////////////////////////////////////////////////////////////////////////////
///
/// \file    StitchBlobsMPIDebug.h
/// \author  Hongyu Chen
/// \version January 21st, 2026
///
/// \remarks
///   Copyright 2000-2026 Paul Ullrich, Hongyu Chen
///
///   This file is distributed as part of the Tempest source code package.
///   Permission is granted to use, copy, modify and distribute this
///   source code and its documentation under the terms of the GNU General
///   Public License. This software is provided "as is" without express
///   or implied warranty.
///
///////////////////////////////////////////////////////////////////////////////


// Enable with -DSTITCHBLOBS_MPI_DEBUG=1
#if !defined(STITCHBLOBS_MPI_DEBUG)
#define STITCHBLOBS_MPI_DEBUG 0
#endif

#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <fstream>
#include <sstream>
#include <limits>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif


#if STITCHBLOBS_MPI_DEBUG


// Current RSS in bytes (Linux)
static inline size_t GetRSSBytes() {
    std::ifstream s("/proc/self/status");
    std::string key;
    size_t kb = 0;
    while (s >> key) {
        if (key == "VmRSS:") {
            if (s >> kb) {
                std::string unit; s >> unit; // expect "kB"
                return kb * 1024ULL;
            }
            return 0ULL;
        }
        s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0ULL;
}

static inline double GB(size_t bytes) {
    return double(bytes) / (1024.0 * 1024.0 * 1024.0);
}

template <class T>
inline double VecCapGiB(const std::vector<T>& v) {
    return double(v.capacity()) * sizeof(T)
           / (1024.0 * 1024.0 * 1024.0);
}

// printf-style formatting helper
inline std::string vformat(const char* fmt, va_list ap) {
    char buf[4096];
    vsnprintf(buf, sizeof(buf), fmt ? fmt : "", ap);
    return std::string(buf);
}

// Safe bytes calculators
template <class T>
inline double BytesToGiB(size_t n) {
    return double(n) * sizeof(T)
           / (1024.0 * 1024.0 * 1024.0);
}

template <class T>
inline double VecBytesGiB(const std::vector<T>& v) {
    return BytesToGiB<T>(v.size());
}

// Bounds-check helper for slice walkers
inline void AssertSliceBounds(
    const char* where,
    size_t off,
    size_t len,
    size_t total
) {
    if (off > total || len > total || off + len > total) {
        std::fprintf(
            stderr,
            "FATAL %s: slice OOB (off=%zu len=%zu total=%zu)\n",
            where, off, len, total
        );
        std::fflush(stderr);
        std::abort(); // debug-only abort
    }
}

static inline double EstimateSetOverheadGiB(size_t nPoints) {
    // conservative average per std::set node
    const double per = 40.0; // bytes
    return (nPoints * per) / (1024.0 * 1024.0 * 1024.0);
}


// ========================= MPI debug helpers =========================

#if defined(TEMPEST_MPIOMP)

inline void MpiDbgAtF(
    MPI_Comm comm,
    const char* where,
    const char* fmt,
    ...
) {
    int r = -1;
    MPI_Comm_rank(comm, &r);

    double rssGB = GB(GetRSSBytes());

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(1);

    oss << "[rank " << r << "] "
        << (where ? where : "")
        << " | RSS=" << rssGB << " GB";

    if (fmt) {
        va_list ap;
        va_start(ap, fmt);
        oss << " : " << vformat(fmt, ap);
        va_end(ap);
    }

    oss << "\n";
    std::cout << oss.str() << std::flush;
}

inline void MpiDbgF(const char* fmt, ...) {
    int r = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    double rssGB = GB(GetRSSBytes());

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(1);

    oss << "[rank " << r << "] RSS=" << rssGB << " GB : ";

    va_list ap;
    va_start(ap, fmt);
    oss << vformat(fmt, ap);
    va_end(ap);

    oss << "\n";
    std::cout << oss.str() << std::flush;
}

inline void MpiDbgAtFOrdered(
    MPI_Comm comm,
    const char* where,
    const char* fmt,
    ...
) {
    int me = -1, n = -1;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &n);

    va_list ap;
    va_start(ap, fmt);
    std::string msg = vformat(fmt, ap);
    va_end(ap);

    for (int r = 0; r < n; ++r) {
        MPI_Barrier(comm);
        if (r == me) {
            double rssGB = GB(GetRSSBytes());
            std::ostringstream oss;
            oss.setf(std::ios::fixed);
            oss.precision(1);

            oss << "[rank " << me << "] "
                << (where ? where : "")
                << " | RSS=" << rssGB << " GB : "
                << msg << "\n";

            std::cout << oss.str() << std::flush;
        }
    }
    MPI_Barrier(comm);
}

#else  // TEMPEST_MPIOMP not defined

using MPI_Comm = int;
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif

inline void MpiDbgAtF(MPI_Comm, const char*, const char*, ...) {}
inline void MpiDbgF(const char*, ...) {}
inline void MpiDbgAtFOrdered(MPI_Comm, const char*, const char*, ...) {}

#endif // TEMPEST_MPIOMP


#else  // STITCHBLOBS_MPI_DEBUG == 0

// ========================= No-op stubs =========================

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#else
using MPI_Comm = int;
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
#endif
#endif

inline void MpiDbgAtF(MPI_Comm, const char*, const char*, ...) {}
inline void MpiDbgF(const char*, ...) {}
inline void MpiDbgAtFOrdered(MPI_Comm, const char*, const char*, ...) {}

#endif // STITCHBLOBS_MPI_DEBUG

#endif // STITCHBLOBS_MPI_DEBUG_H
