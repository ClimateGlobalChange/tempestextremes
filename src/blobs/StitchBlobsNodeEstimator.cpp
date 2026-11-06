///////////////////////////////////////////////////////////////////////////////
///
///	\file    StitchBlobsNodeEstimator.cpp
///	\version November 6th, 2026
///
///	<remarks>
///		Copyright 2000-2026 Paul Ullrich, Hongyu Chen
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

// Serial pre-flight helper for the MPI build of StitchBlobs.  It samples tagged
// points from the input files, applies StitchBlobs' file-to-rank distribution,
// and reports the minimum number of nodes and an estimated range of total ranks.
// It uses no MPI and is safe to run on a login node.

#include "netcdfcpp.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double BytesPerGiB = 1024.0 * 1024.0 * 1024.0;

// Perlmutter CPU node: memory a Slurm job can use, in GiB (975604 MiB per 2 nodes).
const double DefaultNodeMemGiB = 476.0;

// Target ranks on a node before rounding to whole files; the density of the validated runs.
const int DefaultRanksPerNode = 16;

// Perlmutter CPU node: physical cores, the most ranks placed on a node.
const int CoresPerNode = 128;

// Heap cost of one std::set<int> node: 40-byte tree node in a 48-byte malloc chunk.
const double DefaultBytesPerPoint = 48.0;

// Per-rank cost of the grid and connectivity (2.27 GiB run of one ERA5 file).
const double RankBytesPerGridPoint = 200.0;

// Extra rank-0 cost of gathering tags and the graph (about 10 GiB per 4.96e10 points).
const double RootBytesPerPoint = 0.22;

// Largest share of node memory the fullest node may use at the minimum node count.
const double MinimumNodeFraction = 0.95;

struct Options {
	std::string strInputFileList;
	std::string strVariable;
	double dNodeMemGiB = DefaultNodeMemGiB;
	int nRanksPerNode = DefaultRanksPerNode;
	int nSampleFiles = 24;
	int nSampleTimes = 6;
	double dBytesPerPoint = DefaultBytesPerPoint;
};

struct FileInfo {
	std::string strVarPath;
	long nTimes = 0;
	double dPointsPerTime = 0.0;
};

struct RankRange {
	int nLow = 0;
	int nHigh = 0;
};

void PrintUsage() {
	printf(
		"Usage: StitchBlobsNodeEstimator --in_list <file> --var <name> [options]\n"
		"  --node_mem_gib <x>     usable memory per node (default %.0f, Perlmutter CPU)\n"
		"  --ranks_per_node <n>   target MPI ranks per node (default %d)\n"
		"  --sample_files <n>     files to sample (default 24)\n"
		"  --sample_times <n>     time steps to sample per file (default 6)\n"
		"  --bytes_per_point <x>  memory per tagged point (default %.0f)\n",
		DefaultNodeMemGiB, DefaultRanksPerNode, DefaultBytesPerPoint);
}

Options ParseArguments(int argc, char** argv) {
	Options opt;
	for (int i = 1; i < argc; i++) {
		const std::string strArg(argv[i]);
		if ((strArg == "--help") || (strArg == "-h")) {
			PrintUsage();
			std::exit(0);
		}
		if (i + 1 >= argc) {
			throw std::runtime_error("Missing value for argument " + strArg);
		}
		const char* szValue = argv[++i];
		if (strArg == "--in_list") {
			opt.strInputFileList = szValue;
		} else if (strArg == "--var") {
			opt.strVariable = szValue;
		} else if (strArg == "--node_mem_gib") {
			opt.dNodeMemGiB = std::atof(szValue);
		} else if (strArg == "--ranks_per_node") {
			opt.nRanksPerNode = std::atoi(szValue);
		} else if (strArg == "--sample_files") {
			opt.nSampleFiles = std::atoi(szValue);
		} else if (strArg == "--sample_times") {
			opt.nSampleTimes = std::atoi(szValue);
		} else if (strArg == "--bytes_per_point") {
			opt.dBytesPerPoint = std::atof(szValue);
		} else {
			throw std::runtime_error("Unknown argument " + strArg);
		}
	}

	if (opt.strInputFileList.empty() || opt.strVariable.empty()) {
		throw std::runtime_error("--in_list and --var are required");
	}
	if ((opt.dNodeMemGiB <= 0.0) || (opt.nRanksPerNode <= 0) ||
	    (opt.nSampleFiles <= 0) || (opt.nSampleTimes <= 0) ||
	    (opt.dBytesPerPoint <= 0.0)
	) {
		throw std::runtime_error("Numeric arguments must be positive");
	}
	return opt;
}

std::vector<std::string> ReadInputList(const std::string & strPath) {
	std::ifstream ifs(strPath);
	if (!ifs) {
		throw std::runtime_error("Unable to open input list " + strPath);
	}
	std::vector<std::string> vecLines;
	std::string strLine;
	while (std::getline(ifs, strLine)) {
		if (!strLine.empty() && (strLine.back() == '\r')) {
			strLine.pop_back();
		}
		if (!strLine.empty()) {
			vecLines.push_back(strLine);
		}
	}
	if (vecLines.empty()) {
		throw std::runtime_error("Input list is empty: " + strPath);
	}
	return vecLines;
}

// Each line of an input list may hold several files separated by semicolons.
// Returns the first of them that contains the variable.
FileInfo ReadFileInfo(const std::string & strLine, const std::string & strVariable) {
	std::stringstream ss(strLine);
	std::string strPath;
	while (std::getline(ss, strPath, ';')) {
		NcFile ncfile(strPath.c_str(), NcFile::ReadOnly);
		if (!ncfile.is_valid()) {
			throw std::runtime_error("Unable to open " + strPath);
		}
		if (ncfile.get_var(strVariable.c_str()) == NULL) {
			continue;
		}
		NcDim* dimTime = ncfile.get_dim("time");
		if (dimTime == NULL) {
			throw std::runtime_error("No dimension \"time\" in " + strPath);
		}
		FileInfo info;
		info.strVarPath = strPath;
		info.nTimes = dimTime->size();
		return info;
	}
	throw std::runtime_error("Variable " + strVariable + " not found in: " + strLine);
}

// Counts points StitchBlobs would tag (nonzero values) at evenly spaced times.
double SamplePointsPerTime(
	const FileInfo & info,
	const std::string & strVariable,
	int nSampleTimes,
	long & nGridPoints
) {
	if (info.nTimes == 0) {
		return 0.0;
	}
	NcFile ncfile(info.strVarPath.c_str(), NcFile::ReadOnly);
	NcVar* var = ncfile.get_var(strVariable.c_str());
	if (!ncfile.is_valid() || (var == NULL)) {
		throw std::runtime_error("Unable to read " + strVariable + " from " + info.strVarPath);
	}
	const int nDims = var->num_dims();
	if ((nDims < 2) || (std::string(var->get_dim(0)->name()) != "time")) {
		throw std::runtime_error(
			"Variable " + strVariable + " must have \"time\" as its first dimension");
	}

	std::vector<long> vecCount(nDims);
	std::vector<long> vecStart(nDims, 0);
	vecCount[0] = 1;
	nGridPoints = 1;
	for (int d = 1; d < nDims; d++) {
		vecCount[d] = var->get_dim(d)->size();
		nGridPoints *= vecCount[d];
	}

	const int nSamples = static_cast<int>(std::min<long>(nSampleTimes, info.nTimes));
	std::vector<float> vecData(nGridPoints);
	double dTotal = 0.0;
	for (int s = 0; s < nSamples; s++) {
		vecStart[0] = ((2 * s + 1) * info.nTimes) / (2 * nSamples);
		if (!var->set_cur(vecStart.data()) || !var->get(vecData.data(), vecCount.data())) {
			throw std::runtime_error("Read failed for " + info.strVarPath);
		}
		for (const float dValue : vecData) {
			if (dValue != 0.0f) {
				dTotal += 1.0;
			}
		}
	}
	return dTotal / nSamples;
}

// Same distribution as StitchBlobs: the first (nFiles % nRanks) ranks get one extra file.
void GetRankFileRange(int nFiles, int nRanks, int iRank, int & iBegin, int & iEnd) {
	const int nAverage = nFiles / nRanks;
	const int nRemainder = nFiles % nRanks;
	iBegin = iRank * nAverage + std::min(iRank, nRemainder);
	iEnd = iBegin + nAverage + ((iRank < nRemainder) ? 1 : 0);
}

// Rounds up to an even rank count, or down when that would exceed one rank per file.
int RoundToEvenRanks(int nRanks, int nFiles) {
	if (nRanks % 2 == 0) {
		return nRanks;
	}
	return (nRanks < nFiles) ? nRanks + 1 : std::max(1, nRanks - 1);
}

// A file cannot be split, so the slowest rank is the one with ceil(nFiles / nRanks)
// files.  The range ends are the smallest even rank counts that reach the two whole
// numbers of files per rank around the target; counts between them lower nothing.
RankRange GetRankRange(int nFiles, int nNodes, int nTargetRanksPerNode) {
	const int nMaxRanks = nNodes * CoresPerNode;
	const int nTarget = std::min(nNodes * nTargetRanksPerNode, nMaxRanks);
	const int nFilesAtHigh = std::max(1, nFiles / nTarget);
	const int nFilesAtLow = nFilesAtHigh + 1;

	RankRange range;
	range.nLow = RoundToEvenRanks(
		std::max((nFiles + nFilesAtLow - 1) / nFilesAtLow, std::min(nNodes, nFiles)), nFiles);
	range.nHigh = RoundToEvenRanks((nFiles + nFilesAtHigh - 1) / nFilesAtHigh, nFiles);
	if (range.nHigh > nMaxRanks) {
		range.nHigh = range.nLow;
	}
	return range;
}

// Memory of the fullest node, in GiB.
double EstimateFullestNodeGiB(
	const std::vector<FileInfo> & vecFiles,
	long nGridPoints,
	double dTotalPoints,
	int nNodes,
	int nRanks,
	const Options & opt
) {
	const int nFiles = static_cast<int>(vecFiles.size());
	const double dRankBaseBytes = RankBytesPerGridPoint * static_cast<double>(nGridPoints);

	double dMaxNodeGiB = 0.0;

	// Slurm's block distribution: the first (nRanks % nNodes) nodes hold one extra rank.
	int iRank = 0;
	for (int n = 0; n < nNodes; n++) {
		const int nNodeRanks = nRanks / nNodes + ((n < nRanks % nNodes) ? 1 : 0);

		double dNodeBytes = dRankBaseBytes * nNodeRanks;
		for (int r = 0; r < nNodeRanks; r++, iRank++) {
			int iBegin = 0;
			int iEnd = 0;
			GetRankFileRange(nFiles, nRanks, iRank, iBegin, iEnd);
			for (int f = iBegin; f < iEnd; f++) {
				dNodeBytes += opt.dBytesPerPoint
					* vecFiles[f].dPointsPerTime * static_cast<double>(vecFiles[f].nTimes);
			}
		}
		if (n == 0) {
			dNodeBytes += RootBytesPerPoint * dTotalPoints;
		}
		dMaxNodeGiB = std::max(dMaxNodeGiB, dNodeBytes / BytesPerGiB);
	}
	return dMaxNodeGiB;
}

int Run(int argc, char** argv) {
	const Options opt = ParseArguments(argc, argv);

	// Report NetCDF errors through return values instead of aborting.
	NcError error(NcError::silent_nonfatal);

	const std::vector<std::string> vecLines = ReadInputList(opt.strInputFileList);
	const int nFiles = static_cast<int>(vecLines.size());

	std::vector<FileInfo> vecFiles;
	vecFiles.reserve(nFiles);
	for (const std::string & strLine : vecLines) {
		vecFiles.push_back(ReadFileInfo(strLine, opt.strVariable));
	}

	// Sample evenly spaced files; every other file takes the nearest sample.
	const int nSamples = std::min(opt.nSampleFiles, nFiles);
	std::vector<int> vecSampleFile(nSamples);
	std::vector<double> vecSamplePoints(nSamples);
	long nGridPoints = 0;
	for (int s = 0; s < nSamples; s++) {
		vecSampleFile[s] = static_cast<int>(
			((2 * static_cast<long>(s) + 1) * nFiles) / (2 * nSamples));
		vecSamplePoints[s] = SamplePointsPerTime(
			vecFiles[vecSampleFile[s]], opt.strVariable, opt.nSampleTimes, nGridPoints);
	}
	if (nGridPoints == 0) {
		throw std::runtime_error("Sampled files contain no time steps");
	}

	double dTotalPoints = 0.0;
	for (int f = 0; f < nFiles; f++) {
		int iNearest = 0;
		for (int s = 1; s < nSamples; s++) {
			if (std::abs(vecSampleFile[s] - f) < std::abs(vecSampleFile[iNearest] - f)) {
				iNearest = s;
			}
		}
		vecFiles[f].dPointsPerTime = vecSamplePoints[iNearest];
		dTotalPoints += vecFiles[f].dPointsPerTime * static_cast<double>(vecFiles[f].nTimes);
	}

	// Smallest node count whose fullest node stays under the limit at both range ends.
	const double dLimitGiB = MinimumNodeFraction * opt.dNodeMemGiB;
	for (int n = 1; n <= nFiles; n++) {
		const RankRange range = GetRankRange(nFiles, n, opt.nRanksPerNode);
		const double dFullestGiB = std::max(
			EstimateFullestNodeGiB(vecFiles, nGridPoints, dTotalPoints, n, range.nLow, opt),
			EstimateFullestNodeGiB(vecFiles, nGridPoints, dTotalPoints, n, range.nHigh, opt));
		if (dFullestGiB <= dLimitGiB) {
			printf("Minimum nodes:   %d (fullest node about %.0f GiB of %.0f GiB, %.0f%%)\n",
				n, dFullestGiB, opt.dNodeMemGiB, 100.0 * dFullestGiB / opt.dNodeMemGiB);
			if (range.nLow == range.nHigh) {
				printf("Estimated ranks: %d in total\n", range.nLow);
			} else {
				printf("Estimated ranks: %d to %d in total\n", range.nLow, range.nHigh);
			}
			return 0;
		}
	}
	throw std::runtime_error("A single file does not fit on a node");
}

} // namespace

int main(int argc, char** argv) {
	try {
		return Run(argc, argv);
	} catch (const std::exception & e) {
		fprintf(stderr, "ERROR: %s\n", e.what());
		PrintUsage();
		return 1;
	}
}

///////////////////////////////////////////////////////////////////////////////
