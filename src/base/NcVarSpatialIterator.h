///////////////////////////////////////////////////////////////////////////////
///
///	\file    NcVarSpatialIterator.h
///	\author  Paul Ullrich
///	\version August 7, 2026
///
///	<remarks>
///		Copyright 2026 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _NCVARSPATIALITERATOR_H_
#define _NCVARSPATIALITERATOR_H_

#include "Exception.h"
#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of open NetCDF files.
///	</summary>
class NcVarSpatialIterator {

private:
	///	<summary>
	///		Get or put.
	///	</summary>
	enum OP_GETPUT {
		OP_GET,
		OP_PUT
	};

	///	<summary>
	///		1D SetCur / GetPut operation (shifts get/put specification to argument)
	///	</summary>
	void SetCurGetPut(
		OP_GETPUT op,
		NcVar * var,
		float * pData,
		long l0,
		long s0
	) {
		var->set_cur((long)l0);
		if (op == OP_GET) {
			var->get(pData, s0);
		} else {
			var->put(pData, s0);
		}
	}

	///	<summary>
	///		1D SetCur / GetPut operation (shifts get/put specification to argument)
	///	</summary>
	void SetCurGetPut(
		OP_GETPUT op,
		NcVar * var,
		float * pData,
		long l0,
		long l1,
		long s0,
		long s1
	) {
		if (l0 == (-1)) {
			return SetCurGetPut(op, var, pData, l1, s1);
		}
		var->set_cur(l0, l1);
		if (op == OP_GET) {
			var->get(pData, s0, s1);
		} else {
			var->put(pData, s0, s1);
		}
	}

	///	<summary>
	///		2D GetPut operation (shifts get/put specification to argument)
	///	</summary>
	void SetCurGetPut(
		OP_GETPUT op,
		NcVar * var,
		float * pData,
		long l0,
		long l1,
		long l2,
		long s0,
		long s1,
		long s2
	) {
		if (l0 == (-1)) {
			return SetCurGetPut(op, var, pData, l1, l2, s1, s2);
		}
		var->set_cur(l0, l1, l2);
		if (op == OP_GET) {
			//printf("GET %li %li %li", lTime, lY, lX);
			var->get(pData, s0, s1, s2);
		} else {
			var->put(pData, s0, s1, s2);
		}
	}

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	NcVarSpatialIterator(
		size_t sSpatialDim,
		size_t sBlockSize,
		bool fWholeChunksOnly = true
	) :
		m_sSpatialDim(sSpatialDim),
		m_sBlockSize(sBlockSize),
		m_sCurrentPos(0),
		m_fWholeChunksOnly(fWholeChunksOnly),
		m_fAllBlocksLoaded(false)
	{
		if ((sSpatialDim < 1) || (sSpatialDim > 2)) {
			_EXCEPTIONT("sSpatialDim must be 1 or 2");
		}
	}

private:
	///	<summary>
	///		GetPut the next block of data.
	///	</summary>
	size_t GetPutBlockFloat(
		OP_GETPUT op,
		NcVar * var,
		long lTime,
		float * pData
	) {
		_ASSERT(var != NULL);
		_ASSERT((m_sSpatialDim == 1) || (m_sSpatialDim == 2));
		_ASSERT(pData != NULL);

		long lNumDims = var->num_dims();

		// Read/write in all data at time slice
		if (m_sBlockSize == 0) {
			if ((op == OP_GET) && (m_sCurrentPos != 0)) {
				return 0;
			}
			if ((op == OP_PUT) && (m_sCurrentPos != 0)) {
				_EXCEPTIONT("Logic error");
			}

			m_fAllBlocksLoaded = true;

			if (m_sSpatialDim == 1) {
				_ASSERT(lNumDims == 2);
				long lDimXSize = var->get_dim(lNumDims-1)->size();
				SetCurGetPut(op, var, pData, lTime, 0, 1, lDimXSize);
				return static_cast<size_t>(lDimXSize);

			} else {
				_ASSERT(lNumDims == 3);
				long lDimYSize = var->get_dim(lNumDims-2)->size();
				long lDimXSize = var->get_dim(lNumDims-1)->size();
				SetCurGetPut(op, var, pData, lTime, 0, 0, 1, lDimYSize, lDimXSize);
				return static_cast<size_t>(lDimYSize * lDimXSize);
			}

		// 1D data, limited block size
		} else if (m_sSpatialDim == 1) {
			_ASSERT(((lNumDims == 2) && (lTime >= 0)) || ((lNumDims == 1) && (lTime == -1)));
			long lDimXSize = var->get_dim(lNumDims-1)->size();

			if (m_sCurrentPos * m_sBlockSize >= static_cast<size_t>(lDimXSize)) {
				return 0;
			}
			long lPosX = static_cast<long>(m_sBlockSize * m_sCurrentPos);
			long lSizeX = static_cast<long>(m_sBlockSize);
			if (lPosX + static_cast<long>(m_sBlockSize) > static_cast<size_t>(lDimXSize)) {
				lSizeX = lDimXSize - lPosX;
			}
			SetCurGetPut(op, var, pData, lTime, lPosX, 1, lSizeX);

			if (lPosX + lSizeX == lDimXSize) {
				m_fAllBlocksLoaded = true;
			}

			return lSizeX;

		// 2D data, limited block size
		} else {
			_ASSERT(((lNumDims == 3) && (lTime >= 0)) || ((lNumDims == 2) && (lTime == -1)));
			long lDimYSize = var->get_dim(lNumDims-2)->size();
			long lDimXSize = var->get_dim(lNumDims-1)->size();

			std::vector<size_t> vecChunkSize = var->get_chunk_size();
			_ASSERT((vecChunkSize.size() == 0) || (vecChunkSize.size() == lNumDims));
			for (size_t s = 0; s < vecChunkSize.size(); s++) {
				_ASSERT(vecChunkSize[s] > 0);
			}

			// Data is stored contiguous
			if (vecChunkSize.size() == 0) {
				long lRowsPerBlock = static_cast<long>(m_sBlockSize) / lDimXSize;
				long lPos1D;
				if (m_fWholeChunksOnly) {
					if (lRowsPerBlock == 0) {
						_EXCEPTIONT("Cannot fit full row from NetCDF variable into block; "
							"increase block size");
					}
					lPos1D = static_cast<long>(m_sCurrentPos) * lRowsPerBlock * lDimXSize;
				} else {
					lPos1D = static_cast<long>(m_sBlockSize * m_sCurrentPos);
				}

				long lPosX = lPos1D % lDimXSize;
				long lPosY = lPos1D / lDimXSize;

				if (lPosY >= lDimYSize) {
					m_fAllBlocksLoaded = true;
					return 0;
				}

				// Number of points read
				long lReadPoints = 0;

				// Read/write points starting partially along this row
				if (!m_fWholeChunksOnly && (lPosX != 0)) {
					lReadPoints = std::min(lDimXSize - lPosX, static_cast<long>(m_sBlockSize));
					//printf("Y,X: %li,%li\n", lPosY, lPosX);
					SetCurGetPut(op, var, pData, lTime, lPosY, lPosX, 1, 1, lReadPoints);

					if (lReadPoints == static_cast<long>(m_sBlockSize)) {
						return static_cast<size_t>(lReadPoints);
					} else {
						pData += static_cast<size_t>(lReadPoints);
						lPosX = 0;
						lPosY++;
					}
					if (lPosY == lDimYSize) {
						m_fAllBlocksLoaded = true;
						return static_cast<size_t>(lReadPoints);
					}
				}

				// Read/write full rows
				_ASSERT(static_cast<long>(m_sBlockSize) > lReadPoints);
				long lToReadYSize = (static_cast<long>(m_sBlockSize) - lReadPoints) / lDimXSize;
				if (lToReadYSize != 0) {
					if (lPosY + lToReadYSize >= lDimYSize) {
						lToReadYSize = lDimYSize - lPosY;
					}
					//printf("Reading %li rows", lToReadYSize);
					//printf("..Y,X: %li,0\n", lPosY);
					SetCurGetPut(op, var, pData,
						lTime, lPosY, 0,
						1, lToReadYSize, lDimXSize);
					lReadPoints += lToReadYSize * lDimXSize;
					pData += static_cast<size_t>(lToReadYSize * lDimXSize);
					lPosY += lToReadYSize;
				}

				// Read/write remaining elements in new rows
				if (!m_fWholeChunksOnly && 
				    (lPosY != lDimYSize) &&
				    (lReadPoints < static_cast<long>(m_sBlockSize))
				) {
					_ASSERT(static_cast<long>(m_sBlockSize) - lReadPoints <= lDimXSize)
					//printf("..Y,X: %li,0\n", lPosY);
					SetCurGetPut(op, var, pData,
						lTime, lPosY, 0,
						1, 1, static_cast<long>(m_sBlockSize) - lReadPoints);
					lReadPoints = static_cast<long>(m_sBlockSize);
				}

				//printf(" %li + %li == %li\n", lPos1D, lReadPoints, lDimXSize * lDimYSize);
				if (lPos1D + lReadPoints == lDimXSize * lDimYSize) {
					m_fAllBlocksLoaded = true;
				}

				return static_cast<size_t>(lReadPoints);

			// Data is chunked
			} else {
				long lDimYSize = var->get_dim(lNumDims-2)->size();
				long lDimXSize = var->get_dim(lNumDims-1)->size();

				size_t sChunkSizeY = vecChunkSize[vecChunkSize.size()-2];
				size_t sChunkSizeX = vecChunkSize[vecChunkSize.size()-1];
				size_t sChunkSize1D = sChunkSizeX * sChunkSizeY;

				size_t sChunksX = (lDimXSize + sChunkSizeX - 1) / sChunkSizeX;
				size_t sChunksY = (lDimYSize + sChunkSizeY - 1) / sChunkSizeY;
				size_t sChunks1D = sChunksX * sChunksY;

				// Read whole chunks only
				if (m_fWholeChunksOnly) {
					if (sChunkSize1D > m_sBlockSize) {
						_EXCEPTIONT("Cannot fit full chunk from NetCDF variable into block; "
							"increase block size");
					}

					//TODO: Optimize first chunk identification logic
					// Determine where to start reading/writing this block
					size_t sAccumulatedPoints = 0;
					size_t sCurrentBlock = 0;

					size_t sC1D = 0;
					for (; sC1D < sChunks1D; sC1D++) {
						size_t sCY = sC1D / sChunksX;
						size_t sCX = sC1D % sChunksX;

						size_t sCurrentChunkSizeX =
							std::min(sChunkSizeX,
								static_cast<size_t>(lDimXSize) - sCX * sChunkSizeX);
						size_t sCurrentChunkSizeY =
							std::min(sChunkSizeY,
								static_cast<size_t>(lDimYSize) - sCY * sChunkSizeY);
						size_t sCurrentChunkSize1D =
							sCurrentChunkSizeX * sCurrentChunkSizeY;

						if (sAccumulatedPoints + sCurrentChunkSize1D > m_sBlockSize) {
							sCurrentBlock++;
							sAccumulatedPoints = sCurrentChunkSize1D;
						} else {
							sAccumulatedPoints += sCurrentChunkSize1D;
						}
						if (sCurrentBlock == m_sCurrentPos) {
							break;
						}
					}

					// Read/write chunks until no more chunks fit in block
					sAccumulatedPoints = 0;
					for (; sC1D < sChunks1D; sC1D++) {
						size_t sCY = sC1D / sChunksX;
						size_t sCX = sC1D % sChunksX;

						size_t sCurrentChunkSizeX =
							std::min(sChunkSizeX,
								static_cast<size_t>(lDimXSize) - sCX * sChunkSizeX);
						size_t sCurrentChunkSizeY =
							std::min(sChunkSizeY,
								static_cast<size_t>(lDimYSize) - sCY * sChunkSizeY);
						size_t sCurrentChunkSize1D =
							sCurrentChunkSizeX * sCurrentChunkSizeY;

						if (sAccumulatedPoints + sCurrentChunkSize1D > m_sBlockSize) {
							break;
						}

						// TODO: Read in multiple chunks simultaneously?
/*
						printf("Chunk (%i) %lu / %lu : %li %lu %lu : %lu %lu\n",
							op,
							sC1D, sChunks1D,
							lTime, sCY * sChunkSizeY, sCX * sChunkSizeX,
							sCurrentChunkSizeY, sCurrentChunkSizeX);
*/
						SetCurGetPut(op, var, pData,
							lTime,
							static_cast<long>(sCY * sChunkSizeY),
							static_cast<long>(sCX * sChunkSizeX),
							1,
							sCurrentChunkSizeY,
							sCurrentChunkSizeX);

						sAccumulatedPoints += sCurrentChunkSizeX * sCurrentChunkSizeY;
					}

					if (sC1D == sChunks1D) {
						m_fAllBlocksLoaded = true;
					}

					return sAccumulatedPoints;

				// Read/write partial chunks
				} else {
					_EXCEPTIONT("Not implemented");
/*
					// Determine where to start reading this block
					size_t sCurrentPoint = m_sCurrentPos * m_sBlockSize;

					size_t sAccumulatedPoints = 0;
					size_t sCurrentBlock = 0;

					size_t sCY = 0;
					size_t sCX = 0;
					for (; sCY < sChunksY; sCY++) {
						for (; sCX < sChunksX; sCX++) {
							size_t sCurrentChunkSizeX = sChunkSizeX;
							if (sCX == sChunksX-1) {
								sCurrentChunkSizeX = sEdgeChunkSizeX;
							}
							size_t sCurrentChunkSizeY = sChunkSizeY;
							if (sCY == sChunksY-1) {
								sCurrentChunkSizeY = sEdgeChunkSizeY;
							}
							size_t sCurrentChunkSize1D =
								sCurrentChunkSizeX * sCurrentChunkSizeY;

							if (sAccumulatedPoints + sCurrentChunkSize1D > sCurrentPoint) {
								sCurrentBlock = m_sCurrentPos;
								break;
							}

							sAccumulatedPoints += sCurrentChunkSize1D;
						}
						if (sCurrentBlock == m_sCurrentPos) {
							break;
						}
					}

					size_t sToReadPoints = m_sBlockSize;
					for (; sCY < sChunksY; sCY++) {
						for (; sCX < sChunksX; sCX++) {
							size_t sCurrentChunkSizeX = sChunkSizeX;
							if (sCX == sChunksX-1) {
								sCurrentChunkSizeX = sEdgeChunkSizeX;
							}
							size_t sCurrentChunkSizeY = sChunkSizeY;
							if (sCY == sChunksY-1) {
								sCurrentChunkSizeY = sEdgeChunkSizeY;
							}
							size_t sCurrentChunkSize1D =
								sCurrentChunkSizeX * sCurrentChunkSizeY;

							if (sToReadPoints > sCurrentChunkSize1D) {
								var->set_cur(1,
									static_cast<long>(sCX * sChunkSizeX), 
									static_cast<long>(sCY * sChunkSizeY));
							} else {
							}
						}
						if (sCurrentBlock == m_sCurrentPos) {
							break;
						}
					}
*/
				}
			}
		}

		return 0;
	}

public:
	///	<summary>
	///		Get the next block of data.
	///	</summary>
	size_t GetBlockFloat(
		NcVar * varIn,
		long lTime,
		float * dData
	) {
		return GetPutBlockFloat(OP_GET, varIn, lTime, dData);
	}

	///	<summary>
	///		Get the next block of data.
	///	</summary>
	size_t PutBlockFloat(
		NcVar * varIn,
		long lTime,
		float * dData
	) {
		return GetPutBlockFloat(OP_PUT, varIn, lTime, dData);
	}

	///	<summary>
	///		Go to next block.
	///	</summary>
	bool NextBlock() {
		if (m_fAllBlocksLoaded) {
			return false;
		}
		m_sCurrentPos++;
		return true;
	}

	///	<summary>
	///		Reset the block counter.
	///	</summary>
	void Reset() {
		m_sCurrentPos = 0;
		m_fAllBlocksLoaded = false;
	}

	///	<summary>
	///		Get the current position.
	///	</summary>
	size_t GetCurrentBlockIx() const {
		return m_sCurrentPos;
	}

private:
	///	<summary>
	///		Spatial dimension of grid.
	///	</summary>
	size_t m_sSpatialDim;

	///	<summary>
	///		Number of data points per block.
	///	</summary>
	size_t m_sBlockSize;

	///	<summary>
	///		Current position.
	///	</summary>
	size_t m_sCurrentPos;

	///	<summary>
	///		Only read whole chunks.
	///	</summary>
	bool m_fWholeChunksOnly;

	///	<summary>
	///		All blocks loaded.
	///	</summary>
	bool m_fAllBlocksLoaded;

};

#endif

