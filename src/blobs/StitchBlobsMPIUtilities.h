///////////////////////////////////////////////////////////////////////////////
///
///	\file    BlobUtilities.h
///	\author  Hongyu Chen
///	\version January 21st, 2026
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


#ifndef STITCHBLOBS_MPI_UTILITIES_H
#define STITCHBLOBS_MPI_UTILITIES_H

// This header is intended to be PRIVATE to StitchBlobs.cpp.
// It contains MPI-only helpers and MPI op classes.
// It should only be included when TEMPEST_MPIOMP is enabled.

#if defined(TEMPEST_MPIOMP)

#include <mpi.h>

// Standard library headers that are commonly needed by the MPI helpers.
// Keep this minimal; add includes only when a compile error tells you to.
#include <algorithm>   // std::min, std::max
#include <array>
#include <cmath>       // std::ceil, std::log2
#include <map>
#include <set>
#include <string>
#include <utility>     // std::pair
#include <vector>
#include <mutex>

#include "Exception.h"   // _EXCEPTION1, _EXCEPTIONT

// IMPORTANT: This header must be included only from StitchBlobs.cpp
// and only AFTER struct Tag is defined and GetMPI_Tag_type() is declared.


namespace stitchblobs {
namespace mpi {
///	<summary>
///		A class used for MPI collective communication for vecAllBlobTags
///		(gather and scatter) between processors)
///		This class will serialize the vecAllBlobTags into an 1-D array
///		and deserialize the received 1-D vecAllBlobTags into the original
///		2-D vectors after the communication. 
///		This class can also generate and return the setAllTags based on the local vecAllBlobTags.
///	</summary>
class TagCollectiveOP {
private:
    /// <summary> Tag for MPI communication for sending/receiving Tag </summary>
    int gather_tag = 107;

    /// <summary> Tag for MPI communication for sending/receiving Tag index </summary>
    int gather_tag_index = 108;

    /// <summary> The MPI Communicator </summary>
    MPI_Comm m_comm;


    /// <summary>
    ///   Flag marking whether the flat buffers are currently in serialized state.
    ///   (false initially; true after Serialize(); false after Deserialize()).
    /// </summary>
    bool serializedFlag;

    /// <summary>
    ///   Optional non-owning view of the input 2-D tags. If set, Serialize() reads from *in_;
    ///   otherwise it reads from the owning _vecAllBlobTags.
    /// </summary>
    const std::vector<std::vector<Tag>>* in_ = nullptr;

    /// <summary>
    ///   Serialize the std::vector<std::vector<Tag>> into a 1-D array (tags)
    ///   and an index array (prefix offsets, length = ncols+1).
    /// </summary>
    void Serialize() {
        const auto& src = (in_ ? *in_ : _vecAllBlobTags);

        // Pre-size to minimize reallocations
        size_t total = 0;
        serialVecAllBlobTags.clear();
        serialVecAllBlobTags_index.clear();
        serialVecAllBlobTags_index.reserve(src.size() + 1);
        for (const auto& v : src) total += v.size();
        serialVecAllBlobTags.reserve(total);

        int curIndx = 0;
        serialVecAllBlobTags_index.push_back(curIndx);
        for (size_t i = 0; i < src.size(); ++i) {
            // append this column
            for (size_t j = 0; j < src[i].size(); ++j) {
                serialVecAllBlobTags.push_back(src[i][j]);
                ++curIndx;
            }
            serialVecAllBlobTags_index.push_back(curIndx);
        }
        serializedFlag = true;
    }

    /// <summary>
    ///   Deserialize the local 1-D arrays into the 2-D vector + set.
    ///   After deserialization, the flat buffers are cleared and released.
    /// </summary>
    void Deserialize() {
        desirialVecAllBlobTags.clear();
        if (serialVecAllBlobTags_index.size() > 0) {
            desirialVecAllBlobTags.reserve(serialVecAllBlobTags_index.size() - 1);
        }

        for (size_t i = 0; i + 1 < serialVecAllBlobTags_index.size(); ++i) {
            const int startIndx = serialVecAllBlobTags_index[i];
            const int endIndx   = std::min(serialVecAllBlobTags_index[i+1],
                                           static_cast<int>(serialVecAllBlobTags.size()));
            std::vector<Tag> curVecBlobTags;
            curVecBlobTags.reserve(endIndx - startIndx);
            for (int k = startIndx; k < endIndx; ++k) {
                curVecBlobTags.push_back(serialVecAllBlobTags[k]);
                deserialSetAllTags.insert(serialVecAllBlobTags[k]);
            }
            desirialVecAllBlobTags.push_back(std::move(curVecBlobTags));
        }

        // release flat buffers
        serialVecAllBlobTags.clear();            serialVecAllBlobTags.shrink_to_fit();
        serialVecAllBlobTags_index.clear();      serialVecAllBlobTags_index.shrink_to_fit();
        serializedFlag = false;
    }

protected:
    /// <summary> Owning fallback storage for the input 2-D vecAllBlobTags (used if in_ == nullptr) </summary>
    std::vector<std::vector<Tag>> _vecAllBlobTags;

    /// <summary> Serialized 1-D tags buffer </summary>
    std::vector<Tag> serialVecAllBlobTags;

    /// <summary> Serialized 1-D index buffer (length = ncols+1) </summary>
    std::vector<int> serialVecAllBlobTags_index;

    /// <summary> Output setAllTags after deserialization (only used on processor 0) </summary>
    std::set<Tag> deserialSetAllTags;

    /// <summary> Deserialized 2-D vecAllBlobTags after communication (only used on processor 0) </summary>
    std::vector<std::vector<Tag>> desirialVecAllBlobTags;

    /// <summary> Per-rank total Tag counts (sum of sizes of all local columns) </summary>
    std::vector<int> vecScatterCounts;

    /// <summary>
    ///  Per-rank column counts (number of columns), NOT (ncols+1).
    ///  The +1 sentinel is handled locally after Scatter.
    /// </summary>
    std::vector<int> vecScatterCounts_index;

public:
    /// <summary>
    ///   Constructor: stores a non-owning view of vecAllBlobTags. Uses GetMPI_Tag_type() (cached process-wide datatype) when communicating.
    ///   No deep copies are made here.
    /// </summary>
    TagCollectiveOP(MPI_Comm communicator, const std::vector<std::vector<Tag>>& vecAllBlobTags)
    : m_comm(communicator), serializedFlag(false), in_(&vecAllBlobTags) {}

    /// <summary> Destructor. </summary>
    ~TagCollectiveOP() = default;


    /// <summary>
    ///   Hypercube gather that collects all local vecAllBlobTags on rank 0.
    ///   Rank!=0 ranks exit early once they’ve sent their payloads.
    /// </summary>
    void Gather() {
        int rank = -1, size = 1;
        MPI_Comm_rank(m_comm, &rank);
        MPI_Comm_size(m_comm, &size);

        Serialize(); // serialize current local input (non-owning view)

        // ceil(log2(size)) stages
        const int d = (size <= 1) ? 0 : static_cast<int>(std::ceil(std::log2(static_cast<double>(size))));

        for (int j = 0; j < d; ++j) {
            const int bit = (1 << j);

            if ((rank & bit) != 0) {
                // sender: send both payloads and return
                const int destRank = rank ^ bit;
                if (destRank >= size) continue;

                MPI_Send(serialVecAllBlobTags.data(),
                        static_cast<int>(serialVecAllBlobTags.size()),
                        GetMPI_Tag_type(), destRank, gather_tag, m_comm);


                MPI_Send(serialVecAllBlobTags_index.data(),
                         static_cast<int>(serialVecAllBlobTags_index.size()),
                         MPI_INT, destRank, gather_tag_index, m_comm);

                return; // done on this rank after sending
            } else {
                // receiver: append sender's payloads
                const int sourceRank = rank ^ bit;
                if (sourceRank >= size) continue;

                // --- Tags: probe -> resize once -> receive directly into the back
                MPI_Status status;
                int recvCount = 0;
                MPI_Probe(sourceRank, gather_tag, m_comm, &status);
                MPI_Get_count(&status, GetMPI_Tag_type(), &recvCount);


                const int offT = static_cast<int>(serialVecAllBlobTags.size());
                serialVecAllBlobTags.resize(offT + recvCount);
                MPI_Recv(serialVecAllBlobTags.data() + offT, recvCount, GetMPI_Tag_type(),
                sourceRank, gather_tag, m_comm, MPI_STATUS_IGNORE);


                // --- Index: probe -> receive into tmp -> append [1..end) with offset
                MPI_Status status_index;
                int recvCount_index = 0;
                MPI_Probe(sourceRank, gather_tag_index, m_comm, &status_index);
                MPI_Get_count(&status_index, MPI_INT, &recvCount_index);

                std::vector<int> tmp(recvCount_index);
                MPI_Recv(tmp.data(), recvCount_index, MPI_INT,
                         sourceRank, gather_tag_index, m_comm, MPI_STATUS_IGNORE);

                const int curLocalTagSize =
                    serialVecAllBlobTags_index.empty()
                      ? 0
                      : serialVecAllBlobTags_index.back();

                // append indices [1..end) with base offset
                const size_t oldSz = serialVecAllBlobTags_index.size();
                serialVecAllBlobTags_index.resize(oldSz + std::max(0, recvCount_index - 1));
                for (int k = 1; k < recvCount_index; ++k) {
                    serialVecAllBlobTags_index[oldSz + (k - 1)] = tmp[k] + curLocalTagSize;
                }
            }
        }
    }

    /// <summary> Return the gathered setAllTags (only rank 0 should call). </summary>
    std::set<Tag> GetGatheredSetAllTags() {
        int rank = -1;
        MPI_Comm_rank(m_comm, &rank);
        if (rank == 0) {
            if (serializedFlag) {
                Deserialize();
            }
            return std::move(deserialSetAllTags);
        } else {
            _EXCEPTIONT("Only processor 0 should call GetGatheredSetAllTags().");
        }
    }

    /// <summary>
    ///   Return the gathered/scattered vecAllBlobTags.
    ///   GatheredFlag = 1 → returning global (only rank 0 valid).
    ///   GatheredFlag = 0 → returning local (all ranks valid).
    /// </summary>
    std::vector<std::vector<Tag>> GetUnserialVecAllTags(int GatheredFlag) {
        int rank = -1;
        MPI_Comm_rank(m_comm, &rank);
        if ((GatheredFlag == 1 && rank == 0) || GatheredFlag == 0) {
            if (serializedFlag) {
                Deserialize();
            }
            return std::move(desirialVecAllBlobTags);
        } else {
            _EXCEPTIONT("Only processor 0 should call GetUnserialVecAllTags().");
        }
    }

    /// <summary>
    ///   Gather per-rank sizes from the ORIGINAL (unexchanged) layout.
    ///   On rank 0, fills vecScatterCounts (tags) and vecScatterCounts_index (ncols).
    ///   On other ranks, just updates the non-owning view to the provided input (no copy).
    /// </summary>
    void GatherTagCounts(const std::vector<std::vector<Tag>>& vecAllBlobTags) {
        int rank = -1, size = 1;
        MPI_Comm_rank(m_comm, &rank);
        MPI_Comm_size(m_comm, &size);

        int curSize = 0;
        for (const auto& vecBlobTags : vecAllBlobTags) { // no copies
            curSize += static_cast<int>(vecBlobTags.size());
        }
        const int localCols = static_cast<int>(vecAllBlobTags.size());

        vecScatterCounts.resize(size);
        vecScatterCounts_index.resize(size);

        if (rank == 0) {
            MPI_Gather(&curSize,   1, MPI_INT, vecScatterCounts.data(),       1, MPI_INT, 0, m_comm);
            MPI_Gather(&localCols, 1, MPI_INT, vecScatterCounts_index.data(), 1, MPI_INT, 0, m_comm);
        } else {
            // point our non-owning view at the caller's ORIGINAL layout (no deep copy)
            in_ = &vecAllBlobTags;
            MPI_Gather(&curSize,   1, MPI_INT, nullptr, 0, MPI_INT, 0, m_comm);
            MPI_Gather(&localCols, 1, MPI_INT, nullptr, 0, MPI_INT, 0, m_comm);
        }
    }

    /// <summary>
    ///   Scatter the (updated) global vecAllBlobTags on rank 0 back to all ranks.
    ///   Counts for tags and columns come from GatherTagCounts().
    ///   Column indices (ncols entries per rank) are scattered, and each rank fixes the +1 sentinel locally.
    /// </summary>
    void Scatter() {
        int rank = -1, size = 1;
        MPI_Comm_rank(m_comm, &rank);
        MPI_Comm_size(m_comm, &size);

        // Serialize current view on every rank:
        //  - rank 0 serializes the UPDATED GLOBAL
        //  - others serialize their ORIGINAL local (pointed to by in_)
        Serialize();

        if (rank == 0) {
            // Build counts / displacements for tags (elements) and index (ncols)
            std::vector<int> arrayScatterCounts(size);
            std::vector<int> arrayScatterDisplacements(size);
            _ASSERT(vecScatterCounts.size() > 0);
            arrayScatterCounts[0]        = vecScatterCounts[0];
            arrayScatterDisplacements[0] = 0;
            for (int i = 1; i < static_cast<int>(vecScatterCounts.size()); ++i) {
                arrayScatterCounts[i]        = vecScatterCounts[i];
                arrayScatterDisplacements[i] = arrayScatterDisplacements[i - 1] + vecScatterCounts[i - 1];
            }

            std::vector<int> arrayScatterCounts_index(size);
            std::vector<int> arrayScatterDisplacements_index(size);
            _ASSERT(vecScatterCounts_index.size() > 0);
            arrayScatterCounts_index[0]        = vecScatterCounts_index[0]; // ncols (NOT ncols+1)
            arrayScatterDisplacements_index[0] = 0;
            for (int i = 1; i < static_cast<int>(vecScatterCounts_index.size()); ++i) {
                arrayScatterCounts_index[i]        = vecScatterCounts_index[i];
                arrayScatterDisplacements_index[i] = arrayScatterDisplacements_index[i - 1] + vecScatterCounts_index[i - 1];
            }

            // Scatter without cloning global buffers (MPI_IN_PLACE on root)
            // Tags
            MPI_Scatterv(
                serialVecAllBlobTags.data(),
                arrayScatterCounts.data(),
                arrayScatterDisplacements.data(),
                GetMPI_Tag_type(),
                MPI_IN_PLACE, 0, GetMPI_Tag_type(), 0, m_comm);


            // Index (ncols per rank; each rank will append/set the +1 sentinel)
            MPI_Scatterv(
                serialVecAllBlobTags_index.data(),
                arrayScatterCounts_index.data(),
                arrayScatterDisplacements_index.data(),
                MPI_INT,
                MPI_IN_PLACE, 0, MPI_INT, 0, m_comm);

            // Keep only rank 0 slice locally
            serialVecAllBlobTags.resize(arrayScatterCounts[0]);
            serialVecAllBlobTags_index.resize(arrayScatterCounts_index[0]);
        } else {
            // Non-root: receive only; sizes come from our local serialized view
            const int localTagSize        = static_cast<int>(serialVecAllBlobTags.size());
            const int localTagSize_index  = static_cast<int>(serialVecAllBlobTags_index.size()); // this is (ncols+1)

            serialVecAllBlobTags.clear();
            serialVecAllBlobTags.resize(localTagSize);
            MPI_Scatterv(nullptr, nullptr, nullptr, GetMPI_Tag_type(),
                        serialVecAllBlobTags.data(), localTagSize, GetMPI_Tag_type(),
                        0, m_comm);


            serialVecAllBlobTags_index.clear();
            serialVecAllBlobTags_index.resize(localTagSize_index);
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT,
                         serialVecAllBlobTags_index.data(), localTagSize_index, MPI_INT,
                         0, m_comm);
            // NOTE: root sends only ncols indices; we have an extra slot for the sentinel we set below.
        }

        // Localize index to [0 .. localTagCount] and fix sentinel (UNCHANGED semantics)
        if (!serialVecAllBlobTags_index.empty()) {
            const int base = serialVecAllBlobTags_index[0];
            serialVecAllBlobTags_index[0] = 0;
            for (size_t i = 1; i + 1 < serialVecAllBlobTags_index.size(); ++i) {
                serialVecAllBlobTags_index[i] -= base;
            }
            if (rank == 0) {
                serialVecAllBlobTags_index.push_back(static_cast<int>(serialVecAllBlobTags.size()));
            } else {
                serialVecAllBlobTags_index.back() = static_cast<int>(serialVecAllBlobTags.size());
            }
        }

        // hygiene
        serialVecAllBlobTags.shrink_to_fit();
        serialVecAllBlobTags_index.shrink_to_fit();

        // still serialized (ready for caller to Deserialize() via Get*() accessors)
        serializedFlag = true;
    }
};


template <class ColumnVec2D>
void StripHalosInPlace(ColumnVec2D& cols, int rank, int size) {
    // Only even ranks received halos. Odd ranks were unchanged by the exchange ops.
    if (rank % 2 != 0) return;

    if (cols.empty()) return;

    // interior even ranks: [Lhalo] local... [Rhalo]
    if (rank != 0 && rank != size - 1) {
        if (cols.size() >= 2) {
            cols.erase(cols.begin());                 // drop left halo
            cols.pop_back();                          // drop right halo
        }
    }
    // rank 0 (even): local... [Rhalo]
    else if (rank == 0) {
        if (cols.size() >= 1) {
            cols.pop_back();                          // drop right halo only
        }
    }
    // rank size-1 (could be even if size is even): [Lhalo] local...
    else if (rank == size - 1) {
        if (cols.size() >= 1) {
            cols.erase(cols.begin());                // drop left halo only
        }
    }
}



///	<summary>
///		Enumerator of exchange directions.
///	</summary>
typedef enum {
	DIR_LEFT = 0, 
	DIR_RIGHT = 1
} CommDirection;

///	<summary>
///		A Class used for exchanging vecAllBlobsTag between processors
///		The exchange process will first update the exchanged Tags.time to the actual global time
///		And then start the exchange process.
///	</summary>
class TagExchangeOP {
	private:

		///	<summary>
		///		Tag for MPI communication for blob Tag
		///	</summary>
		int tag = 100;

		
		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

		///	<summary>
		///		The tool function that uses the prefix sum algorithmn to assign global time for each Tag
		///	</summary>		
		void UpdateTime(){
			int err, rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			gloablTimeIndx.clear();                         
			gloablTimeIndx.resize(size);                    
			const int localEndTime = static_cast<int>(in_->size());   
			int prefixTotal = 0;
			int rc = MPI_Scan(&localEndTime, &prefixTotal, 1, MPI_INT, MPI_SUM, m_comm);
			if (rc != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Scan failed (code %i)", rc); 
			}

			// do NOT copy sizes repeatedly; directly size the per-rank window
			// and compute prefix sums to get global start/end for THIS rank.
			const int globalStart = prefixTotal - localEndTime;  
			const int globalEnd   = prefixTotal;
			gloablTimeIndx[rank].resize(2);
			gloablTimeIndx[rank][0] = globalStart;
			gloablTimeIndx[rank][1] = globalEnd;

			// Gather (start,end) for ALL ranks so downstream routines
			// that rely on gloablTimeIndx having full coverage are correct.
			int mypair[2] = { globalStart, globalEnd };                 
			std::vector<int> allpairs(2*size);                          
			rc = MPI_Allgather(mypair, 2, MPI_INT,                      
				allpairs.data(), 2, MPI_INT, m_comm);                   
			if (rc != MPI_SUCCESS) {                                    
				_EXCEPTION1("MPI_Allgather failed in UpdateTime (code %i)", rc);
			}
			for (int r = 0; r < size; ++r) {                            
				gloablTimeIndx[r].resize(2);
				gloablTimeIndx[r][0] = allpairs[2*r+0];
				gloablTimeIndx[r][1] = allpairs[2*r+1];
			}

			// mutate in-place instead of updating a private copy:
			// set Tag.time for each local column to its global index.
			int globalTime = globalStart;
			for (int i = 0; i < static_cast<int>(in_->size()); ++i) {
				std::vector<Tag> & col = (*in_)[static_cast<size_t>(i)];   
				for (size_t j = 0; j < col.size(); ++j) {
					col[j].time = globalTime;
				}
				++globalTime;
			}
			_ASSERT(globalTime == globalEnd);    
		}


	protected:
		///	<summary>
		///		The initial vecAllBlobTags that need to be exchanged
		///		make this a non-owning pointer (avoids deep copy)
		///	</summary>
    	std::vector<std::vector<Tag>>* in_ = nullptr; 



		///	<summary>
		///		The vecAllBlobTags after the exchange. it's column number is _vecAllBlobTags' column number plus 2 (left and right) 
		///		except for p0 and pn-1 (for these two processors, the column number is _vecAllBlobTags' column number plus 1).
		///	</summary>
		std::vector< std::vector<Tag>> exchangedvecAllBlobTags;


		///	<summary>
		///		The buffer for vecAllBlobTags that will be received
		///		recvTags[0] is the left vector and recvTags[1] is the right vector
		///	</summary>
		std::vector<std::vector<Tag>> recvTags;

		///	<summary>
		///		The array (size is nMPISize) records the start(inclusive) and end (exclusive) global time index in each processor
		///		start: gloablTimeIndx[p_i][0]
		///		end: gloablTimeIndx[p_i][1]
		///	</summary>	
		std::vector<std::vector<int>> gloablTimeIndx;	

	public:

		///	<summary>
		///		Construct the Operator with vecAllBlobTags
		///     take a NON-const reference so we can update Tag.time in-place (no deep copy).
		///     Matches the memory-minimizing style used in BlobsExchange (non-owning view).
		///     Also constructs the derived MPI_Datatype for Tag and commits it.
		///	</summary>
        TagExchangeOP(MPI_Comm communicator,
                    std::vector<std::vector<Tag>> & vecAllBlobTags)
        : tag(100),
        m_comm(communicator),
        in_(&vecAllBlobTags)
        {
        }



		///	<summary>
		///		Destructor
		///	</summary>
		~TagExchangeOP(){
			MPIrequests.clear();
			MPIstatuses.clear();
		}

		/// <summary>
		///     Return the original unexchanged vecAllBlobTags
		///     return a const view (no copy) 
		/// </summary>
		const std::vector<std::vector<Tag>>& ViewOriginalVecAllBlobTags() const noexcept { 
			return *in_;                                                                     
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the TagExchangeOP should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			
			//First update all processors' Tag.time to the actual global time
			this->UpdateTime();

			//----------------------Send Tags data first----------------------

			// size recv vectors here 
			recvTags.resize(2);

			if (rank % 2 == 0) {                           
				if (rank == 0 || rank == size - 1) {
					exchangedvecAllBlobTags.resize(in_->size() + 1);
				} else {
					exchangedvecAllBlobTags.resize(in_->size() + 2);
				}
			} else {
				// odd ranks: delay allocation; EndExchange() will just swap with *in_
			}

			// Send data
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						//Rank 0 Do Nothing
						destRank = -1;
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						//Rank n-1 Do Nothing
						destRank = -1;
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
						continue;
				}

				//----------------------Send sendBlobs----------------------

				// Only the odd number processors will send out the data
				if ((rank % 2) != 0) {
					MPI_Request request;
					const std::vector<Tag>& src =
						(dir == DIR_LEFT) ? in_->front() : in_->back();   
					const int count = static_cast<int>(src.size());
					// zero-count safety: pass nullptr when count==0 
                    int result = MPI_Isend(count ? const_cast<Tag*>(src.data()) : nullptr,
                                        count, GetMPI_Tag_type(),
                                        destRank, tag, m_comm, &request);

					if (result != MPI_SUCCESS) {
						_EXCEPTION1("The MPI routine MPI_Isend failed (code %i)", result);
					}
					MPIrequests.emplace_back(request);
					MPIstatuses.emplace_back();
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank = -1;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					// Rank 0 will not receive from the left		 
					if (rank == 0) {			
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					// rank n-1 will not receive from the right
					if (rank == size - 1) {
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}
				if (sourceRank > size - 1) {
						continue;
				}



				//----------------------Receive----------------------	
				// Only the even-numbered processors receive the blobs
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;

					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count(&status, GetMPI_Tag_type(), &recvCount);
					recvTags[dir].resize(recvCount);

					int result =
                        MPI_Irecv(recvTags[dir].data(), static_cast<int>(recvTags[dir].size()),
                                GetMPI_Tag_type(),
                                sourceRank, tag, m_comm, &request);

					if (result != MPI_SUCCESS) {
						_EXCEPTION1("The MPI routine MPI_Isend failed (code %i)", result);
					}
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());
				}
			}

		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the TagExchangeOP without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for all Irecv to complete
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);


			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}

			MPIrequests.clear();
			MPIstatuses.clear();

			// Pack the data into the vecAllBlobTags
			// Only need to pack data for the prime number processors 
			if ((rank % 2) == 0) {
				if (rank == 0) {
					// [FIX] copy local columns so *in_ stays intact for downstream
					for (int i = 0; i < static_cast<int>(exchangedvecAllBlobTags.size()) - 1; ++i) {
						exchangedvecAllBlobTags[i] = (*in_)[static_cast<size_t>(i)];               // [FIX]
					}
					// move the received right boundary in (we own recvTags)
					exchangedvecAllBlobTags.back() = std::move(recvTags[DIR_RIGHT]); 
				} else if (rank == size - 1) {
					// move the received left boundary in
					exchangedvecAllBlobTags.front() = std::move(recvTags[DIR_LEFT]); 
					// [FIX] copy local columns (shifted by one)
					for (int i = 1; i < static_cast<int>(exchangedvecAllBlobTags.size()); ++i) {
						exchangedvecAllBlobTags[i] = (*in_)[static_cast<size_t>(i - 1)];           // [FIX]
					}
				} else {
					//  move both received edges in
					exchangedvecAllBlobTags.front() = std::move(recvTags[DIR_LEFT]); 
					// [FIX] copy local middle columns
					for (int i = 1; i < static_cast<int>(exchangedvecAllBlobTags.size()) - 1; ++i) {
						exchangedvecAllBlobTags[i] = (*in_)[static_cast<size_t>(i - 1)];           // [FIX]
					}
					exchangedvecAllBlobTags.back() = std::move(recvTags[DIR_RIGHT]);  
				}
			} else {
				// Odd ranks: unchanged view; keep original intact
				// [FIX] copy instead of swap, so *in_ remains available later
				exchangedvecAllBlobTags = *in_;                                                   // [FIX]
			}

			recvTags.clear();      recvTags.shrink_to_fit();   
			gloablTimeIndx.clear();                            
			gloablTimeIndx.shrink_to_fit();                    

		}

		/// <summary>
		///     Move-only handoff to avoid copies 
		/// </summary>
		std::vector<std::vector<Tag>> TakeExchangedVecAllBlobTags() && noexcept { 
			return std::move(this->exchangedvecAllBlobTags);                       
		}

};


///	<summary>
///		Class for exchanging vecAllBlobBoxesDeg among processors
///	</summary>
class BlobBoxesDegExchangeOP {
	private:
		///	<summary>
		///		Tag for MPI communication for LatlonBox<double>
		///	</summary>
		int tag = 103;

		// ///	<summary>
		// ///		The MPI Datatype for LatlonBox<double>
		// ///	</summary>
		// MPI_Datatype MPI_LatonBox_double_type;

		// ///	<summary>
		// ///		The MPI Datatype for double[2]
		// ///	</summary>
		// MPI_Datatype MPI_doubleArray;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

	protected:
		///	<summary>
		///		The initial vecAllBlobBoxesDeg that will be sent, non-owning pointer
		///	</summary>
		std::vector<std::vector<LatLonBox<double>>>* in_ = nullptr;

		///	<summary>
		///		The vecAllBlobBoxesDeg that is after the exchange. it's column number is _vecAllBlobBoxesDegs' column number plus 2 (left and right) 
		///		except for p0 and pn-1 (for these two processors, the column number is _vecAllBBlobBoxesDegs' column number plus 1).
		///	</summary>
		std::vector<std::vector<LatLonBox<double>>> exchangedvecAllBlobBoxesDeg;

		///	<summary>
		///		The buffer for vecAllBlobBoxesDeg that will be received
		///		sendBlobBoxesDeg[0] is the left vector and sendBlobBoxesDeg[1] is the right vector
		///	</summary>
		std::vector<std::vector<LatLonBox<double>>> recvBlobBoxesDeg;

	public:

		///	<summary>
		///		Construct the Operator with vecAllBlobBoxesDeg
		///		It will contruct the this->m_comm and this->_vecAllBlobBoxesDeg based on the input communicator and vecAllBlobBoxesDeg
		///	</summary>
		BlobBoxesDegExchangeOP(
			MPI_Comm communicator,
			std::vector<std::vector<LatLonBox<double>>>& vecAllBlobBoxesDeg // [CHANGE] const& view
		) noexcept {
			this->in_ = &vecAllBlobBoxesDeg;
			this->m_comm = communicator;

		//########################### Notes for  derived MPI datatype of the LatLonBox<double>(Hongyu Chen) ############################################################################ 
		//1.  "Because vector<bool> holds bits instead of bools, it can't return a bool& from its indexing operator or iterator dereference" (src: https://isocpp.org/blog/2012/11/on-vectorbool)
		//    Therefore, creating the userdefined datatype for LatLonBox<double> and then use vector.data() to send/recv like TagExhangeOP is not working here
		//2.  Now we're using the MPI_BYTE to manually calculate the send/recv buffer size in byte here, which is working currently
		//3.  If the program breaks again here, please consider going to the BlobUtilities.h and modify the constructer at line 232 and line 236 according to the description there.

			// //Create an MPI datatype for the LatLonBox<double>:
			// //First create the datatype for double[2]
			// MPI_Type_contiguous	(2,	MPI_DOUBLE,&MPI_doubleArray);	
			// MPI_Type_commit	(&MPI_doubleArray);
			// //Then use this doubleArray to construct LatlonBox<double>
			// LatLonBox<double> sampleBox;
			// int LatlonBoxFieldCount = 5;
			// MPI_Datatype LatlonBox_typesig[5] = {MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_DOUBLE, MPI_doubleArray, MPI_doubleArray};
			// int LatlonBox_block_lengths[5] = {1,1,1,1,1};
			// MPI_Aint LatlongBox_displacements[5];
			// MPI_Aint base_address;
			// MPI_Get_address(&sampleBox.is_null, &LatlongBox_displacements[0]);
			// MPI_Get_address(&sampleBox.lon_periodic, &LatlongBox_displacements[1]);
			// MPI_Get_address(&sampleBox.lon_width, &LatlongBox_displacements[2]);
			// MPI_Get_address(&sampleBox.lon, &LatlongBox_displacements[3]);
			// MPI_Get_address(&sampleBox.lat, &LatlongBox_displacements[4]);
			// LatlongBox_displacements[0] = MPI_Aint_diff(LatlongBox_displacements[0], base_address);
			// LatlongBox_displacements[1] = MPI_Aint_diff(LatlongBox_displacements[1], base_address);
			// LatlongBox_displacements[2] = MPI_Aint_diff(LatlongBox_displacements[2], base_address);
			// LatlongBox_displacements[3] = MPI_Aint_diff(LatlongBox_displacements[3], base_address);
			// LatlongBox_displacements[4] = MPI_Aint_diff(LatlongBox_displacements[4], base_address);
			// MPI_Type_create_struct(LatlonBoxFieldCount, LatlonBox_block_lengths, LatlongBox_displacements, LatlonBox_typesig, &MPI_LatonBox_double_type);
			// MPI_Type_commit(&MPI_LatonBox_double_type);

			//########################### End Notes for  derived MPI datatype of the LatLonBox<double>(Hongyu Chen) ############################################################################ 

			recvBlobBoxesDeg.resize(2);
			MPIrequests.reserve(4);
			MPIstatuses.reserve(4);
		}


		///	<summary>
		///		Destructor.
		///	</summary>
		~BlobBoxesDegExchangeOP(){
			// MPI_Type_free(&MPI_LatonBox_double_type);
			// MPI_Type_free(&MPI_doubleArray);
			MPIrequests.clear();
			MPIstatuses.clear();
		}


		// ergonomic view 
		const std::vector<std::vector<LatLonBox<double>>>&
		ViewOriginalVecAllBlobBoxesDeg() const noexcept {
			return *this->in_;
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the BlobBoxesDegExchangeOP should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);

			//----------------------Send data first----------------------
			// Pack data into the send buffer  do NOT pre-clone edge columns for everyone.
			// We build the per-direction send buffer ONLY inside the
			// branch where this rank actually sends for that direction.
			// sendBlobBoxesDeg[0] = _vecAllBlobBoxesDeg[0];
			// sendBlobBoxesDeg[1] = _vecAllBlobBoxesDeg[_vecAllBlobBoxesDeg.size()-1];

			if (rank % 2 == 0) {
				if (rank == 0 || rank == size - 1) {
					exchangedvecAllBlobBoxesDeg.resize(in_->size() + 1);
				} else {
					exchangedvecAllBlobBoxesDeg.resize(in_->size() + 2);
				}
			} else {
				// odd ranks delay allocation; EndExchange() will just swap with *in_
			}




			// Send data
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						// Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						//Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}

				if (destRank > size - 1) {
					continue;
				}

				//----------------------Send sendBlobBoxesDeg----------------------
				// Only the odd number processors will send out data
				if (rank % 2 != 0) {
					MPI_Request request;

					// Send directly from the existing edge column, no staging clone.
					const std::vector<LatLonBox<double>>& src =
						(dir == DIR_LEFT) ? in_->front() : in_->back();
					
					const int nbytes = static_cast<int>(src.size() * sizeof(LatLonBox<double>));

					int rc = MPI_Isend(nbytes ? (void*)src.data() : nullptr,
									nbytes, MPI_BYTE,
									destRank, tag, m_comm, &request);
					if (rc != MPI_SUCCESS) {
						_EXCEPTION1("MPI_Isend failed (code %i)", rc);
					}
					MPIrequests.emplace_back(request);
					MPIstatuses.emplace_back();
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}

				if (sourceRank > size - 1) {
					continue;
				}

				// Only the even number processors will receive data 
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvBytes = 0;
					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count(&status, MPI_BYTE, &recvBytes);
					_ASSERT(recvBytes % sizeof(LatLonBox<double>) == 0);
					const int nElem = recvBytes / static_cast<int>(sizeof(LatLonBox<double>));

					recvBlobBoxesDeg[dir].resize(std::max(0, nElem));

					int rc = MPI_Irecv(
						recvBytes ? static_cast<void*>(recvBlobBoxesDeg[dir].data()) : nullptr,
						recvBytes, MPI_BYTE, sourceRank, tag, m_comm, &request);
					if (rc != MPI_SUCCESS) {
						_EXCEPTION1("MPI_Irecv failed (code %i)", rc);
					}
					MPIrequests.emplace_back(request);
					MPIstatuses.emplace_back();
				}
			}


		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the BlobBoxesDegExchangeOP data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for BOTH sends and receives (we tracked send requests above)
			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}


			MPIrequests.clear();
			MPIstatuses.clear();
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);



			// Build exchanged vector like Blobs/Tags:
			// - local columns are copied (cannot move from non-owning input)
			// - received edge columns are moved (cheap) into the ends.
			if ((rank % 2) == 0) {
				if (rank == 0) {
					// [0 .. N-1] from local input (swap to avoid copy)
					for (size_t i = 0; i < in_->size(); ++i) {
						exchangedvecAllBlobBoxesDeg[i].swap((*in_)[i]); 
					}
					// last column is the RIGHT edge received
					exchangedvecAllBlobBoxesDeg.back() = std::move(recvBlobBoxesDeg[DIR_RIGHT]);
				} else if (rank == size - 1) {
					// first column is the LEFT edge received
					exchangedvecAllBlobBoxesDeg.front() = std::move(recvBlobBoxesDeg[DIR_LEFT]);
					// [1 .. end] from local input (swap to avoid copy)
					for (size_t i = 1; i < exchangedvecAllBlobBoxesDeg.size(); ++i) {
						exchangedvecAllBlobBoxesDeg[i].swap((*in_)[i - 1]); 
					}
				} else {
					// interior: left recv, middle local, right recv
					exchangedvecAllBlobBoxesDeg.front() = std::move(recvBlobBoxesDeg[DIR_LEFT]);
					for (size_t i = 1; i + 1 < exchangedvecAllBlobBoxesDeg.size(); ++i) {
						exchangedvecAllBlobBoxesDeg[i].swap((*in_)[i - 1]);
					}
					exchangedvecAllBlobBoxesDeg.back()  = std::move(recvBlobBoxesDeg[DIR_RIGHT]);
				}
			} else {
				// Odd ranks: nothing changes; just hand back the original without copies.
				exchangedvecAllBlobBoxesDeg.swap(*in_); 
			}

			// We have moved from _vecAllBlobBoxesDeg where applicable; drop original & staging to reduce live memory
			for (auto& v : recvBlobBoxesDeg) v.shrink_to_fit();
			recvBlobBoxesDeg.clear();
			recvBlobBoxesDeg.shrink_to_fit();


		}

		///	<summary>
		///		Return the exchanged vecAllBlobTags
		///	</summary>
		std::vector<std::vector<LatLonBox<double>>>
		TakeExchangedVecAllBlobBoxesDeg() && noexcept {
			return std::move(this->exchangedvecAllBlobBoxesDeg);
		}
};

#endif 



#if defined(TEMPEST_MPIOMP) 

///	<summary>
///		A Class used for exchanging vecAllBlobs between processors
///	</summary>
class BlobsExchangeOp {
	private:

		///	<summary>
		///		Tag for MPI communication for sending/receiving blobs
		///	</summary>
		int blob_tag = 101;

		///	<summary>
		///		Tag for MPI communication for sending/receiving serialized index
		///	</summary>
		int indx_tag = 102;


		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;


		///	<summary>
		///		Serialize the vector<IndicatorSet> and generate the sendBlobsIndx array
		///	</summary>
		void Serialize(){

			sendBlobs.clear();sendBlobsIndx.clear();
			
			sendBlobs.resize(2);sendBlobsIndx.resize(2);

			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {

				// The vector<IndicatorSet> to that need to be serialized and sent
				// Avoid copying entire boundary column; bind by const reference
				const std::vector<IndicatorSet> & sendVecBlobs = (dir == DIR_LEFT) ? in_->front() : in_->back();


				// Roughly estimate total ints to reserve (exact sum keeps memory stable)
				size_t totalInts = 0;
				for (const auto & s : sendVecBlobs) totalInts += s.size();
				sendBlobs[dir].reserve(totalInts);
				sendBlobsIndx[dir].reserve(sendVecBlobs.size() + 1);


				int curIndx = 0;  //Point to the next empty slot for inserting a new set

				sendBlobsIndx[dir].push_back(curIndx);

				for (int i = 0; i < static_cast<int>(sendVecBlobs.size()); i++) {
					// Avoid copying each set; bind by const reference
					const IndicatorSet & curSet = sendVecBlobs[i];

					for (auto it = curSet.begin(); it != curSet.end(); it++) {
						sendBlobs[dir].push_back(*it);
						curIndx++;
					}
					sendBlobsIndx[dir].push_back(curIndx);//Now it records the starting position of the next IndicatorSet
				}
			}
		}

	protected:

		///	<summary>
		///		The initial local vecAllBlobs before the exchange, Non-owning pointer to the caller's original vecAllBlobs (no deep copy)
		///	</summary>
		std::vector<std::vector<IndicatorSet>>* in_ = nullptr;

		///	<summary>
		///		The vecAllBlobs that is after the exchange. it's column number is _vecAllAllBlobs' column number plus 2 (left and right) 
		///		except for p0 and pn-1 (for these two processors, the column number is _vecAllAllBlobs' column number plus 1).
		///	</summary>
		std::vector<std::vector<IndicatorSet>> exchangedVecAllBlobs;

		///	<summary>
		///		The buffer for vecBlobs that is serialized and will be sent
		///		sendBlobs[0] is the left vector and sendBlobs[1] is the right vector
		///	</summary>
		std::vector<std::vector<int>> sendBlobs;


		///	<summary>
		///		The Array recording the starting index for each set that is serialized
		///		sendBlobsIndx[0] is for the left vector and sendBlobsIndx[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> sendBlobsIndx;

		///	<summary>
		///		The buffer for the received serailized blobs array and need to be unserialized
		///		recvBlobs[0] is for the left vector and recvBlobs[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> recvBlobs;

		///	<summary>
		///		The Buffer  recording the starting index for each set that is in the recvBlobs
		///		recvBlobsIndx[0] is for the left vector and recvBlobsIndx[1] is for the right vector
		///	</summary>
		std::vector<std::vector<int>> recvBlobsIndx;

		///	<summary>
		///		The array recording the unserialize recvBlobs
		///		recvBlobsUnserial[0] is for the left vector and recvBlobsUnserial[1] is for the right vector	
		///	</summary>
		std::vector<std::vector<IndicatorSet>> recvBlobsUnserial;


	public:
		///	<summary>
		///		Construct the Operator with BlobsExchangeOp
		///		It will contruct the this->m_comm and this->_vecAllBlobs based on the input communicator and vecAllBlobs	
		///		Note that the input vecAllBlobs is passed by const reference and no deep copy is made.
		///	</summary>
		BlobsExchangeOp(MPI_Comm communicator, 
						std::vector< std::vector<IndicatorSet> > & vecAllBlobs) noexcept {
			this->in_ = &vecAllBlobs;  // non-owning
			this->m_comm = communicator;
		}

		///	<summary>
		///		Destructor for BlobsExchangeOp
		///	</summary>
		~BlobsExchangeOp(){
			MPIrequests.clear();
			MPIstatuses.clear();
			
		}


		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the BlobsExchangeOp should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {			
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);
			recvBlobs.resize(2);
			recvBlobsIndx.resize(2);

			//Merge the received vecBlobs into the new vecAllBlobs
			if (rank % 2 == 0) {
				if (rank == 0 || rank == size - 1) {
					exchangedVecAllBlobs.resize(in_->size() + 1);
				} else {
					exchangedVecAllBlobs.resize(in_->size() + 2);
				}
			} else {
				// odd ranks: delay allocation; EndExchange() will just swap with *in_
				// (leaving *in_ capacity unchanged). Keep empty here on purpose.
			}

			if ((rank % 2) != 0) {
				this->Serialize(); // builds both DIR_LEFT and DIR_RIGHT buffers
			}


			//----------------------Send sendBlobs/sendBlobsIndx data----------------------
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						//Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						// Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
					continue;

				}

				// only the odd number processors will send out data
				if (rank % 2 != 0) {

					
					//----------------------Send sendBlobs----------------------
					MPI_Request request;

					const int countA = static_cast<int>(sendBlobs[dir].size());
					MPI_Isend(countA ? sendBlobs[dir].data() : nullptr, countA, MPI_INT,
							destRank, blob_tag, m_comm, &request);
					MPIrequests.emplace_back(request);
					MPIstatuses.push_back(MPI_Status());


					//----------------------Send sendBlobsIndx----------------------
					MPI_Request indx_Request;

					const int countB = static_cast<int>(sendBlobsIndx[dir].size());
					MPI_Isend(countB ? sendBlobsIndx[dir].data() : nullptr, countB, MPI_INT,
							destRank, indx_tag, m_comm, &indx_Request);
					MPIrequests.emplace_back(indx_Request);
					MPIstatuses.push_back(MPI_Status());

				}



			}

			//----------------------Receive sendBlobs/recvBlobsIndx Data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;

				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}
				if (sourceRank > size - 1) {
					continue;

				}

				//----------------------Receive the serialized Blobs----------------------

				//only the even number processor will receive the blobs
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount = 0;

					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, blob_tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_INT, &recvCount );
					recvBlobs[dir].resize(recvCount);
					MPI_Irecv(recvCount ? recvBlobs[dir].data() : nullptr, recvCount, MPI_INT,
							sourceRank, blob_tag, m_comm, &request);
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());

					//----------------------Receive the index info for the Blobs----------------------
					MPI_Status indxStatus;
					MPI_Request indxRequest;
					int indxRecvCount = 0;

					// Use a non-blocking probe to know the incoming data size
					int indxFlag = 0;
					while(!indxFlag)
					{
						MPI_Iprobe( sourceRank, indx_tag, m_comm, &indxFlag, &indxStatus);
					}
					MPI_Get_count( &indxStatus, MPI_INT, &indxRecvCount);
					recvBlobsIndx[dir].resize(indxRecvCount);
					MPI_Irecv(indxRecvCount ? recvBlobsIndx[dir].data() : nullptr, indxRecvCount, MPI_INT,
							sourceRank, indx_tag, m_comm, &indxRequest);
					MPIrequests.emplace_back(std::move(indxRequest));
					MPIstatuses.push_back(MPI_Status());

				}

			}




		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the BlobsExchangeOp data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {

			// Wait for all Irecv to complete
			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}


			MPIrequests.clear();
			MPIstatuses.clear();

			sendBlobs.clear();        sendBlobs.shrink_to_fit();
			sendBlobsIndx.clear();    sendBlobsIndx.shrink_to_fit();

			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);








			// Only the even number processors need deserialize
			if (rank % 2 == 0) {

				// Avoid recvBlobsUnserial staging: directly deserialize into exchangedVecAllBlobs halos
				auto direct_deserialize_into = [&](int dir, std::vector<IndicatorSet> &target){
					target.clear();
					if (recvBlobsIndx[dir].empty()) {
						return;
					}
					const int nsets = static_cast<int>(recvBlobsIndx[dir].size()) - 1;
					if (nsets <= 0) {
						return;
					}
					target.reserve(nsets);
					for (int i = 0; i < nsets; ++i) {
						const int startIndx = recvBlobsIndx[dir][i];
						const int endIndx   = std::min(recvBlobsIndx[dir][i+1], int(recvBlobs[dir].size()));

						IndicatorSet curSet;
						for (int j = startIndx; j < endIndx; ++j){
							curSet.insert(recvBlobs[dir][j]);
						}
						target.push_back(std::move(curSet));
					}
				};
				// [CHANGED] end

				if (rank == 0) {
					for (int i = 0; i < static_cast<int>(exchangedVecAllBlobs.size()) - 1; i++) {
						exchangedVecAllBlobs[i].swap((*in_)[i]); 
					}
					// [CHANGED] build right halo directly into place
					direct_deserialize_into(DIR_RIGHT, exchangedVecAllBlobs.back());
				} else if (rank == size - 1) {
					// [CHANGED] build left halo directly into place
					direct_deserialize_into(DIR_LEFT, exchangedVecAllBlobs.front());
					for (int i = 1; i < static_cast<int>(exchangedVecAllBlobs.size()); i++) {
						exchangedVecAllBlobs[i].swap((*in_)[i - 1]);
					}
				} else {
					// [CHANGED] build both halos directly into place
					direct_deserialize_into(DIR_LEFT,  exchangedVecAllBlobs.front());
					for (int i = 1; i < static_cast<int>(exchangedVecAllBlobs.size()) - 1; i++) {
						exchangedVecAllBlobs[i].swap((*in_)[i - 1]);
					}
					direct_deserialize_into(DIR_RIGHT, exchangedVecAllBlobs.back());
				}

				// [CHANGED] aggressively drop incoming flat buffers once consumed
				for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
					std::vector<int>().swap(recvBlobs[dir]);
					std::vector<int>().swap(recvBlobsIndx[dir]);
				}

			} else {
				// For odd number processors, nothing is modified.
				exchangedVecAllBlobs.swap(*in_); 
			}


			recvBlobs.clear();            recvBlobs.shrink_to_fit();
			recvBlobsIndx.clear();        recvBlobsIndx.shrink_to_fit();



		}

		///	<summary>
		///		Return the exchanged VecAllBlobs,  move-only handoff to avoid copies
		///	</summary>
		std::vector<std::vector<IndicatorSet>> TakeExchangedVecAllBlobs() && noexcept {
			return std::move(this->exchangedVecAllBlobs);
		}

		///	<summary>
		///		View the unexchanged VecAllBlobs (non-owning), const view without copying
		///	</summary>
		const std::vector<std::vector<IndicatorSet>>& ViewOriginalVecAllBlobs() const noexcept {
			return *this->in_;
		}

		///	<summary>
		///		Return the exchanged VecAllBlobs, [DEPRECATED] Kept for compatibility; returns a COPY (expensive)
		///	</summary>
		std::vector<std::vector<IndicatorSet>> GetExchangedVecAllBlobs(){
			return this->exchangedVecAllBlobs;
		}

		///	<summary>
		///		Return the unexchanged VecAllBlobs, [DEPRECATED] Kept for compatibility; returns a COPY (expensive)
		///	</summary>
		std::vector<std::vector<IndicatorSet>> GetOriginalVecAllBlobs(){
			return *this->in_;
		}

};

///	<summary>
///		A Class used for exchanging vecGlobalTimes between processors
///		This class will update the local global time accordingly after the exchange process.
///	</summary>	
class GlobalTimesExchangeOp {
	private:
		///	<summary>
		///		Tag for MPI communication for vecGlobalTimes
		///	</summary>
		int tag = 104;

		///	<summary>
		///		The MPI Communicator
		///	</summary>
		MPI_Comm m_comm; 

		///	<summary>
		///		An array of MPI_Request.	
		///	</summary>
		std::vector<MPI_Request> MPIrequests;

		///	<summary>
		///		An array of MPI_Status.	
		///	</summary>
  		std::vector<MPI_Status> MPIstatuses;

	protected:
	///	<summary>
	///		The initial vecGlobalTimes before exchange (non-owning pointer to the caller's original vecGlobalTimes (no deep copy))
	///	</summary>
	 const std::vector<std::vector<Time>>* in_ = nullptr;

	///	<summary>
	///		The vecGlobalTimes that is after the exchange. exchangedVecGlobalTimes[pi]'s size will be  is _vecGlobalTimes[pi]'s size plus 2 (left and right) 
	///		except for p0 and pn-1 (for these two processors, exchangedVecGlobalTimes[pi]'s size is _vecGlobalTimes[pi]'s size plus 1).
	///	</summary>
	std::vector< std::vector<Time> > exchangedVecGlobalTimes;

	///	<summary>
	///		The buffer for vecGlobalTimes that will be sent
	///		sendTimes[0] is the left vector and sendTimes[1] is the right vector
	///	</summary>
	std::array<Time, 2> sendTimes{};

	///	<summary>
	///		The buffer for vecGlobalTimes that will be received
	///		recvTimes[0] is the left vector and sendTimes[1] is the right vector
	///	</summary>
	std::array<Time, 2> recvTimes{};

	///	<summary>
	///		Lowerbound of the file number that this processor is reading (starts with index 0 and inclusive)
	///	</summary>	
	int fileLowerBound = 0;

	///	<summary>
	///		Upperbound of the file number that this processor is reading (starts with index 0 and exclusive)
	///	</summary>
	int fileUpperBound = 0;

	public:

		///	<summary>
		///		It will contruct the this->m_comm this->_vecGlobalTimes this->fileLowerBound and this->fileUpperBound based on the input.
		///	</summary>
		GlobalTimesExchangeOp(
			MPI_Comm communicator,
			const std::vector<std::vector<Time>>& vecGlobalTimes,
			const int& processorResponsibalForFile_LB,
			const int& processorResponsibalForFile_UB
		) noexcept 
		{
			this->in_ = &vecGlobalTimes;                   
			this->m_comm = communicator;
			this->fileLowerBound = processorResponsibalForFile_LB;
			this->fileUpperBound = processorResponsibalForFile_UB;
		}

		///	<summary>
		///		Destructor
		///	</summary>
		~GlobalTimesExchangeOp(){
			MPIrequests.clear();
			MPIstatuses.clear();
		}

		///	<summary>
		///		Start the exchange process.
		/// 	this function is non-blocking and the data values in the GlobalTimesExchangeOp should not be modified
		/// 	The exchange values are not guaranteed to be current when this function returns and need to be used with the EndExchange()
		///	</summary>
		void StartExchange() {
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);


			//----------------------Send Data----------------------
			for (auto dir: {DIR_LEFT, DIR_RIGHT}) {
				int destRank;//Destination Rank
				if (dir == DIR_LEFT) {
					// Sending Data to the left
					if (rank == 0) {
						// Rank 0 Do Nothing
						continue;
					} else {
						destRank = rank - 1;
					}

				} else {
					// Sending  Data to the right
					if (rank == size - 1) {
						// Rank n-1 Do Nothing
						continue;
					} else {
						destRank = rank + 1;
					}
				}
				if (destRank > size - 1) {
					continue;
				}


				// Only the odd number processors will send out data
				if (rank % 2 != 0) {
					// Read from *in_ (non-owning), no staging vector copies.
					// Left boundary time = first time in our first local file
					sendTimes[0] = (*in_)[fileLowerBound].front();
					// Right boundary time = last time in our last local file
					sendTimes[1] = (*in_)[fileUpperBound - 1].back();

					MPI_Request request;
					// Send as raw bytes (size known): exactly one Time
					MPI_Isend(
						reinterpret_cast<void*>(&sendTimes[dir]),
						sizeof(Time),
						MPI_BYTE,
						destRank, tag, m_comm, &request
					);
					MPIrequests.emplace_back(request);
					MPIstatuses.emplace_back();
				}
			}

			//----------------------Then Receive data----------------------
			for (auto dir : {DIR_LEFT, DIR_RIGHT}) {
				int sourceRank;				
				if (dir == DIR_LEFT) {
					// Receive Data From the left.
					if (rank == 0) {// Rank 0 will not receive from the left
						continue;
					} else {
						sourceRank = rank - 1;
					}

				} else {
					// Receive Data From the right.
					if (rank == size - 1) {// rank n-1 will not receive from the right
						continue;

					} else {
						sourceRank = rank + 1;
					}

				}

				if (sourceRank > size - 1) {
					continue;
				}

				// Only the even number processors will receive data 
				if (rank % 2 == 0) {
					MPI_Status status;
					MPI_Request request;
					int recvCount;
					// Use a non-blocking probe to know the incoming data size
					int flag = 0;
					while(!flag)
					{
						MPI_Iprobe( sourceRank, tag, m_comm, &flag, &status );
					}
					MPI_Get_count( &status, MPI_BYTE, &recvCount );						
					MPI_Irecv(&recvTimes[dir], recvCount, MPI_BYTE,
							sourceRank, tag, m_comm, &request);
					MPIrequests.emplace_back(std::move(request));
					MPIstatuses.push_back(MPI_Status());
				}
			}




		}

		///	<summary>
		///		End the exchange process.
		// 		this function is blocking until:
		// 		- it is safe to modify the values in the GlobalTimesExchangeOp data without
		//   		affecting the exchange values for other processes
		// 		- the exchange values can be read: they contain to up-to-date values
		//   		from other processes
		///	</summary>
		void EndExchange() {
			// Wait for all Irecv to complete
			int result = MPI_Waitall( MPIrequests.size(), MPIrequests.data(), MPIstatuses.data());
			if (result != MPI_SUCCESS) {
				_EXCEPTION1("The MPI routine MPI_Waitall failed (code %i)", result);
			}


			MPIrequests.clear();
			MPIstatuses.clear();
			int rank, size;
			MPI_Comm_size(m_comm, &size);
			MPI_Comm_rank(m_comm, &rank);



			// [CHANGED] Build the result by moving from input to avoid deep copy,
			//           then mutate only the two touched columns in place.
			//           DO NOT swap with *in_, we return exchangedVecGlobalTimes.
			exchangedVecGlobalTimes = std::move(*in_);  // steal columns (O(outer) pointer moves)

			// Pack the data into the vecGlobalTimes 
			// Only the even number processors need to pack the data
			if ((rank % 2) == 0) {
				if (rank == 0) {
					// append right neighbor's first time to our last local file
					auto & lastCol = exchangedVecGlobalTimes[fileUpperBound - 1];
					lastCol.reserve(lastCol.size() + 1);                  // [CHANGED] pre-reserve
					lastCol.push_back(recvTimes[DIR_RIGHT]);
				} else if (rank == size - 1) {
					// prepend left neighbor's last time to our first local file
					auto & firstCol = exchangedVecGlobalTimes[fileLowerBound];
					firstCol.reserve(firstCol.size() + 1);                // [CHANGED] pre-reserve
					firstCol.insert(firstCol.begin(), recvTimes[DIR_LEFT]);
				} else {
					// interior: both sides
					auto & firstCol = exchangedVecGlobalTimes[fileLowerBound];
					auto & lastCol  = exchangedVecGlobalTimes[fileUpperBound - 1];

					firstCol.reserve(firstCol.size() + 1);                // [CHANGED] pre-reserve
					firstCol.insert(firstCol.begin(), recvTimes[DIR_LEFT]);

					lastCol.reserve(lastCol.size() + 1);                  // [CHANGED] pre-reserve
					lastCol.push_back(recvTimes[DIR_RIGHT]);
				}
			} else {
				// Odd ranks: no edits needed; exchangedVecGlobalTimes already holds original via move
			}



		}


		///	<summary>
		///		Return the exchanged vecGlobalTimes
		///	</summary>
		std::vector<std::vector<Time>> TakeExchangedVecGlobalTimes() && noexcept {
			return std::move(this->exchangedVecGlobalTimes);
		}
};
#endif 

#if defined(TEMPEST_MPIOMP) 


/// <summary>
///   A Class for MPI communication for sending/receiving Tag pairs (MapGraph)
/// </summary>
class MapGraphGatherOp {
private:
    /// <summary> Tag for MPI communication for sending/receiving Tag pairs </summary>
    int reduce_tag = 106;

    /// <summary> The MPI Communicator </summary>
    MPI_Comm m_comm;

    /// <summary> An array of MPI_Request. </summary>
    std::vector<MPI_Request> MPIrequests;

    /// <summary> An array of MPI_Status. </summary>
    std::vector<MPI_Status>  MPIstatuses;

    /// <summary>
    ///   Non-owning view of the input multimap (avoids ctor deep copy).
    ///   Falls back to _multimapTagGraph only if in_ == nullptr.
    /// </summary>
    const MapGraph* in_ = nullptr;

    /// <summary>
    ///   Serialize the MapGraph and generate the local std::pair<Tag,Tag> array.
    ///   Uses a single reserve() to avoid growth copies.
    /// </summary>
    void Serialize() {
        const MapGraph& src = (in_ ? *in_ : _multimapTagGraph);

        localPairs.clear();
        localPairs.reserve(src.size());
        for (const auto& kv : src) {
            // kv is pair<const Tag, Tag>; copy into pair<Tag, Tag>
            localPairs.emplace_back(kv.first, kv.second);
        }
    }

    /// <summary>
    ///   Deserialize local std::pair<Tag,Tag> array back into MapGraph
    ///   (only processor 0 will call it). Frees localPairs afterwards.
    /// </summary>
    void Deserialize() {
        outputTagGraph.clear();
        outputTagGraph.insert(localPairs.begin(), localPairs.end());

        // hygiene: drop the big buffer after building the graph
        std::vector<std::pair<Tag,Tag>>().swap(localPairs);
    }

protected:
    /// <summary>
    ///   Fallback owning storage (unused in normal path; kept for drop-in compatibility).
    /// </summary>
    MapGraph _multimapTagGraph;

    /// <summary> Serialized local multimapTagGraph. </summary>
    std::vector<std::pair<Tag, Tag>> localPairs;

    /// <summary> Output multimapTagGraph. (Only used for processor 0) </summary>
    MapGraph outputTagGraph;

public:
    /// <summary>
    ///   Construct the Operator with a non-owning view of multimapTagGraph.
    ///   Avoids deep copy; only stores the communicator and a pointer.
    /// </summary>
    MapGraphGatherOp(MPI_Comm communicator, const MapGraph& multimapTagGraph) {
        // Set communicator first so we can log with it
        this->m_comm = communicator;

        // Non-owning view (no deep copy)
        this->in_ = &multimapTagGraph;

    }

    /// <summary> Destructor </summary>
    ~MapGraphGatherOp() {
        MPIrequests.clear();
        MPIstatuses.clear();
    }

    /// <summary>
    ///   The MPI gather that will gather each local MapGraph to processor 0.
    ///   Hypercube pattern. Senders return immediately after sending.
    /// </summary>
    void Gather() {
        int rank=-1, size=1;
        MPI_Comm_rank(m_comm, &rank);
        MPI_Comm_size(m_comm, &size);

        Serialize();

        const int d = (size <= 1)
                      ? 0
                      : static_cast<int>(std::ceil(std::log2(static_cast<double>(size))));

        for (int j = 0; j < d; ++j) {
            const int bit = (1 << j);

            if ((rank & bit) != 0) {
                // send to partner then return
                const int destRank = rank ^ bit;
                if (destRank >= size) continue;

                const int nbytes = static_cast<int>(localPairs.size() * sizeof(std::pair<Tag,Tag>));
                MPI_Send(localPairs.data(), nbytes, MPI_BYTE, destRank, reduce_tag, m_comm);
                return;
            } else {
                const int sourceRank = rank ^ bit;
                if (sourceRank >= size) continue;

                // probe to get byte count
                MPI_Status status;
                int byteCount = 0;
                MPI_Probe(sourceRank, reduce_tag, m_comm, &status);
                MPI_Get_count(&status, MPI_BYTE, &byteCount);

                // resize once and receive directly into the tail (no temp buffer)
                const int recvPairs = byteCount / static_cast<int>(sizeof(std::pair<Tag,Tag>));
                const size_t off    = localPairs.size();
                localPairs.resize(off + recvPairs);

                MPI_Recv(localPairs.data() + off, byteCount, MPI_BYTE,
                         sourceRank, reduce_tag, m_comm, MPI_STATUS_IGNORE);
            }
        }
    }

    /// <summary>
    ///   Return the gathered multi-graph (only processor 0 should call).
    ///   Builds from localPairs and returns by move to avoid extra copy.
    /// </summary>
    MapGraph GetGatheredTagGraph() {
        int rank=-1;
        MPI_Comm_rank(m_comm, &rank);
        if (rank == 0) {
            Deserialize();
            return std::move(outputTagGraph);
        } else {
            _EXCEPTIONT("Only processor 0 should call GetGatheredTagGraph().");
        }
    }
};

} // namespace mpi
} // namespace stitchblobs

#endif // defined(TEMPEST_MPIOMP)

#endif // STITCHBLOBS_MPI_UTILITIES_H
