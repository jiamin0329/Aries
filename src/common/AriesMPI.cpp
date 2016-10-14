/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Simple utility class for interfacing with MPI
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    07-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */

#include "AriesMPI.hpp"

/* Aries includes */
#include "Logger.hpp"
#include "Utilities.hpp"
/* C++ includes */
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

using namespace std;

#ifdef __INSURE__
/*
 *  These are defined in mpich mpi.h and break the insure compile.
 *  This may impact Globus in some way, at least from the comments
 *  in the mpi.h header file.  Why mpich externs something that is
 *  not defined in the mpich is confusing and probably just broken.
 */
int MPICHX_TOPOLOGY_DEPTHS;
int MPICHX_TOPOLOGY_COLORS;
int MPICHX_PARALLELSOCKETS_PARAMETERS;
#endif

/*
 *================================================================================
 *    Class namespaces
 *================================================================================
 */
namespace ARIES
{
    bool AriesMPI::d_mpiIsInitialized = false;
    bool AriesMPI::d_weStartedMpi(false);
    AriesMPI AriesMPI::d_ariesWorld(MPI_COMM_NULL);

    bool AriesMPI::d_callAbortInSerialInsteadOfExit = false;
    bool AriesMPI::d_callAbortInParallelInsteadOfMpiAbort = false;
    int AriesMPI::d_invalidRank = -1;
    
    AriesMPI::AriesMPI(const Comm& comm):
            d_comm(comm),
            d_rank(-1),
            d_size(-1)
    {
        if (comm != MPI_COMM_NULL)
        {
#ifdef ARIES_HAVE_MPI
            if (d_mpiIsInitialized)
            {
                MPI_Comm_rank(d_comm, &d_rank);
                MPI_Comm_size(d_comm, &d_size);
            }
#else
            d_rank = 0;
            d_size = 1;
#endif
        }
    }
    
    AriesMPI::AriesMPI(const AriesMPI& other):
            d_comm(other.d_comm),
            d_rank(other.d_rank),
            d_size(other.d_size)
    {
    }
    
    void AriesMPI::SetCommunicator(const Comm& comm)
    {
        d_comm = comm;
        d_rank = 0;
        d_size = 1;
#ifdef ARIES_HAVE_MPI
        if (d_mpiIsInitialized)
        {
            if (d_comm != MPI_COMM_NULL)
            {
                MPI_Comm_rank(d_comm, &d_rank);
                MPI_Comm_size(d_comm, &d_size);
            }
        }
#endif
    }

    void AriesMPI::DupCommunicator(const AriesMPI& r)
    {
#ifdef ARIES_HAVE_MPI
        int rval = r.Comm_dup(&d_comm);
        if (rval != MPI_SUCCESS)
        {
            ARIES_ERROR("AriesMPI::dupCommunicator: Error duplicating\n" << "communicator.");
        }
        MPI_Comm_rank(d_comm, &d_rank);
        MPI_Comm_size(d_comm, &d_size);
        ARIES_ASSERT(d_rank == r.d_rank);
        ARIES_ASSERT(d_size == r.d_size);
#else
        d_comm = r.d_comm;
        d_rank = r.d_rank;
        d_size = r.d_size;
#endif
    }

    void AriesMPI::FreeCommunicator()
    {
#ifdef ARIES_HAVE_MPI
        if (d_comm != MPI_COMM_NULL)
        {
            ARIES_ASSERT(AriesMPI::usingMPI());
            Comm_free(&d_comm);
            // d_comm is now set to MPI_COMM_NULL;
        }
#else
        d_comm = MPI_COMM_NULL;
#endif
        d_rank = 0;
        d_size = 1;
    }
    
    int AriesMPI::CompareCommunicator(const AriesMPI& r) const
    {
#ifdef ARIES_HAVE_MPI
        int compare_result;
        int mpi_err = Comm_compare(d_comm, r.d_comm, &compare_result);
        if (mpi_err != MPI_SUCCESS)
        {
            ARIES_ERROR("AriesMPI::compareCommunicator: Error comparing two communicators.");
        }
        return compare_result;
#else
        NULL_USE(r);
        return d_comm == r.d_comm ? MPI_IDENT : MPI_CONGRUENT;
#endif
    }

    // Static MPI wrappers matching MPI interfaces.
    int AriesMPI::Comm_rank(Comm comm, int* rank)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(comm);
#endif
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Comm_rank is a no-op without run-time MPI!");
        }
        int rval = MPI_SUCCESS;
#ifdef ARIES_HAVE_MPI
        rval = MPI_Comm_rank(comm, rank);
#else
        *rank = 0;
#endif
        return rval;
    }

    int AriesMPI::Comm_size(Comm comm, int* size)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(comm);
#endif
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Comm_size is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        return MPI_Comm_size(comm, size);

#else
        *size = 1;
        return MPI_SUCCESS;
#endif
    }
    
    int AriesMPI::Comm_compare(Comm comm1, Comm comm2, int* result)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(comm1);
        NULL_USE(comm2);
        NULL_USE(result);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Comm_compare is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Comm_compare(comm1, comm2, result);
        }
#endif
        return rval;
    }

    int AriesMPI::Comm_free(Comm* comm)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(comm);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Comm_free is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Comm_free(comm);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Finalized(int* flag)
    {
        int rval = MPI_SUCCESS;
#ifdef ARIES_HAVE_MPI
        rval = MPI_Finalized(flag);
#else
        *flag = true;
#endif
        return rval;
    }

    int AriesMPI::Get_count(Status* status, Datatype datatype, int* count)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(status);
        NULL_USE(datatype);
        NULL_USE(count);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Get_count is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Get_count(status, datatype, count);
        }
#endif
        return rval;
    }

    int AriesMPI::Request_free(Request* request)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(request);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Get_count is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Request_free(request);
        }
#endif
        return rval;
    }

    int AriesMPI::Test(Request* request, int* flag, Status* status)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(request);
        NULL_USE(flag);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Test is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Test(request, flag, status);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Test_cancelled(Status* status, int* flag)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(status);
        NULL_USE(flag);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Test_canceled is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Test_cancelled(status, flag);
        }
#endif
        return rval;
    }

    int AriesMPI::Wait(Request* request, Status* status)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(request);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Wait is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Wait(request, status);
        }
#endif
        return rval;
    }

    int AriesMPI::Waitall(int count, Request* reqs, Status* stats)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(count);
        NULL_USE(reqs);
        NULL_USE(stats);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Waitall is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Waitall(count, reqs, stats);
        }
#endif
        return rval;
    }

    int AriesMPI::Waitany(int count, Request* array_of_requests, int* index, Status* status)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(count);
        NULL_USE(array_of_requests);
        NULL_USE(index);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Waitany is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Waitany(count, array_of_requests, index, status);
        }
#endif
        return rval;
    }

    int AriesMPI::Waitsome(int incount, Request* array_of_requests, int* outcount, int* array_of_indices, Status* array_of_statuses)
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(incount);
        NULL_USE(array_of_requests);
        NULL_USE(outcount);
        NULL_USE(array_of_indices);
        NULL_USE(array_of_statuses);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Waitsome is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
        }
#endif
        return rval;
    }
    
    double AriesMPI::Wtime()
    {
        double rval = 0.0;
        if (!d_mpiIsInitialized)
        {
            struct tms tmp_tms;
            clock_t clock_ticks_since_reference = times(&tmp_tms);
            const double clock_ticks_per_second = double(sysconf(_SC_CLK_TCK));
            rval = double(clock_ticks_since_reference) / clock_ticks_per_second;
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Wtime();
        }
#endif
        return rval;
    }
    //* Static MPI wrappers matching MPI interfaces.

    // MPI wrappers matching MPI interfaces.
    int AriesMPI::Allgather(void* sendbuf, int sendcount, Datatype sendtype, void* recvbuf, int recvcount, Datatype recvtype) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(sendcount);
        NULL_USE(sendtype);
        NULL_USE(recvbuf);
        NULL_USE(recvcount);
        NULL_USE(recvtype);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Allgather is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Allgatherv(void* sendbuf, int sendcounts, Datatype sendtype,
                             void* recvbuf, int* recvcounts, int* displs, Datatype recvtype) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(sendcounts);
        NULL_USE(sendtype);
        NULL_USE(recvbuf);
        NULL_USE(recvcounts);
        NULL_USE(displs);
        NULL_USE(recvtype);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Algatherv is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Allgatherv(sendbuf,
                                  sendcounts,
                                  sendtype,
                                  recvbuf,
                                  recvcounts,
                                  displs,
                                  recvtype,
                                  d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Allreduce(void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(recvbuf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(op);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Allreduce is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, d_comm);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Attr_get(int keyval, void* attribute_val, int* flag) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(keyval);
        NULL_USE(attribute_val);
        NULL_USE(flag);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Attr_get is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Attr_get(d_comm, keyval, attribute_val, flag);
        }
#endif
        return rval;
    }

    int AriesMPI::Barrier() const
    {
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            // do nothing
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Barrier(d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Bcast(void* buffer, int count, Datatype datatype, int root) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(buffer);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(root);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            // do nothing
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Bcast(buffer, count, datatype, root, d_comm);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Comm_dup(Comm* newcomm) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(newcomm);
#endif
        *newcomm = MPI_COMM_NULL;
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Comm_dup is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Comm_dup(d_comm, newcomm);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Comm_rank(int* rank) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(rank);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            *rank = 0;
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Comm_rank(d_comm, rank);
        }
#endif
        return rval;
    }

    int AriesMPI::Comm_size(int* size) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(size);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            *size = 1;
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Comm_size(d_comm, size);
        }
#endif
        return rval;
    }

    int AriesMPI::Gather(void* sendbuf, int sendcount, Datatype sendtype,
                         void* recvbuf, int recvcount, Datatype recvtype, int root) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(sendcount);
        NULL_USE(sendtype);
        NULL_USE(recvbuf);
        NULL_USE(recvcount);
        NULL_USE(recvtype);
        NULL_USE(root);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Gather is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Gatherv(void* sendbuf, int sendcount, Datatype sendtype,
                          void* recvbuf, int* recvcounts, int* displs, Datatype recvtype, int root) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(sendcount);
        NULL_USE(sendtype);
        NULL_USE(recvbuf);
        NULL_USE(recvcounts);
        NULL_USE(displs);
        NULL_USE(recvtype);
        NULL_USE(root);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Gatherv is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Gatherv(sendbuf,
                               sendcount,
                               sendtype,
                               recvbuf,
                               recvcounts,
                               displs,
                               recvtype,
                               root,
                               d_comm);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Iprobe(int source, int tag, int* flag, Status* status) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(source);
        NULL_USE(tag);
        NULL_USE(flag);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Iprobe is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Iprobe(source, tag, d_comm, flag, status);
        }
#endif
        return rval;
    }

    int AriesMPI::Isend(void* buf, int count, Datatype datatype,
                        int dest, int tag, Request* req) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(buf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(dest);
        NULL_USE(tag);
        NULL_USE(req);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Isend is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Isend(buf, count, datatype, dest, tag, d_comm, req);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Irecv(void* buf, int count, Datatype datatype,
                        int source, int tag, Request* request) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(buf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(source);
        NULL_USE(tag);
        NULL_USE(request);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Irecv is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Irecv(buf, count, datatype, source, tag, d_comm, request);
        }
#endif
        return rval;
    }

    int AriesMPI::Probe(int source, int tag, Status* status) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(source);
        NULL_USE(tag);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Probe is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Probe(source, tag, d_comm, status);
        }
#endif
        return rval;
    }

    int AriesMPI::Recv(void* buf, int count, Datatype datatype, int source, int tag, Status* status) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(buf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(source);
        NULL_USE(tag);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Recv is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Recv(buf, count, datatype, source, tag, d_comm, status);
        }
#endif
        return rval;
    }

    int AriesMPI::Reduce(void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op, int root) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(recvbuf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(op);
        NULL_USE(root);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Reduce is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Send(void* buf, int count, Datatype datatype, int dest, int tag) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(buf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(dest);
        NULL_USE(tag);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Send is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Send(buf, count, datatype, dest, tag, d_comm);
        }
#endif
        return rval;
    }

    int AriesMPI::Sendrecv(void* sendbuf, int sendcount, Datatype sendtype, int dest, int sendtag,
                           void* recvbuf, int recvcount, Datatype recvtype, int source, int recvtag, Status* status) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(sendcount);
        NULL_USE(sendtype);
        NULL_USE(dest);
        NULL_USE(sendtag);
        NULL_USE(recvbuf);
        NULL_USE(recvcount);
        NULL_USE(recvtype);
        NULL_USE(source);
        NULL_USE(recvtag);
        NULL_USE(status);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Send is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,
                                recvbuf, recvcount, recvtype, source, recvtag, d_comm, status);
        }
#endif
        return rval;
    }
    
    int AriesMPI::Scan(void* sendbuf, void* recvbuf, int count, Datatype datatype, Op op) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(sendbuf);
        NULL_USE(recvbuf);
        NULL_USE(count);
        NULL_USE(datatype);
        NULL_USE(op);
#endif
        int rval = MPI_SUCCESS;
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::Scan is a no-op without run-time MPI!");
        }
#ifdef ARIES_HAVE_MPI
        else
        {
            rval = MPI_Scan(sendbuf, recvbuf, count, datatype, op, d_comm);
        }
#endif
        return rval;
    }

    /*
     * Methods named like MPI's native interfaces (without the MPI_ prefix)
     * are wrappers for the native interfaces.  The AriesMPI versions
     * introduce a flag to determine whether MPI is really used at run time.
     * When the run-time flag is on, these wrappers are identical to the MPI
     * versions.  When the flag is off, most of these methods are no-ops
     * (which is not necessarily the same as calling the MPI functions with
     * only 1 process in the communicator).
     */
    int AriesMPI::AllReduce(int* x, int count, Op op, int* ranks_of_extrema) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(x);
        NULL_USE(count);
        NULL_USE(op);
        NULL_USE(ranks_of_extrema);
#endif
        if ((op == MPI_MINLOC || op == MPI_MAXLOC) && ranks_of_extrema == 0)
        {
            ARIES_ERROR("AriesMPI::AllReduce: If you specify reduce\n"
                        << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
                        << "provide space for the ranks in the 'ranks_of_extrema'\n"
                        << "argument.");
        }
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::AllReduce is a no-op without run-time MPI!");
        }

        int rval = MPI_SUCCESS;
        /*
         * Get ranks of extrema if user operation specified it or user
         * specified min/max operation and provides space for rank.
         */
        bool get_ranks_of_extrema =
                op == MPI_MINLOC ? true :
                op == MPI_MAXLOC ? true :
                ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

        if (!get_ranks_of_extrema)
        {
            std::vector<int> recv_buf(count);
            rval = Allreduce(x, &recv_buf[0], count, MPI_INT, op);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c];
            }
        }
        else
        {
            Op locop =
                    op == MPI_MIN ? MPI_MINLOC :
                    op == MPI_MAX ? MPI_MAXLOC :
                    op;
            IntIntStruct* send_buf = new IntIntStruct[count];
            IntIntStruct* recv_buf = new IntIntStruct[count];
            for (int c = 0; c < count; ++c)
            {
                send_buf[c].j = x[c];
                send_buf[c].i = d_rank;
            }
            rval = Allreduce(send_buf, recv_buf, count, MPI_2INT, locop);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c].j;
                ranks_of_extrema[c] = recv_buf[c].i;
            }

            delete[] send_buf;
            delete[] recv_buf;
        }

        return rval;
    }

    int AriesMPI::AllReduce(double* x, int count, Op op, int* ranks_of_extrema) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(x);
        NULL_USE(count);
        NULL_USE(op);
        NULL_USE(ranks_of_extrema);
#endif
        if ((op == MPI_MINLOC || op == MPI_MAXLOC) &&
            ranks_of_extrema == 0)
        {
            ARIES_ERROR("AriesMPI::AllReduce: If you specify reduce\n"
                        << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
                        << "provide space for the ranks in the 'ranks_of_extrema'\n"
                        << "argument.");
        }
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::AllReduce is a no-op without run-time MPI!");
        }

        int rval = MPI_SUCCESS;
        /*
         * Get ranks of extrema if user operation specified it or user
         * specified min/max operation and provides space for rank.
         */
        bool get_ranks_of_extrema =
                op == MPI_MINLOC ? true :
                op == MPI_MAXLOC ? true :
                ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

        if (!get_ranks_of_extrema)
        {
            std::vector<double> recv_buf(count);
            rval = Allreduce(x, &recv_buf[0], count, MPI_DOUBLE, op);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c];
            }
        }
        else
        {
            Op locop =
                    op == MPI_MIN ? MPI_MINLOC :
                    op == MPI_MAX ? MPI_MAXLOC :
                    op;
            DoubleIntStruct* send_buf = new DoubleIntStruct[count];
            DoubleIntStruct* recv_buf = new DoubleIntStruct[count];
            for (int c = 0; c < count; ++c)
            {
                send_buf[c].d = x[c];
                send_buf[c].i = d_rank;
            }
            rval = Allreduce(send_buf, recv_buf, count, MPI_DOUBLE_INT, locop);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c].d;
                ranks_of_extrema[c] = recv_buf[c].i;
            }

            delete[] send_buf;
            delete[] recv_buf;
        }

        return rval;
    }
    
    int AriesMPI::AllReduce( float* x, int count, Op op, int* ranks_of_extrema) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(x);
        NULL_USE(count);
        NULL_USE(op);
        NULL_USE(ranks_of_extrema);
#endif
        if ((op == MPI_MINLOC || op == MPI_MAXLOC) &&
            ranks_of_extrema == 0)
        {
            ARIES_ERROR("AriesMPI::AllReduce: If you specify reduce\n"
                        << "operation MPI_MINLOC or MPI_MAXLOC, you must\n"
                        << "provide space for the ranks in the 'ranks_of_extrema'\n"
                        << "argument.");
        }
        if (!d_mpiIsInitialized)
        {
            ARIES_ERROR("AriesMPI::AllReduce is a no-op without run-time MPI!");
        }

        int rval = MPI_SUCCESS;
        /*
         * Get ranks of extrema if user operation specified it or user
         * specified min/max operation and provides space for rank.
         */
        bool get_ranks_of_extrema =
                op == MPI_MINLOC ? true :
                op == MPI_MAXLOC ? true :
                ranks_of_extrema != 0 && (op == MPI_MIN || op == MPI_MAX);

        if (!get_ranks_of_extrema)
        {
            std::vector<float> recv_buf(count);
            rval = Allreduce(x, &recv_buf[0], count, MPI_FLOAT, op);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c];
            }
        }
        else
        {
            Op locop =
                    op == MPI_MIN ? MPI_MINLOC :
                    op == MPI_MAX ? MPI_MAXLOC :
                    op;
            FloatIntStruct* send_buf = new FloatIntStruct[count];
            FloatIntStruct* recv_buf = new FloatIntStruct[count];
            for (int c = 0; c < count; ++c)
            {
                send_buf[c].f = x[c];
                send_buf[c].i = d_rank;
            }
            rval = Allreduce(send_buf, recv_buf, count, MPI_FLOAT_INT, locop);
            for (int c = 0; c < count; ++c)
            {
                x[c] = recv_buf[c].f;
                ranks_of_extrema[c] = recv_buf[c].i;
            }

            delete[] send_buf;
            delete[] recv_buf;
        }

        return rval;
    }

    int AriesMPI::ParallelPrefixSum(int* x, int count, int tag) const
    {
#ifndef ARIES_HAVE_MPI
        NULL_USE(x);
        NULL_USE(count);
        NULL_USE(tag);
#endif
        // Scratch data.
        std::vector<int> send_scr(count), recv_scr(count);

        Request send_req, recv_req;
        Status send_stat, recv_stat;
        int mpi_err = MPI_SUCCESS;

        for (int distance = 1; distance < d_size; distance *= 2)
        {
            const int recv_from = d_rank - distance;
            const int send_to = d_rank + distance;

            if (recv_from >= 0)
            {
                mpi_err = Irecv(&recv_scr[0], count, MPI_INT, recv_from, tag, &recv_req);
                if (mpi_err != MPI_SUCCESS)
                {
                    return mpi_err;
                }
            }

            if (send_to < d_size)
            {
                send_scr.clear();
                send_scr.insert(send_scr.end(), x, x + count);
                mpi_err = Isend(&send_scr[0], count, MPI_INT, send_to, tag, &send_req);
                if (mpi_err != MPI_SUCCESS)
                {
                    return mpi_err;
                }
            }

            if (recv_from >= 0)
            {
                mpi_err = Wait(&recv_req, &recv_stat);
                if (mpi_err != MPI_SUCCESS)
                {
                    return mpi_err;
                }
                for (int i = 0; i < count; ++i)
                {
                    x[i] += recv_scr[i];
                }
            }

            if (send_to < d_size)
            {
                mpi_err = Wait(&send_req, &send_stat);
                if (mpi_err != MPI_SUCCESS)
                {
                    return mpi_err;
                }
            }

        }

        return MPI_SUCCESS;
    }

    bool AriesMPI::HasReceivableMessage(Status* status, int source, int tag) const
    {
        int flag = false;
        if (d_mpiIsInitialized)
        {
            AriesMPI::Status tmp_status;
            Barrier();
            int mpi_err = Iprobe(source, tag, &flag, status ? status : &tmp_status);
            if (mpi_err != MPI_SUCCESS)
            {
                ARIES_ERROR("AriesMPI::hasReceivableMessage: Error probing for message." << std::endl);
            }
            Barrier();
        }
        return flag == true;
    }
    //* MPI wrappers matching MPI interfaces.
    
    void AriesMPI::Abort()
    {
#ifdef ARIES_HAVE_MPI
        const AriesMPI& mpi(AriesMPI::GetAriesWorld());
        if (mpi.getSize() > 1)
        {
            if (d_callAbortInParallelInsteadOfMpiAbort)
            {
                ::abort();
            }
            else
            {
                MPI_Abort(mpi.getCommunicator(), -1);
            }
        }
        else
        {
            if (d_callAbortInSerialInsteadOfExit)
            {
                ::abort();
            }
            else
            {
                exit(-1);
            }
        }
#else
        if (d_callAbortInSerialInsteadOfExit)
        {
            ::abort();
        }
        else
        {
            exit(-1);
        }
#endif
    }

    void AriesMPI::Init(int* argc, char** argv[])
    {
#ifdef ARIES_HAVE_MPI
        MPI_Init(argc, argv);
        d_mpiIsInitialized = true;
        d_weStartedMpi = true;

        Comm dup_comm;
        MPI_Comm_dup(MPI_COMM_WORLD, &dup_comm);
        d_ariesWorld.setCommunicator(dup_comm);
#else
        NULL_USE(argc);
        NULL_USE(argv);
        d_ariesWorld.d_comm = MPI_COMM_WORLD;
        d_ariesWorld.d_size = 1;
        d_ariesWorld.d_rank = 0;
#endif

        if (getenv("ARIES_ABORT_ON_ERROR"))
        {
            AriesMPI::SetCallAbortInSerialInsteadOfExit(true);
            AriesMPI::SetCallAbortInParallelInsteadOfMPIAbort(true);
        }
    }
    
    void AriesMPI::Init(Comm comm)
    {
        if (comm == MPI_COMM_NULL)
        {
            std::cerr << "AriesMPI::Init: invalid initializing Communicator."
                      << std::endl;
        }
#ifdef ARIES_HAVE_MPI
        d_mpiIsInitialized = true;
        d_weStartedMpi = false;

        Comm dup_comm;
        MPI_Comm_dup(comm, &dup_comm);
        d_ariesWorld.setCommunicator(dup_comm);
#endif
        if (getenv("ARIES_ABORT_ON_ERROR"))
        {
            AriesMPI::SetCallAbortInSerialInsteadOfExit(true);
            AriesMPI::SetCallAbortInParallelInsteadOfMPIAbort(true);
        }
    }
    
    void AriesMPI::InitMPIDisabled()
    {
        d_mpiIsInitialized = false;
        d_weStartedMpi = false;

        d_ariesWorld.d_comm = MPI_COMM_WORLD;
        d_ariesWorld.d_size = 1;
        d_ariesWorld.d_rank = 0;

        if (getenv("ARIES_ABORT_ON_ERROR"))
        {
            AriesMPI::SetCallAbortInSerialInsteadOfExit(true);
            AriesMPI::SetCallAbortInParallelInsteadOfMPIAbort(true);
        }
    }

    void AriesMPI::Finalize()
    {
#ifdef ARIES_HAVE_MPI
        if (d_mpiIsInitialized)
        {
            MPI_Comm_free(&d_ariesWorld.d_comm);
        }
        else
        {
            d_ariesWorld.d_comm = MPI_COMM_NULL;
        }

        if (d_weStartedMpi)
        {
            MPI_Finalize();
        }
#endif
    }
    
#ifndef ARIES_HAVE_MPI
    AriesMPI::Status::Status():
            MPI_SOURCE(-1),
            MPI_TAG(-1),
            MPI_ERROR(-1)
    {
    }
#endif

} // end namespace ARIES

