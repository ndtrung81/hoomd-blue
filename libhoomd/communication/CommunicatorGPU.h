/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2008-2011 Ames Laboratory
Iowa State University and The Regents of the University of Michigan All rights
reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of HOOMD-blue, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* All publications and presentations based on HOOMD-blue, including any reports
or published results obtained, in whole or in part, with HOOMD-blue, will
acknowledge its use according to the terms posted at the time of submission on:
http://codeblue.umich.edu/hoomd-blue/citations.html

* Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
http://codeblue.umich.edu/hoomd-blue/

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of HOOMD-blue's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Maintainer: jglaser

/*! \file CommunicatorGPU.h
    \brief Defines the CommunicatorGPU class
*/

#ifndef __COMMUNICATOR_GPU_H__
#define __COMMUNICATOR_GPU_H__

#ifdef ENABLE_MPI
#ifdef ENABLE_CUDA

#include "Communicator.h"

#include "CommunicatorGPU.cuh"

#include "GPUFlags.h"
#include "GPUArray.h"

/*! \ingroup communication
*/

//! Class that handles MPI communication (GPU version)
/*! CommunicatorGPU is the GPU implementation of the base communication class.
*/
class CommunicatorGPU : public Communicator
    {
    public:
        //! Constructor
        /*! \param sysdef system definition the communicator is associated with
         *  \param decomposition Information about the decomposition of the global simulation domain
         */
        CommunicatorGPU(boost::shared_ptr<SystemDefinition> sysdef,
                        boost::shared_ptr<DomainDecomposition> decomposition);
        virtual ~CommunicatorGPU();

        //! \name communication methods
        //@{

        /*! Perform ghosts update
         *
         * \param timestep The time step
         */
        virtual void beginUpdateGhosts(unsigned int timestep);

        /*! Finish ghost update
         *
         * \param timestep The time step
         */
        virtual void finishUpdateGhosts(unsigned int timestep);

        //! Transfer particles between neighboring domains
        virtual void migrateParticles();

        //! Build a ghost particle list, exchange ghost particle data with neighboring processors
        virtual void exchangeGhosts();

        //@}

        //! Set maximum number of communication stages
        /*! \param max_stages Maximum number of communication stages
         */
        void setMaxStages(unsigned int max_stages)
            {
            m_max_stages = max_stages;
            initializeCommunicationStages();
            forceMigrate();
            }

    protected:
        //! Helper class to perform the communication tasks related to bonded groups
        template<class group_data>
        class GroupCommunicatorGPU
            {
            public:
                typedef struct rank_element<typename group_data::ranks_t> rank_element_t;
                typedef typename group_data::packed_t group_element_t;

                //! Constructor
                GroupCommunicatorGPU(CommunicatorGPU& gpu_comm, boost::shared_ptr<group_data> gdata);

                //! Migrate groups
                /*! \param incomplete If true, mark all groups that have non-local members and update local
                 *         member rank information. Otherwise, mark only groups flagged for communication
                 *         in particle data
                 *
                 * A group is marked for sending by setting its rtag to GROUP_NOT_LOCAL, and by updating
                 * the rank information with the destination ranks (or the local ranks if incomplete=true)
                 */
                void migrateGroups(bool incomplete);

                //! Mark ghost particles
                /* All particles that need to be sent as ghosts because they are members
                 * of incomplete groups are marked, and destination ranks are compute accordingly.
                 *
                 * \param plans Array of particle plans to write to
                 * \param mask Mask for allowed sending directions
                 */
                void markGhostParticles(const GPUArray<unsigned int>& plans, unsigned int mask);

            private:
                CommunicatorGPU& m_gpu_comm;                            //!< The outer class
                boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //< The execution configuration
                boost::shared_ptr<group_data> m_gdata;                  //!< The group data

                GPUVector<unsigned int> m_rank_mask;                    //!< Bitfield for every group to keep track of updated rank fields
                GPUVector<unsigned int> m_scan;                         //!< Temporary array for exclusive scan of group membership information

                GPUVector<rank_element_t> m_ranks_out;                  //!< Packed ranks data
                GPUVector<rank_element_t> m_ranks_sendbuf;              //!< Send buffer for ranks information
                GPUVector<rank_element_t> m_ranks_recvbuf;              //!< Recv buffer for ranks information

                GPUVector<group_element_t> m_groups_out;                //!< Packed group data
                GPUVector<unsigned int> m_rank_mask_out;                //!< Output buffer for rank update bitfields
                GPUVector<group_element_t> m_groups_sendbuf;            //!< Send buffer for groups
                GPUVector<group_element_t> m_groups_recvbuf;            //!< Recv buffer for groups
                GPUVector<group_element_t> m_groups_in;                 //!< Input buffer of unique groups
            };

    private:
        /* General communication */
        unsigned int m_max_stages;                     //!< Maximum number of (dependent) communication stages
        unsigned int m_num_stages;                     //!< Number of stages
        std::vector<unsigned int> m_comm_mask;         //!< Communication mask per stage
        std::vector<int> m_stages;                     //!< Communication stage per unique neighbor

        /* Particle migration */
        GPUVector<pdata_element> m_gpu_sendbuf;        //!< Send buffer for particle data
        GPUVector<pdata_element> m_gpu_recvbuf;        //!< Receive buffer for particle data
        GPUVector<unsigned int> m_comm_flags;          //!< Output buffer for communication flags

        GPUVector<unsigned int> m_send_keys;           //!< Destination rank for particles

        /* Communication of bonded groups */
        GroupCommunicatorGPU<BondData> m_bond_comm;    //!< Communication helper for bonds
        friend class GroupCommunicatorGPU<BondData>;

        GroupCommunicatorGPU<AngleData> m_angle_comm;  //!< Communication helper for angles
        friend class GroupCommunicatorGPU<AngleData>;

        GroupCommunicatorGPU<DihedralData> m_dihedral_comm;  //!< Communication helper for dihedrals
        friend class GroupCommunicatorGPU<DihedralData>;

        GroupCommunicatorGPU<ImproperData> m_improper_comm;  //!< Communication helper for impropers
        friend class GroupCommunicatorGPU<ImproperData>;

        /* Ghost communication */
        GPUVector<unsigned int> m_tag_ghost_sendbuf;   //!< List of ghost particles tags per stage, ordered by neighbor
        GPUVector<unsigned int> m_tag_ghost_recvbuf;   //!< Buffer for recveiving particle tags
        GPUVector<Scalar4> m_pos_ghost_sendbuf;        //<! Buffer for sending ghost positions
        GPUVector<Scalar4> m_pos_ghost_recvbuf;        //<! Buffer for receiving ghost positions

        GPUVector<Scalar4> m_vel_ghost_sendbuf;        //<! Buffer for sending ghost velocities
        GPUVector<Scalar4> m_vel_ghost_recvbuf;        //<! Buffer for receiving ghost velocities

        GPUVector<Scalar> m_charge_ghost_sendbuf;      //!< Buffer for sending ghost charges
        GPUVector<Scalar> m_charge_ghost_recvbuf;      //!< Buffer for sending ghost charges

        GPUVector<Scalar> m_diameter_ghost_sendbuf;    //!< Buffer for sending ghost charges
        GPUVector<Scalar> m_diameter_ghost_recvbuf;    //!< Buffer for sending ghost charges

        GPUVector<Scalar4> m_orientation_ghost_sendbuf;//<! Buffer for sending ghost orientations
        GPUVector<Scalar4> m_orientation_ghost_recvbuf;//<! Buffer for receiving ghost orientations

        GPUVector<unsigned int> m_ghost_begin;          //!< Begin index for every stage and neighbor in send buf
        GPUVector<unsigned int> m_ghost_end;            //!< Begin index for every and neighbor in send buf

        GPUVector<uint2> m_ghost_idx_adj;             //!< Indices and adjacency relationships of ghosts to send
        GPUVector<unsigned int> m_ghost_neigh;        //!< Neighbor ranks for every ghost particle
        GPUVector<unsigned int> m_ghost_plan;         //!< Plans for every particle
        std::vector<unsigned int> m_idx_offs;         //!< Per-stage offset into ghost idx list

        GPUVector<unsigned int> m_neigh_counts;       //!< List of number of neighbors to send ghost to (temp array)

        std::vector<std::vector<unsigned int> > m_n_send_ghosts; //!< Number of ghosts to send per stage and neighbor
        std::vector<std::vector<unsigned int> > m_n_recv_ghosts; //!< Number of ghosts to receive per stage and neighbor
        std::vector<std::vector<unsigned int> > m_ghost_recv_offs;//!< Begin offset in recv buf per stage and neighbor
        std::vector<std::vector<unsigned int> > m_ghost_send_offs;//!< Begin offset in send buf per stage and neighbor
        std::vector<std::vector<unsigned int> > m_n_max_send_ghosts; //!< Maximum number of ghosts to send per stage and neighbor
        std::vector<std::vector<unsigned int> > m_n_max_recv_ghosts; //!< Maximum number of ghosts to receive per stage and neighbor
        std::vector<GPUVector<unsigned int> > m_ghost_recv_idx; //!< Indices of ghosts in recv buffer
        std::vector<GPUVector<unsigned int> > m_ghost_send_idx; //!< Indices of ghosts in send buffer
 
        float m_ghost_resize_factor;                   //!< Factor to use for amortized buffer resizing
 
        std::vector<unsigned int> m_n_send_ghosts_tot; //!< Total number of sent ghosts per stage
        std::vector<unsigned int> m_n_recv_ghosts_tot; //!< Total number of received ghosts per stage

        mgpu::ContextPtr m_mgpu_context;              //!< MGPU context
        cudaEvent_t m_event;                          //!< CUDA event for synchronization

        //! Helper function to allocate various buffers
        void allocateBuffers();

        //! Helper function to set up communication stages
        void initializeCommunicationStages();
    };

//! Export CommunicatorGPU class to python
void export_CommunicatorGPU();

#endif // ENABLE_CUDA
#endif // ENABLE_MPI
#endif // __COMMUNICATOR_GPU_H
