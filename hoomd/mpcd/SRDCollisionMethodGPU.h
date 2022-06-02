// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/SRDCollisionMethodGPU.h
 * \brief Declaration of mpcd::SRDCollisionMethodGPU
 */

#ifndef MPCD_SRD_COLLISION_METHOD_GPU_H_
#define MPCD_SRD_COLLISION_METHOD_GPU_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "SRDCollisionMethod.h"
#include "hoomd/Autotuner.h"

namespace hoomd
    {
namespace mpcd
    {
class PYBIND11_EXPORT SRDCollisionMethodGPU : public mpcd::SRDCollisionMethod
    {
    public:
    //! Constructor
    SRDCollisionMethodGPU(std::shared_ptr<mpcd::SystemData> sysdata,
                          unsigned int cur_timestep,
                          unsigned int period,
                          int phase,
                          uint16_t seed,
                          std::shared_ptr<mpcd::CellThermoCompute> thermo);

    protected:
    //! Randomly draw cell rotation vectors
    virtual void drawRotationVectors(uint64_t timestep);

    //! Apply rotation matrix to velocities
    virtual void rotate(uint64_t timestep);

    private:
    std::unique_ptr<Autotuner> m_tuner_rotvec; //!< Tuner for drawing rotation vectors
    std::unique_ptr<Autotuner> m_tuner_rotate; //!< Tuner for rotating velocities
    };

namespace detail
    {
//! Export SRDCollisionMethodGPU to python
void export_SRDCollisionMethodGPU(pybind11::module& m);
    } // end namespace detail

    }  // end namespace mpcd
    }  // end namespace hoomd
#endif // MPCD_SRD_COLLISION_METHOD_GPU_H_
