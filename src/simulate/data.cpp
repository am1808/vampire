//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "sim.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------
   uint64_t time         = 0; // time step counter
   uint64_t total_time   = 10000; // total time steps (non-loop code)
   uint64_t loop_time    = 10000; // loop time steps (hysteresis/temperature loops)
   uint64_t partial_time = 1000; // same as time-step-increment
   uint64_t equilibration_time = 0; // equilibration time steps

   namespace internal{

      //----------------------------------------------------------------------------
      // Shared variables used within sim module
      //---------------------------------------------------------------------------
      std::vector<sim::internal::mp_t> mp; // array of material properties
      std::vector<double> slonczewski_aj; // array of adiabatic spin torques
      std::vector<double> slonczewski_bj; // array of non-adiabatic spin torques
      std::vector<double> slonczewski_spin_polarization_unit_vector(3,0.0); // spin polarization direction

      int num_monte_carlo_preconditioning_steps = 0;

   } // end of internal namespace

} // end of sim namespace
