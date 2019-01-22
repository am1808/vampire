#ifndef SIM_INTERNAL_H_
#define SIM_INTERNAL_H_
//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

//---------------------------------------------------------------------
// Defines shared internal data structures and functions for the
// simulation methods implementation. These functions should
// not be accessed outside of the simulate module.
//---------------------------------------------------------------------

namespace sim{

   extern double applied_field_pulse_duration; // Squared field pulse

   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal shared variables used for the simulation
      //-----------------------------------------------------------------------------
      // simple initialised class for set variables
      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      struct mp_t{
         set_double_t slonczewski_aj; // adiabatic torque
         set_double_t slonczewski_bj; // field torque (non-adiabatic)
         set_double_t spin_hall_angle; // Spin Hall angle
         set_double_t prefactor_DL;    // prefactor for amplitude of damping-like SOT
         set_double_t prefactor_FL;    // prefactor for amplitude of field-like SOT
      };

      extern std::vector<sim::internal::mp_t> mp; // array of material properties
      extern std::vector<double> slonczewski_aj; // array of adiabatic spin torques
      extern std::vector<double> slonczewski_bj; // array of non-adiabatic spin torques
      extern std::vector<double> spin_hall_angle; // Array of spin Hall angles
      extern std::vector<double> prefactor_DL;   // Array of sot damping-like prefactors
      extern std::vector<double> prefactor_FL;   // Array of sot field-like prefactors
      extern std::vector<double> slonczewski_spin_polarization_unit_vector; // spin polarization direction
      extern std::vector<double> spin_orbit_torque_polarization_unit_vector; // SOT current polarization direction
      extern std::vector<double> spin_orbit_torque_H_applied_unit_vector;   // SOT applied field direction
      extern double spin_orbit_torque_polarization_magnitude; // magnitude of injected current for SOT
      extern double spin_orbit_torque_pulse_duration; // pulse duration of SOT
      extern double spin_orbit_torque_H_applied; //Applied field strenght for SOT

      extern int num_monte_carlo_preconditioning_steps;

      // internal function declarations
      extern void monte_carlo_preconditioning();

   } // end of internal namespace
} // end of sim namespace

#endif //SIM_INTERNAL_H_
