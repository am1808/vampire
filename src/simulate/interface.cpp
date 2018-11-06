//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      // set maximum allowable value for time steps (10^12)
      const uint64_t max_time = 1000000000000;
      const std::string max_time_str = "1,000,000,000,000";

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="slonczewski-spin-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::slonczewski_spin_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::spin_orbit_torque_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-polarization-magnitude";
      if(word==test){
         double je = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(je, word, line, prefix, unit, "", -1.0e18, 1.0e18,"input","-1e18 - 1e18");
         sim::internal::spin_orbit_torque_polarization_magnitude = je;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-pulse-duration";
      if(word==test){
         double dt = atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-3,"input","0.01 attosecond - 1 microsecond");
         sim::internal::spin_orbit_torque_pulse_duration = dt;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-applied-field-strength";
      if(word==test){
         double H = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(H, word, line, prefix, unit, "field", -1.e4, 1.0e4,"input","+/- 10,000 T");
         sim::internal::spin_orbit_torque_H_applied=H;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-applied-field-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::spin_orbit_torque_H_applied_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="applied-field-pulse-duration";
      if(word==test){
         double dt = atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-3,"input","0.01 attosecond - 1 microsecond");
         sim::applied_field_pulse_duration = dt;
         return true;
      }
      //-------------------------------------------------------------------
      test="preconditioning-steps";
      if(word==test){
         int n = atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(n, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         sim::internal::num_monte_carlo_preconditioning_steps = n;
         return true;
      }
      //-------------------------------------------------------------------
      test="time-step";
      if(word==test){
         double dt = atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
         mp::dt_SI = dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="total-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::total_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="loop-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::loop_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::equilibration_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="time-steps-increment";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 1, max_time,"input","1 - "+max_time_str);
         sim::partial_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > sim::internal::mp.size() && super_index + 1 < 101) sim::internal::mp.resize(super_index + 1);

      //------------------------------------------------------------
      std::string test="slonczewski-adiabatic-spin-torque";
      /*
         aj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_aj.set(aj);
         return true;
      }
      //------------------------------------------------------------
      test="slonczewski-non-adiabatic-spin-torque";
      /*
         bj parameter for material in slonczewski torque calculation
         */
      if(word==test){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", 0.0, 1.0e2,"input","0 - 100T");
         sim::internal::mp[super_index].slonczewski_bj.set(bj);
         return true;
      }
      //------------------------------------------------------------
      test="spin-hall-angle";
      /*
         Spin Hall angle
         */
      if(word==test){
         double th_SH=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(th_SH, word, line, prefix, unit, "", -1.0e18, 1.0e18,"input","-1e18 - 1e18");
         sim::internal::mp[super_index].spin_hall_angle.set(th_SH);
         return true;
      }
      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace sim
