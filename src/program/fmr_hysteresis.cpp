//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2015 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains the Time Series program
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Function to calculate hystereis loop when applying an AC field
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
///
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///
void fmr_hysteresis(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::hysteresis-loop has been called" << std::endl;

	//double temp=sim::temperature;

	// disable fmr fields for equilibration
	sim::enable_fmr = false;

/*------------------- Beginnign equilibration region added from hysteresis.cpp -----------------*/
	// Setup min and max fields and increment (uT)
	int iHmax=vmath::iround(double(sim::Hmax)*1.0E6);
	int miHmax=-iHmax;
	int parity_old;
	int iH_old;
	int start_time;

	// Equilibrate system in saturation field, i.e. the largest between equilibration and maximum field set by the user
   if(sim::Heq >= sim::Hmax){
	   sim::H_applied=sim::Heq;
   }
   else{
   	sim::H_applied=sim::Hmax;
   }

	// Initialise sim::integrate only if it not a checkpoint
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else sim::integrate(sim::equilibration_time);

   // Hinc must be positive
	int iHinc=vmath::iround(double(fabs(sim::Hinc))*1.0E6);

   int Hfield;
   int iparity=sim::parity;
	parity_old=iparity;

   // Save value of iH from previous simulation
	if(sim::load_checkpoint_continue_flag) iH_old=int(sim::iH);
/*------------------- End equilibration region added from hysteresis.cpp -----------------*/

	//sim::temperature=temp;

   // Reset mean magnetisation counters
   stats::mag_m_reset();

	// enable fmr fields
	sim::enable_fmr = true;

	// Initialize direction along z if not already set
	if(sim::fmr_field_unit_vector.size() != 3){
		sim::fmr_field_unit_vector.resize(3,0.0);
		sim::fmr_field_unit_vector[2] = 1.0;
	}

/*--------------- Hystersis loop simulation (from hysteresis.cpp) -------------*/
	// Perform Field Loop -parity
	while(iparity<2){
		// If checkpoint is loaded with continue flag, then set up correctly max,min field values
		if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag)
		{
			//necessary to upload value of iH_old when loading the checkpoint !!!
			iH_old=int(sim::iH);
			//Setup min and max fields and increment (uT)
			if(parity_old<0){
				if(iparity<0) miHmax=iH_old;
				else if(iparity>0 && iH_old<=0) miHmax=iH_old; //miHmax=(iHmax-iHinc);
				else if(iparity>0 && iH_old>0) miHmax=-(iHmax);
			}
			else if(parity_old>0) miHmax=iH_old;
			Hfield=miHmax;
		}
		else	Hfield=miHmax;

		// Perform Field Loop -field
		while(Hfield<=iHmax){

			// Set applied field (Tesla)
			sim::H_applied=double(Hfield)*double(iparity)*1.0e-6;

			// Reset start time
			start_time=sim::time;

			// Reset mean magnetisation counters
			stats::mag_m_reset();

			// Integrate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);

				// Calculate mag_m, mag
				stats::mag_m();

			} // End of integration loop

			// Increment of iH
			Hfield+=iHinc;
			sim::iH=int64_t(Hfield); //sim::iH+=iHinc;

			// Output to screen and file after each field
			vout::data();

		} // End of field loop

		// Increment of parity
		iparity+=2;
		sim::parity=int64_t(iparity);

	} // End of parity loop
/*--------------- End of Hystersis loop simulation (from hysteresis.cpp) -------------*/

}

}//end of namespace program
