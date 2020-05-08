//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@nanohpc.com
/// @version 1.0
/// @date    10/06/2011
/// @internal
///	Created:		10/06/2011
///	Revision:	01/05/2020 by Andrea Meo, am1808@york.ac.uk
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

/// @brief Program to simulate Heat Assisted Magnetic Recording (HAMR)
///
/// @details Performs a time series with moving localised head field and temperature pulse
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    10/06/2011
///
/// @internal
///	Created:		10/06/2011
///	Revision:	01/05/2020 by Andrea Meo, am1808@york.ac.uk
///=====================================================================================
///
void hamr(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hamr has been called" << std::endl;}

		// Set equilibration temperature and field
		sim::temperature=sim::Teq;
		
		// Disable laser
		sim::head_laser_on=false;

		// Equilibrate system
		while(sim::time<sim::equilibration_time){
			
			sim::integrate(sim::partial_time);
			
			// Calculate magnetisation statistics
			stats::mag_m();
			
			// Output data
			vout::data();
		}

		// now enable laser
		sim::head_laser_on=true;
		
		/* ------------------------------------------------------------------  /
		/  Pulsed HAMR process                                                 /
		/  T(x,y,t) = (Tmax-Tmin)* exp( -(t-3tp)*(t-3tp)/(tp*tp) )             /
		/                        * exp( -(x-Xc) * (x-Xc)/(2*sigmax*sigmax) )   /
		/                        * exp( -(y-Yc) * (y-Yc)/(2*sigmay*sigmay) )   /
		/                                                                      /
		/  tp = pulse time; (Xc,Yc) = coordinates of bit centre                /
		/  sigmax,sigmay = FWHM in x,y-directions                              /
		/  ------------------------------------------------------------------ */
		if(sim::harm_continuous==false){
			std::cout << " >>>> Performing HAMR pulsed simulation" << std::endl;
			zlog << zTs() << " >>>> Performing HAMR pulsed simulation" << std::endl;
			// Calibrate head region to be not larger than system size
			if(sim::hamr_H_bounds_x > cs::system_dimensions[0]){
				 sim::hamr_H_bounds_x=cs::system_dimensions[0];
				 sim::hamr_bit_spacing_x = 0.0;
			}
			if(sim::hamr_H_bounds_y > cs::system_dimensions[1]){
				sim::hamr_H_bounds_y=cs::system_dimensions[1];
				sim::hamr_bit_spacing_y = 0.0;
			}

			// Initialise number of bits so that at least one bit is written
			int N_bit_x = 1;
			int N_bit_y = 1;

			// If simulation of single bit, set the head in the centre of the system
			if(sim::hamr_single_bit==true){
				std::cout << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
				zlog << zTs() << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
				sim::head_position[0] = cs::system_dimensions[0]*0.5;
				sim::head_position[1] = cs::system_dimensions[1]*0.5;
			}
			// Otherwise, determine the initial position of the head
			else{
				if( floor(cs::system_dimensions[0]/(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x + 1.0))>0 ){
					 N_bit_x = floor(cs::system_dimensions[0]/(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x + 1.0)); // number of bits in x
				}
				if( floor(cs::system_dimensions[1]/(sim::hamr_H_bounds_y+sim::hamr_bit_spacing_y + 1.0))>0 ){
					 N_bit_y = floor(cs::system_dimensions[1]/(sim::hamr_H_bounds_y+sim::hamr_bit_spacing_y + 1.0)); // number of bits in y
				}
				// Initialise head position ensuring that if only one bit, head is in the centre of the system
				if(N_bit_x > 1){
					 sim::head_position[0] = sim::hamr_H_bounds_x * 0.5 + sim::hamr_bit_spacing_x;
				}
				else{
					 sim::head_position[0] = cs::system_dimensions[0]*0.5;
				}
				if(N_bit_y > 1){
					 sim::head_position[1] = sim::hamr_H_bounds_y * 0.5 + sim::hamr_bit_spacing_y;
				}
				else{
					 sim::head_position[1] = cs::system_dimensions[1]*0.5;
				}
			}

			const uint64_t write_time = 6.0*ceil(sim::hamr_laser_peak_time/mp::dt_SI);
			const uint64_t ramp_time = int(sim::hamr_H_ramp_time/mp::dt_SI);
			const double Hmax = sim::Hmax; // max field
			const double Hmin = sim::Hmin;
			sim::H_applied = Hmin;

			uint64_t start_time=sim::time;

			/*std::cout << "Hmax,Hmin\t" << Hmax << "\t" << Hmin << std::endl;
			std::cout << "cooling time\t" << sim::hamr_laser_peak_time << std::endl;
			std::cout << cs::system_dimensions[0] << "\t" << cs::system_dimensions[1] << std::endl;
			std::cout << sim::hamr_H_bounds_x << "\t" << sim::hamr_H_bounds_y << std::endl;
			std::cout << sim::head_position[0] << "\t" << sim::head_position[1] << std::endl;
			std::cout <<N_bit_x << "\t" << N_bit_y << std::endl;
			std::cout << "tramp\t" << sim::hamr_H_ramp_time << "\t" << ramp_time << std::endl;
			std::cout << sim::total_time << "\t" << sim::partial_time << "\t" << write_time << std::endl;
			std::cout << " >>>> New total number of steps\t" << write_time*N_bit_x*N_bit_y << std::endl;*/
			std::cout << " >>>> Initial Head poistion:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;
			std::cout << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
			std::cout << " >>>> New total simulated time:\t" << write_time*N_bit_x*N_bit_y*mp::dt_SI << "\ts" << std::endl;
			zlog << zTs() << " >>>> Initial Head poistion:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;
			zlog << zTs() << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
			zlog << zTs() << " >>>> New total simulated time:\t" << write_time*N_bit_x*N_bit_y*mp::dt_SI << "\ts" << std::endl;

			// Perform HAMR simulation with static application of temperature and field pulses along x

			// Displace Head in y-direction
			for(int j_bit=0; j_bit<N_bit_y; j_bit++){
				// If more than one bit in x, set head starting from x-edge
				if(N_bit_x>1){
					sim::head_position[0] = sim::hamr_H_bounds_x * 0.5 + sim::hamr_bit_spacing_x;
				}
				// if only one bit in x, set the head in the centre of the x-system
				else{
					sim::head_position[0] = cs::system_dimensions[0]*0.5;
				}

				for(int i_bit=0; i_bit<N_bit_x; i_bit++){

					const int delta_write_time = (i_bit+j_bit*(N_bit_x));

//					std::cout << "\tNew HEAD position\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\t" << sim::time << "\t"
//							    << delta_write_time << "\t" << write_time*delta_write_time  << "\t" << write_time + write_time*delta_write_time <<  "\n" << std::endl;
					std::cout << " >>>> Moving HEAD **** New HEAD position\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;
					zlog << zTs() << " >>>> Moving HEAD **** New HEAD position\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;

					while(sim::time<write_time + write_time*delta_write_time){
						// loop over partial time
						for(int tt=0; tt < sim::partial_time; tt++){
							const uint64_t ramp_time_init = write_time*delta_write_time + ramp_time;
							const uint64_t ramp_time_end  = (write_time + write_time*delta_write_time) - ramp_time;
							// Determine max magnitude of external field to be applied in writing process
							if(sim::time <= ramp_time_init && abs(sim::H_applied) <= abs(Hmax)){
								sim::H_applied += (Hmax-Hmin) /sim::hamr_H_ramp_time * mp::dt_SI;
							}
							// In the centra region the max magnitude of the field is Hmax
							else if(sim::time > ramp_time_init && sim::time < ramp_time_end){
								 sim::H_applied = Hmax;
							}
							// Decrease field in the final region
							else if(sim::time >= ramp_time_end && sim::time <= write_time + write_time*delta_write_time){
								sim::H_applied -= (Hmax-Hmin) /sim::hamr_H_ramp_time * mp::dt_SI;
							}
							// Apply Gaussian temperature profile
							double centre_time = (0.5*write_time + write_time*delta_write_time)*mp::dt_SI;
							double time_from_centre=mp::dt_SI*double(sim::time-start_time)-centre_time;
							sim::temperature = sim::Tmin + (sim::Tmax-sim::Tmin)*exp(-(time_from_centre)*(time_from_centre)/(2.0*(sim::hamr_laser_peak_time)*(sim::hamr_laser_peak_time)));
							// Integrate system
							sim::integrate(1);
						} // end loop over time-step-increment=partial_time
						// Calculate magnetisation statistics
						stats::mag_m();
						// Output data
						vout::data();
					} // end of integration loop
					std::cout << " >>>> Finished writing bit\t(" << i_bit << "," << j_bit << ")" <<std::endl;
					zlog << zTs() << " >>>> Finished writing bit\t(" << i_bit << "," << j_bit << ")" <<std::endl;
					// Move head along x
					sim::head_position[0] += sim::hamr_H_bounds_x + sim::hamr_bit_spacing_x;
				}
				sim::head_position[1] += sim::hamr_H_bounds_y + sim::hamr_bit_spacing_y;
			} // End of displacement of Head along y
		}
		/* ------------------------------------------------------------------  /
		/  Continuous HAMR process                                             /
		/  T(x,y,t) = (Tmax-Tmin)* exp( -(x-v*t)*(x-v*t)/(2*sigmax*sigmax) )   /
		/                        * exp( -(y-Yc)*(y-Yc)  /(2*sigmay*sigmay) )   /
		/                                                                      /
		/  v=Head speed; Yc=Head y-coordinate                                  /
		/  sigmax,sigmay = FWHM in x,y-directions                              /
		/ ------------------------------------------------------------------- */
		else{
			std::cout << " >>>> Performing HAMR continuous simulation" << std::endl;
			zlog << zTs() << " >>>> Performing HAMR continuous simulation" << std::endl;
			// Calibrate head region to be not larger than system size
			if(sim::hamr_H_bounds_x > cs::system_dimensions[0]){
				sim::hamr_H_bounds_x = cs::system_dimensions[0];
				sim::hamr_bit_spacing_x = 0.0;
				std::cout << " >>>> Resizing sim:hamr-head-field-x to\t" << sim::hamr_H_bounds_x << "\tAng" << std::endl;
				zlog << zTs() << " >>>> Resizing sim:hamr-head-field-x to\t" << sim::hamr_H_bounds_x << "\tAng" << std::endl;
			}
			if(sim::hamr_H_bounds_y > cs::system_dimensions[1]){
				sim::hamr_H_bounds_y = cs::system_dimensions[1];
				sim::hamr_bit_spacing_y = 0.0;
				std::cout << " >>>> Resizing sim:hamr-head-field-y to\t" << sim::hamr_H_bounds_y << "\tAng" << std::endl;
				zlog << zTs() << " >>>> Resizing sim:hamr-head-field-y to\t" << sim::hamr_H_bounds_y << "\tAng" << std::endl;
			}

			// Initialise number of bits so that at least one bit is written
			int N_bit_x = 1;
			int N_bit_y = 1;
			int x_bit = 0;
			int y_bit = 0;
			int bit = 0;
			// Initial position is with head away from film
			sim::head_position[0] = 0.0 - sim::hamr_H_bounds_x*0.5;
			sim::hamr_bit_spacing_x = 0.0;
			sim::head_position[1] = sim::hamr_H_bounds_y*0.5 + sim::hamr_bit_spacing_y;

			// If simulation of single bit, set the head in the centre of the y-system dimension
			if(sim::hamr_single_bit==true){
				std::cout << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
				zlog << zTs() << " >>>> Performing HAMR simulation in Single Bit mode" << std::endl;
				sim::hamr_H_osc_amplit = sim::hamr_H_bounds_x;
				std::cout << " >>>> Setting sim:hamr-field-oscillation-frequency to\t" << sim::hamr_H_osc_amplit << "\tAng" << std::endl;
				zlog << zTs() << " >>>> Setting sim:hamr-field-oscillation-frequency to\t" << sim::hamr_H_osc_amplit << "\tAng" << std::endl;
//				sim::head_position[1] = cs::system_dimensions[1]*0.5;
			}
			// Otherwise, determine the total number of bits in x and y-directions
			else{
				// Determine number of bits
				if( floor(cs::system_dimensions[0]/(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x + 1.0))>0 ){
					N_bit_x = floor(cs::system_dimensions[0]/(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x + 1.0)); // number of bits in x
				}
				if( floor(cs::system_dimensions[1]/(sim::hamr_H_bounds_y+sim::hamr_bit_spacing_y + 1.0))>0 ){
					N_bit_y = floor(cs::system_dimensions[1]/(sim::hamr_H_bounds_y+sim::hamr_bit_spacing_y + 1.0)); // number of bits in y
				}
			}

			// Set times
			const uint64_t peak_time = int( (((sim::hamr_H_bounds_x + sim::hamr_bit_spacing_x*0.5)/(sim::head_speed))/6.0) /mp::dt_SI );
			const uint64_t write_time = 6*peak_time;
			const uint64_t pre_write_time = 0.5*write_time;
			const uint64_t ramp_time = int(sim::hamr_H_ramp_time/mp::dt_SI);
			std::cout << " head_speed\t" << sim::head_speed*1e-10 << "\tm/s" << std::endl;
			std::cout << " peak_time\t"  << peak_time*mp::dt_SI   << "\ts  " << std::endl;
			std::cout << " write_time\t" << write_time*mp::dt_SI  << "\ts  " << std::endl;
			std::cout << " ramp_time\t"  << ramp_time*mp::dt_SI   << "\ts  " << std::endl;
			std::cout << " equl time\t"  << sim::equilibration_time*mp::dt_SI << "\ts" << std::endl;
			std::cout << " total time\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time )*mp::dt_SI << "\ts" << std::endl;
			std::cout << " total time + equl time\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;

			const double Hmax = sim::Hmax; // max field
			const double Hmin = sim::Hmin;
			sim::H_applied = Hmin;

			std::cout << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
			std::cout << " >>>> Initial Head position:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;
			std::cout << " >>>> Head velocity:\t" << sim::head_speed*1e-10 << "\tm/s" << std::endl;
			std::cout << " >>>> Time per bit:\t" << (write_time*N_bit_x*N_bit_y)*mp::dt_SI << "\ts" << std::endl;
			std::cout << " >>>> New total simulated time:\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;
			zlog << zTs() << " >>>> Number of bits in x,y:\t" << N_bit_x << "\t" << N_bit_y << std::endl;
			zlog << zTs() << " >>>> Initial Head poistion:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;
			zlog << zTs() << " >>>> Head velocity:\t" << sim::head_speed*1e-10 << "\tm/s" << std::endl;
			zlog << zTs() << " >>>> Time per bit:\t" << (write_time*N_bit_x*N_bit_y)*mp::dt_SI << "\ts" << std::endl;
			zlog << zTs() << " >>>> New total simulated time:\t" << ( (write_time*N_bit_x*N_bit_y) + 2*pre_write_time + sim::equilibration_time )*mp::dt_SI << "\ts" << std::endl;

//			while(sim::time-sim::equilibration_time<write_time*N_bit_x*N_bit_y && y_bit<N_bit_y && sim::head_position[1]<cs::system_dimensions[1]){
			while(sim::time-sim::equilibration_time<write_time*N_bit_x*N_bit_y+pre_write_time*2 &&
					sim::head_position[1]*y_bit<cs::system_dimensions[1]){  /*&& bit<N_bit_x*N_bit_y+1 && y_bit<N_bit_y*/
				// Update head position in y-direction
				sim::head_position[1] = sim::hamr_H_bounds_y*(0.5+y_bit) + sim::hamr_bit_spacing_y;

//				while( sim::head_position[0]-0.5*sim::hamr_H_bounds_x < cs::system_dimensions[0] &&
//						 sim::head_position[0]-0.5*sim::hamr_H_bounds_x < N_bit_x*(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x) &&
//						 bit<(N_bit_x)*N_bit_y+1 ){
				while( sim::head_position[0] < cs::system_dimensions[0]+0.5*sim::hamr_H_bounds_x ){ // &&
//						 sim::head_position[0]  -0.5*sim::hamr_H_bounds_x   < N_bit_x*(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x) ){

					std::cout << "\n >>>> Moving HEAD **** New HEAD position:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng\n" << std::endl;
					zlog << zTs() << " >>>> Moving HEAD **** New HEAD position:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << "\tAng" << std::endl;

//					while( sim::head_position[0] < cs::system_dimensions[0]+0.5*sim::hamr_H_bounds_x &&
//							 sim::head_position[0] < sim::hamr_H_bounds_x*(bit - N_bit_x*(y_bit)+1)    &&
//							 x_bit < N_bit_x ){
					while( sim::head_position[0] < sim::hamr_H_bounds_x*(x_bit+1) &&
							 sim::head_position[0] < cs::system_dimensions[0]+0.5*sim::hamr_H_bounds_x){
						// Update head position in x-direction
						sim::head_position[0] = ( (sim::head_speed) * (sim::time-sim::equilibration_time)*mp::dt_SI - 0.5*sim::hamr_H_bounds_x ) - (N_bit_x*sim::hamr_H_bounds_x)*y_bit;

						// loop over partial time
						for(int tt=0; tt < sim::partial_time; tt++){
							const uint64_t field_time = (sim::time-sim::equilibration_time) - write_time*bit;
							const uint64_t ramp_time_end  = write_time - ramp_time;
							// Set applied field to minimum in pre-write and post-write zones
							if(field_time < pre_write_time || sim::head_position[0]>cs::system_dimensions[0]){
								sim::H_applied = Hmin;
							}
							// Determine max magnitude of external field during initial ramp
							else if(field_time <= ramp_time+pre_write_time && abs(sim::H_applied) <= abs(Hmax)){
								sim::H_applied += (Hmax-Hmin) /sim::hamr_H_ramp_time * mp::dt_SI;
							}
							// In the central region the max magnitude of the field is Hmax
							else if(field_time > ramp_time+pre_write_time && field_time < ramp_time_end+pre_write_time){
								sim::H_applied = Hmax;
							}
							// Decrease field in the final region
							else if(field_time >= ramp_time_end+pre_write_time && field_time <= write_time+pre_write_time){
								sim::H_applied -= (Hmax-Hmin) /sim::hamr_H_ramp_time * mp::dt_SI;
							}
							// Set temperature to be peak temperature.
							sim::temperature = sim::Tmax;

							// Integrate system
							sim::integrate(1);
						} // end loop over time-step-increment=partial_time
						// Calculate magnetisation statistics
						stats::mag_m();
						// Output data
						vout::data();
					}  // End of one single bit
					std::cout << "\n >>>> Finished writing bit\t(" << x_bit << "," << y_bit << ")\t Hz\t"
							    << sim::H_vec[2]*((-1.0)*double(2*(int(sim::head_position[0]/sim::hamr_H_osc_amplit)%2)-1)) << "\n" << std::endl;
					zlog << zTs() << " >>>> Finished writing bit\t(" << x_bit << "," << y_bit << ")\t Hz\t"
							        << sim::H_vec[2]*((-1.0)*double(2*(int(sim::head_position[0]/sim::hamr_H_osc_amplit)%2)-1)) << std::endl;
					++x_bit;
					++bit;
				}  // End of continuous displacement along x
				std::cout << " >>>> Final track head position:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << std::endl;
				std::cout << " \t     N_bit_x (sim::hamr_H_bounds_x + sim::hamr_bit_spacing_x)\t" << N_bit_x*(sim::hamr_H_bounds_x+sim::hamr_bit_spacing_x) << std::endl;
				std::cout << " \t     x_bit\t" << x_bit << std::endl;
				std::cout << " \t     y_bit\t" << y_bit << std::endl;
				std::cout << " \t     bit\t" << bit << std::endl;
				std::cout << " \t     N_bit_x\t" << N_bit_x << std::endl;
				std::cout << " \t     N_bit_y\t" << N_bit_y << std::endl;
				std::cout << "\n >>>> Reset head at the beginning of downtrack and move along offtrack direction\n" << std::endl;
				zlog << zTs() << " >>>> Final track head position:\t" << sim::head_position[0] << "\t" << sim::head_position[1] << std::endl;
				zlog << zTs() << " >>>> Reset head at the beginning of downtrack and move along offtrack direction" << std::endl;
				// Update head position to allow correct evaluation of statement in while() loops
				sim::head_position[0] = 0.0 - 0.5*sim::hamr_H_bounds_x;
				sim::head_position[1] += sim::hamr_H_bounds_y + sim::hamr_bit_spacing_y;
				// Reset counter for bit in down track direction
				x_bit = 0;
				++y_bit;
			}  // End of loop over N_bit_y

//			// force outputting last point of simulation
//			if(sim::time-sim::equilibration_time>=write_time*N_bit_x*N_bit_y){
//				// Disable laser
//				sim::head_laser_on=false;
//				std::cout << "\n>>>> Disable laser and integrate system for 1 time-step" << std::endl;
//				// Integrate
//				sim::integrate(1);
//				// Calculate magnetisation statistics
//				stats::mag_m();
//				// Output data
//				vout::data();
//				std::cout << "\n>>>> Outputting system at the end of HAMR continuous simulations\n" << std::endl;
//				// Reactivate laser
//				sim::head_laser_on=true;
//			}
		}  // End of continuous HAMR process
	} // end of hamr()

}//end of namespace program
