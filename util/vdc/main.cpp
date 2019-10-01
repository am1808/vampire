//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>

// program header
#include "vdc.hpp"

int main(int argc, char* argv[]){

   // process command line arguments
   //command(argc, argv);
   long int min_file_id = 0;        //index of initial file to process
   long int max_file_id = 99999999; //index of final file to process
   // If only one value is specified, then only the minimum is changed
   if(argc>0 && argc<2){
      min_file_id = atoi(argv[1]);
   }
   // If two values are given, both maximum and minimum file indexes are updated
   else if(argc>0 && argc<3){
      min_file_id = atoi(argv[1]);
      max_file_id = atoi(argv[2])+1;
   }

   if(vdc::verbose){
      std::cout << "|------------------------------------------------------------|" << std::endl;
      std::cout << "|              Vampire Data Converter for v5+                |" << std::endl;
      std::cout << "|------------------------------------------------------------|" << std::endl;
   }

   // process coordinates
   vdc::process_coordinates();

   // process spin files
   vdc::process_spins(min_file_id,max_file_id);

   return 0;

}
