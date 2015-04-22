//**********************************************************************************//
//    Copyright (C) 2015  Ovidio Pena <ovidio@bytesfall.com>                        //
//    Copyright (C) 2015  V. R. Iglesias <viglesias@gmail.com>                      //
//                                                                                  //
//    This file is part of scattecc                                                 //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify          //
//    it under the terms of the GNU General Public License as published by          //
//    the Free Software Foundation, either version 3 of the License, or             //
//    (at your option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful,               //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//    GNU General Public License for more details.                                  //
//                                                                                  //
//    The only additional remark is that we expect that all publications            //
//    describing work using this software, or all commercial products               //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "eccmie.h"

const double PI=3.14159265358979323846;

//***********************************************************************************//
// This is the main function of 'scattecc', here we read the parameters as           //
// arguments passed to the program which should be executed with the following       //
// syntaxis:                                                                         //
// ./scattecc -l xi mi.r mi.i xh mh.r mh.i -s dx -i aoi [-t ti tf nt] [-c comment]   //
//                                                                                   //
// When all the parameters were correctly passed we setup xi (xh) and mi (mh),       //
// containing the size parameters and refractive indexes of the inclusion (host,     //
// respectively. We also set dx (shift of the inclusion) and aoi (angle of incidence //
// and then call the function eccMie.                                                //
// If the calculation is successful the results are printed with the following       //
// format:                                                                           //
//                                                                                   //
//    * If no comment was passed:                                                    //
//        'Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                                    //
//                                                                                   //
//    * If a comment was passed:                                                     //
//        'comment, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo'                           //
//***********************************************************************************//
int main(int argc, char *argv[]) {
  try {
    std::vector<std::string> args;
    args.assign(argv, argv + argc);
    std::string error_msg(std::string("Insufficient parameters.\nUsage: ") + args[0]
                          + " -l xi mi.r mi.i xh mh.r mh.i "
                          + " -s dx -i aoi "
                          + "[-t ti tf nt] [-c comment]\n");
    enum mode_states {nothing, read_xi, read_xh, read_mi_r, read_mi_i, read_mh_r, read_mh_i, read_dx, read_aoi, read_ti, read_tf, read_nt, read_comment};
    // for (auto arg : args) std::cout<< arg <<std::endl;
    std::string comment;
    int has_comment = 0;
    int i;
    double xi, xh, dx, aoi;
    std::vector<double> Theta;
    std::complex<double> mi, mh;
    std::vector<std::complex<double> > S1, S2;
    double Qext, Qabs, Qsca, Qbk, Qpr, g, Albedo;

    double ti = 0.0, tf = 90.0;
    int nt = 0;    
    if (argc < 5) throw std::invalid_argument(error_msg);

    int mode = nothing; 
    double tmp_mr;
    for (auto arg : args) {
      // For each arg in args list we detect the change of the current
      // read mode or read the arg. The reading args algorithm works
      // as a finite-state machine.

      // Detecting new read mode (if it is a valid -key) 
      if (arg == "-l") {
        mode = read_xi;
        continue;
      }

      if (arg == "-s") {
        mode = read_dx;
        continue;
      }

      if (arg == "-i") {
        mode = read_aoi;
        continue;
      }

      if (arg == "-t") {
        if ((mode != read_xh) && (mode != read_comment))
          throw std::invalid_argument(std::string("Unfinished layer!\n")
                                                         +error_msg);
        mode = read_ti;
        continue;
      }

      if (arg == "-c") {
        if ((mode != read_xh) && (mode != read_nt))
          throw std::invalid_argument(std::string("Unfinished layer or theta!\n") + error_msg);
        mode = read_comment;
        continue;
      }

      // Reading data. For invalid date the exception will be thrown
      // with the std:: and catched in the end.
      if (mode == read_xi) {
        xi = std::stod(arg);
        mode = read_xh;
        continue;
      }

      if (mode == read_xh) {
        xh = std::stod(arg);
        mode = read_mi_r;
        continue;
      }

      if (mode == read_mi_r) {
        tmp_mr = std::stod(arg);
        mode = read_mi_i;
        continue;
      }

      if (mode == read_mi_i) {
        mi = std::complex<double>(tmp_mr, std::stod(arg));
        mode = read_mh_r;
        continue;
      }

      if (mode == read_mh_r) {
        tmp_mr = std::stod(arg);
        mode = read_mi_i;
        continue;
      }

      if (mode == read_mh_i) {
        mh = std::complex<double>(tmp_mr, std::stod(arg));
        mode = nothing;
        continue;
      }

      if (mode == read_dx) {
        dx = std::stod(arg);
        mode = nothing;
        continue;
      }

      if (mode == read_aoi) {
        aoi = std::stod(arg);
        mode = nothing;
        continue;
      }

      if (mode == read_ti) {
        ti = std::stod(arg);
        mode = read_tf;
        continue;
      }

      if (mode == read_tf) {
        tf = std::stod(arg);
        mode = read_nt;
        continue;
      }

      if (mode == read_nt) {
        nt = std::stoi(arg);
        Theta.resize(nt);
        S1.resize(nt);
        S2.resize(nt);
        mode = nothing;
        continue;
      }

      if (mode ==  read_comment) {
        comment = arg;
        has_comment = 1;
        continue;
      }
    }

    if (nt < 0) {
      printf("Error reading Theta.\n");
      return -1;
    } else if (nt == 1) {
      Theta[0] = ti*PI/180.0;
    } else {
      for (i = 0; i < nt; i++) {
      Theta[i] = (ti + (double)i*(tf - ti)/(nt - 1))*PI/180.0;
      }
    }

    eccmie::eccMie(xi, xh, mi, mh, dx, aoi, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, &Qpr, &g, &Albedo, S1, S2);

    if (has_comment) {
      printf("%6s, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", comment.c_str(), Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
    } else {
      printf("%+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e, %+.5e\n", Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo);
    }
    
    if (nt > 0) {
      printf(" Theta,         S1.r,         S1.i,         S2.r,         S2.i\n");
      
      for (i = 0; i < nt; i++) {
        printf("%6.2f, %+.5e, %+.5e, %+.5e, %+.5e\n", Theta[i]*180.0/PI, S1[i].real(), S1[i].imag(), S2[i].real(), S2[i].imag());
      }
    }


  } catch( const std::invalid_argument& ia ) {
    // Will catch if  multi_layer_mie fails or other errors.
    std::cerr << "Invalid argument: " << ia.what() << std::endl;
    return -1;
  }  
    return 0;
}


