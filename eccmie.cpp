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

//**********************************************************************************//
// This class implements the algorithm for a multilayered sphere described by:      //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp. 1710-1720.  //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include "eccmie.h"
#include <array>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

using namespace std;

namespace eccmie {
  //helpers
  template<class T> inline T pow2(const T value) {return value*value;}
  int round(double x) {
    return x >= 0 ? (int)(x + 0.5):(int)(x - 0.5);
  }


  //**********************************************************************************//
  // This function emulates a C call to calculate the actual scattering parameters    //
  // and amplitudes.                                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   xi: Size parameter of the inclusion sphere                                     //
  //   xh: Size parameter of the host sphere                                          //
  //   pl: Index of PEC layer (None: -1; Inclusion: 0; Host: 1)                       //
  //   mi: Relative refractive index of the inclusion sphere                          //
  //   mh: Relative refractive index of the host sphere                               //
  //   xd: Size parameter of the shift of the inclusion with respect to the center    //
  //       of the host                                                                //
  //   alpha: Angle of incidence of the light with respect to the Z axis              //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int eccMie(double& xi, double& xh, const int pl, complex<double>& mi, complex<double>& mh, double& xd, double& alpha, const unsigned int nTheta, vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2) {

    if (Theta.size() != nTheta)
        throw invalid_argument("Declared number of sample for Theta is not correct!");
    try {
      EccentricMie eccentric_mie;
      eccentric_mie.SetIncSize(xi);
      eccentric_mie.SetHostSize(xh);
      eccentric_mie.SetIncIndex(mi);
      eccentric_mie.SetHostIndex(mh);
      eccentric_mie.SetIncShift(xd);
      eccentric_mie.SetIncAngle(alpha);
      eccentric_mie.SetAngles(Theta);

      eccentric_mie.SetPECLayer(pl);

      eccentric_mie.RunMieCalculation();

      *Qext = eccentric_mie.GetQext();
      *Qsca = eccentric_mie.GetQsca();
      *Qabs = eccentric_mie.GetQabs();
      *Qbk = eccentric_mie.GetQbk();
      *Qpr = eccentric_mie.GetQpr();
      *g = eccentric_mie.GetAsymmetryFactor();
      *Albedo = eccentric_mie.GetAlbedo();
      S1 = eccentric_mie.GetS1();
      S2 = eccentric_mie.GetS2();
    } catch(const invalid_argument& ia) {
      // Will catch if  eccentric_mie fails or other errors.
      cerr << "Invalid argument: " << ia.what() << endl;
      throw invalid_argument(ia);
      return -1;
    }
    return 0;
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'eccMie' function with fewer    //
  // parameters, it is here mainly for compatibility with older versions of the       //
  // program. Also, you can use it if you neither have a PEC layer nor want to define //
  // any limit for the maximum number of terms.                                       //
  //                                                                                  //
  // Input parameters:                                                                //
  //   xi: Size parameter of the inclusion sphere                                     //
  //   xh: Size parameter of the host sphere                                          //
  //   mi: Relative refractive index of the inclusion sphere                          //
  //   mh: Relative refractive index of the host sphere                               //
  //   xd: Size parameter of the shift of the inclusion with respect to the center    //
  //       of the host                                                                //
  //   alpha: Angle of incidence of the light with respect to the Z axis              //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int eccMie(double& xi, double& xh, complex<double>& mi, complex<double>& mh, double& xd, double& alpha, const unsigned int nTheta, vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2) {
    return eccmie::eccMie(xi, xh, -1, mi, mh, xd, alpha, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'eccMie' function with fewer    //
  // parameters, it is useful if you want to include a PEC layer but not a limit      //
  // for the maximum number of terms.                                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   xi: Size parameter of the inclusion sphere                                     //
  //   xh: Size parameter of the host sphere                                          //
  //   pl: Index of PEC layer (None: -1; Inclusion: 0; Host: 1)                       //
  //   mi: Relative refractive index of the inclusion sphere                          //
  //   mh: Relative refractive index of the host sphere                               //
  //   xd: Size parameter of the shift of the inclusion with respect to the center    //
  //       of the host                                                                //
  //   alpha: Angle of incidence of the light with respect to the Z axis              //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int eccMie(double& xi, double& xh, const int pl, complex<double>& mi, complex<double>& mh, double& xd, double& alpha, const unsigned int nTheta, vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2) {
    return eccmie::eccMie(xi, xh, pl, mi, mh, xd, alpha, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function is just a wrapper to call the full 'eccMie' function with fewer    //
  // parameters, it is useful if you want to include a limit for the maximum number   //
  // of terms but not a PEC layer.                                                    //
  //                                                                                  //
  // Input parameters:                                                                //
  //   xi: Size parameter of the inclusion sphere                                     //
  //   xh: Size parameter of the host sphere                                          //
  //   mi: Relative refractive index of the inclusion sphere                          //
  //   mh: Relative refractive index of the host sphere                               //
  //   xd: Size parameter of the shift of the inclusion with respect to the center    //
  //       of the host                                                                //
  //   alpha: Angle of incidence of the light with respect to the Z axis              //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int eccMie(double& xi, double& xh, complex<double>& mi, complex<double>& mh, double& xd, double& alpha, const unsigned int nTheta, vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2) {
    return eccmie::eccMie(xi, xh, -1, mi, mh, xd, alpha, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
  }


  //**********************************************************************************//
  // This function emulates a C call to calculate complex electric and magnetic field //
  // in the surroundings and inside (TODO) the particle.                              //
  //                                                                                  //
  // Input parameters:                                                                //
  //   xi: Size parameter of the inclusion sphere                                     //
  //   xh: Size parameter of the host sphere                                          //
  //   pl: Index of PEC layer (None: -1; Inclusion: 0; Host: 1)                       //
  //   mi: Relative refractive index of the inclusion sphere                          //
  //   mh: Relative refractive index of the host sphere                               //
  //   xd: Size parameter of the shift of the inclusion with respect to the center    //
  //       of the host                                                                //
  //   alpha: Angle of incidence of the light with respect to the Z axis              //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to 0 (zero) and the function will calculate it.       //
  //   ncoord: Number of coordinate points                                            //
  //   Coords: Array containing all coordinates where the complex electric and        //
  //           magnetic fields will be calculated                                     //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic field at the provided coordinates          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  int eccField(double& xi, double& xh, const int pl, complex<double>& mi, complex<double>& mh, double& xd, double& alpha, const int nmax, const unsigned int ncoord, const vector<double>& Xp_vec, const vector<double>& Yp_vec, const vector<double>& Zp_vec, vector<vector<complex<double> > >& E, vector<vector<complex<double> > >& H) {
    if (Xp_vec.size() != ncoord || Yp_vec.size() != ncoord || Zp_vec.size() != ncoord
        || E.size() != ncoord || H.size() != ncoord)
      throw invalid_argument("Declared number of coords do not fit Xp, Yp, Zp, E, or H!");
    for (auto f:E)
      if (f.size() != 3)
        throw invalid_argument("Field E is not 3D!");
    for (auto f:H)
      if (f.size() != 3)
        throw invalid_argument("Field H is not 3D!");
    try {
      EccentricMie eccentric_mie;
      eccentric_mie.SetIncSize(xi);
      eccentric_mie.SetHostSize(xh);
      eccentric_mie.SetIncIndex(mi);
      eccentric_mie.SetHostIndex(mh);
      eccentric_mie.SetIncShift(xd);
      eccentric_mie.SetIncAngle(alpha);

      //eccentric_mie.SetPECLayer(pl); // TODO add PEC layer to field plotting

      eccentric_mie.SetFieldCoords({Xp_vec, Yp_vec, Zp_vec});
      eccentric_mie.RunFieldCalculation();
      E = eccentric_mie.GetFieldE();
      H = eccentric_mie.GetFieldH();
    } catch(const invalid_argument& ia) {
      // Will catch if  eccentric_mie fails or other errors.
      cerr << "Invalid argument: " << ia.what() << endl;
      throw invalid_argument(ia);
      return - 1;
    }
    return 0;
  }


  // ********************************************************************** //
  // Returns previously calculated Qext                                     //
  // ********************************************************************** //
  double EccentricMie::GetQext() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return Qext_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qabs                                     //
  // ********************************************************************** //
  double EccentricMie::GetQabs() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return Qabs_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qsca                                     //
  // ********************************************************************** //
  double EccentricMie::GetQsca() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return Qsca_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qbk                                      //
  // ********************************************************************** //
  double EccentricMie::GetQbk() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return Qbk_;
  }


  // ********************************************************************** //
  // Returns previously calculated Qpr                                      //
  // ********************************************************************** //
  double EccentricMie::GetQpr() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return Qpr_;
  }


  // ********************************************************************** //
  // Returns previously calculated assymetry factor                         //
  // ********************************************************************** //
  double EccentricMie::GetAsymmetryFactor() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return asymmetry_factor_;
  }


  // ********************************************************************** //
  // Returns previously calculated Albedo                                   //
  // ********************************************************************** //
  double EccentricMie::GetAlbedo() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return albedo_;
  }


  // ********************************************************************** //
  // Returns previously calculated S1                                       //
  // ********************************************************************** //
  vector<complex<double> > EccentricMie::GetS1() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return S1_;
  }


  // ********************************************************************** //
  // Returns previously calculated S2                                       //
  // ********************************************************************** //
  vector<complex<double> > EccentricMie::GetS2() {
    if (!isMieCalculated_)
      throw invalid_argument("You should run calculations before result request!");
    return S2_;
  }


  // ********************************************************************** //
  // Modify scattering (theta) angles                                       //
  // ********************************************************************** //
  void EccentricMie::SetAngles(const vector<double>& angles) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    theta_ = angles;
  }


  // ********************************************************************** //
  // Modify size of inclusion                                               //
  // ********************************************************************** //
  void EccentricMie::SetIncSize(const double inc_size) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_inc_ = inc_size;
  }


  // ********************************************************************** //
  // Modify size of host                                                    //
  // ********************************************************************** //
  void EccentricMie::SetHostSize(const double host_size) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_host_ = host_size;
  }


  // ********************************************************************** //
  // Modify inclusion's shift                                               //
  // ********************************************************************** //
  void EccentricMie::SetIncShift(const double inc_shift) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    shift_inc_ = inc_shift;
  }


  // ********************************************************************** //
  // Modify angle of incidence                                              //
  // ********************************************************************** //
  void EccentricMie::SetIncAngle(const double inc_angle) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    angle_inc_ = inc_angle;
  }


  // ********************************************************************** //
  // Modify refractive index of inclusion                                   //
  // ********************************************************************** //
  void EccentricMie::SetIncIndex(const complex<double> inc_index) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    refractive_index_inc_ = inc_index;
  }


  // ********************************************************************** //
  // Modify refractive index of host                                        //
  // ********************************************************************** //
  void EccentricMie::SetHostIndex(const complex<double> host_index) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    refractive_index_host_ = host_index;
  }


  // ********************************************************************** //
  // Modify coordinates for field calculation                               //
  // ********************************************************************** //
  void EccentricMie::SetFieldCoords(const vector< vector<double> >& coords) {
    if (coords.size() != 3)
      throw invalid_argument("Error! Wrong dimension of field monitor points!");
    if (coords[0].size() != coords[1].size() || coords[0].size() != coords[2].size())
      throw invalid_argument("Error! Missing coordinates for field monitor points!");
    coords_ = coords;
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  void EccentricMie::SetPECLayer(int layer_position) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    if (layer_position < -1 || layer_position > 1)
      throw invalid_argument("Error! Valid values for PEC layers are -1, 0 and 1!");
    PEC_layer_position_ = layer_position;
  }


  // ********************************************************************** //
  // Set maximun number of terms to be used                                 //
  // ********************************************************************** //
  void EccentricMie::SetMaxTerms(int nmax) {
    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    nmax_preset_ = nmax;
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  double EccentricMie::GetSizeParameter() {
    if (size_param_host_ > 0)
      return size_param_host_;
    else
      return 0;
  }


  // ********************************************************************** //
  // Clear layer information                                                //
  // ********************************************************************** //
  void EccentricMie::ClearLayers() {
    complex<double> c_z(0.0, 0.0);

    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;
    size_param_inc_ = 0;
    size_param_host_ = 0;
    shift_inc_ = 0;
    refractive_index_inc_ = c_z;
    refractive_index_host_ = c_z;
  }


  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //
  //                         Computational core
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //


  // ********************************************************************** //
  // Calculate calcNstop - equation (17)                                    //
  // ********************************************************************** //
  void EccentricMie::calcNstop() {
    const double& xL = size_param_host_;
    if (xL <= 8) {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 1);
    } else if (xL <= 4200) {
      nmax_ = round(xL + 4.05*pow(xL, 1.0/3.0) + 2);
    } else {
      nmax_ = round(xL + 4.0*pow(xL, 1.0/3.0) + 2);
    }
  }


  // ********************************************************************** //
  // Maximum number of terms required for the calculation                   //
  // ********************************************************************** //
  void EccentricMie::calcNmax() {
    int ri, riM1;
    const vector<double>& x = {size_param_inc_, size_param_host_};
    const vector<complex<double> >& m = {refractive_index_inc_, refractive_index_host_};
    calcNstop();  // Set initial nmax_ value
    for (unsigned int i = 0; i < 2; i++) {
      if (static_cast<int>(i) > PEC_layer_position_)  // static_cast used to avoid warning
        ri = round(abs(x[i]*m[i]));
      else
        ri = 0;
      nmax_ = max(nmax_, ri);
      // first layer is pec, if pec is present
      if ((i > first_layer) && (static_cast<int>(i - 1) > PEC_layer_position_))
        riM1 = round(abs(x[i - 1]* m[i]));
      else
        riM1 = 0;
      nmax_ = max(nmax_, riM1);
    }
    nmax_ += 15;  // Final nmax_ value
  }


  // ********************************************************************** //
  // Calculate an - equation (5)                                            //
  // ********************************************************************** //
  complex<double> EccentricMie::calc_an(int n, double XL, complex<double> Ha, complex<double> mL,
                                              complex<double> PsiXL, complex<double> ZetaXL,
                                              complex<double> PsiXLM1, complex<double> ZetaXLM1) {

    complex<double> Num = (Ha/mL + n/XL)*PsiXL - PsiXLM1;
    complex<double> Denom = (Ha/mL + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }


  // ********************************************************************** //
  // Calculate bn - equation (6)                                            //
  // ********************************************************************** //
  complex<double> EccentricMie::calc_bn(int n, double XL, complex<double> Hb, complex<double> mL,
                                              complex<double> PsiXL, complex<double> ZetaXL,
                                              complex<double> PsiXLM1, complex<double> ZetaXLM1) {

    complex<double> Num = (mL*Hb + n/XL)*PsiXL - PsiXLM1;
    complex<double> Denom = (mL*Hb + n/XL)*ZetaXL - ZetaXLM1;

    return Num/Denom;
  }

  // ********************************************************************** //
  // Calculates S1 - equation (25a)                                         //
  // ********************************************************************** //
  complex<double> EccentricMie::calc_S1(int n, complex<double> an, complex<double> bn,
                                              double Pi, double Tau) {
    return double(n + n + 1)*(Pi*an + Tau*bn)/double(n*n + n);
  }


  // ********************************************************************** //
  // Calculates S2 - equation (25b) (it's the same as (25a), just switches  //
  // Pi and Tau)                                                            //
  // ********************************************************************** //
  complex<double> EccentricMie::calc_S2(int n, complex<double> an, complex<double> bn,
                                              double Pi, double Tau) {
    return calc_S1(n, an, bn, Tau, Pi);
  }


  //**********************************************************************************//
  // This function calculates the logarithmic derivatives of the Riccati-Bessel       //
  // functions (D1 and D3) for a complex argument (z).                                //
  // Equations (16a), (16b) and (18a) - (18d)                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate D1 and D3                                      //
  //   nmax_: Maximum number of terms to calculate D1 and D3                          //
  //                                                                                  //
  // Output parameters:                                                               //
  //   D1, D3: Logarithmic derivatives of the Riccati-Bessel functions                //
  //**********************************************************************************//
  void EccentricMie::calcD1D3(const complex<double> z,
                               vector<complex<double> >& D1,
                               vector<complex<double> >& D3) {

    // Downward recurrence for D1 - equations (16a) and (16b)
    D1[nmax_] = complex<double>(0.0, 0.0);
    const complex<double> zinv = complex<double>(1.0, 0.0)/z;

    for (int n = nmax_; n > 0; n--) {
      D1[n - 1] = static_cast<double>(n)*zinv - 1.0/(D1[n] + static_cast<double>(n)*zinv);
    }

    if (abs(D1[0]) > 100000.0)
      throw invalid_argument("Unstable D1! Please, try to change input parameters!\n");

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_[0] = 0.5*(1.0 - complex<double>(cos(2.0*z.real()), sin(2.0*z.real()))
                     *exp(-2.0*z.imag()));
    D3[0] = complex<double>(0.0, 1.0);
    for (int n = 1; n <= nmax_; n++) {
      PsiZeta_[n] = PsiZeta_[n - 1]*(static_cast<double>(n)*zinv - D1[n - 1])
                                   *(static_cast<double>(n)*zinv - D3[n - 1]);
      D3[n] = D1[n] + complex<double>(0.0, 1.0)/PsiZeta_[n];
    }
  }


  //**********************************************************************************//
  // This function calculates the Riccati-Bessel functions (Psi and Zeta) for a       //
  // complex argument (z).                                                            //
  // Equations (20a) - (21b)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   z: Complex argument to evaluate Psi and Zeta                                   //
  //   nmax: Maximum number of terms to calculate Psi and Zeta                        //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Psi, Zeta: Riccati-Bessel functions                                            //
  //**********************************************************************************//
  void EccentricMie::calcPsiZeta(complex<double> z,
                                  vector<complex<double> >& Psi,
                                  vector<complex<double> >& Zeta) {

    complex<double> c_i(0.0, 1.0);
    vector<complex<double> > D1(nmax_ + 1), D3(nmax_ + 1);

    // First, calculate the logarithmic derivatives
    calcD1D3(z, D1, D3);

    // Now, use the upward recurrence to calculate Psi and Zeta - equations (20a) - (21b)
    Psi[0] = sin(z);
    Zeta[0] = sin(z) - c_i*cos(z);
    for (int n = 1; n <= nmax_; n++) {
      Psi[n]  =  Psi[n - 1]*(static_cast<double>(n)/z - D1[n - 1]);
      Zeta[n] = Zeta[n - 1]*(static_cast<double>(n)/z - D3[n - 1]);
    }
  }


  //**********************************************************************************//
  // This function calculates Pi and Tau for a given value of cos(Theta).             //
  // Equations (26a) - (26c)                                                          //
  //                                                                                  //
  // Input parameters:                                                                //
  //   nmax_: Maximum number of 'n' terms to calculate Pi and Tau                     //
  //   nmax_: Maximum number of 'm' terms to calculate Pi and Tau                     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Pi, Tau: Angular functions Pi and Tau, as defined in equations (26a) - (26c)   //
  //**********************************************************************************//
  void EccentricMie::calcPiTau(const double& costheta,
                                vector<vector<double> >& Pi, vector<vector<double> >& Tau) {

    double sintheta = sqrt(1.0 - costheta*costheta);

    Pi[0][0] = 1.0/sintheta;
    Pi[0][1] = costheta/sintheta;
    Pi[1][0] = 0.0;
    Pi[1][1] = 1.0;

    Tau[0][0] = 0.0;
    Tau[0][1] = costheta*Pi[0][1] - Pi[0][0];
    Tau[1][0] = 0.0;
    Tau[1][1] = costheta*Pi[1][1] - 2*Pi[1][0];

    // If |m| > n => Pi[m][n] = Tau[m][n] = 0
    for (int m = 2; m <= nmax_; m++) {
      Pi[m][0] = 0.0;
      Pi[m][1] = 0.0;

      Tau[m][0] = 0.0;
      Tau[m][1] = 0.0;
    }

    for (int m = 0; m <= nmax_; m++) {
      for (int n = 2; n <= nmax_; n++) {
        if (m > n) {
          Pi[m][n] = 0.0;
          Tau[m][n] = 0.0;
        } else {
          Pi[m][n] = ((n + n - 1)*costheta*Pi[m][n - 1] - (n + m - 1)*Pi[m][n - 2])/(n - m);
          Tau[m][n] = n*costheta*Pi[m][n] - (n + m)*Pi[m][n - 1];
        }
      }
    }
  }  // end of EccentricMie::calcPiTau(...)


  //**********************************************************************************//
  // This function calculates vector spherical harmonics (eq. 4.50, p. 95 BH),        //
  // required to calculate the near-field parameters.                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //   zn: Either the spherical Bessel or Hankel function of first kind               //
  //   dzn: Derivative of zn                                                          //
  //   Pi, Tau: Angular functions Pi and Tau                                          //
  //   n: Order of vector spherical harmonics                                         //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Mo1n, Me1n, No1n, Ne1n: Complex vector spherical harmonics                     //
  //**********************************************************************************//
  void EccentricMie::calcSpherHarm(const double Rho, const double Theta, const double Phi,
                                    const complex<double>& zn, const complex<double>& dzn,
                                    const double& Pi, const double& Tau, const double& n,
                                    vector<complex<double> >& Mo1n, vector<complex<double> >& Me1n,
                                    vector<complex<double> >& No1n, vector<complex<double> >& Ne1n) {

    // using eq 4.50 in BH
    complex<double> c_zero(0.0, 0.0);
    complex<double> deriv = Rho*dzn + zn;

    Mo1n[0] = c_zero;
    Mo1n[1] = cos(Phi)*Pi*zn;
    Mo1n[2] = -sin(Phi)*Tau*zn;
    Me1n[0] = c_zero;
    Me1n[1] = -sin(Phi)*Pi*zn;
    Me1n[2] = -cos(Phi)*Tau*zn;
    No1n[0] = sin(Phi)*(n*n + n)*sin(Theta)*Pi*zn/Rho;
    No1n[1] = sin(Phi)*Tau*deriv/Rho;
    No1n[2] = cos(Phi)*Pi*deriv/Rho;
    Ne1n[0] = cos(Phi)*(n*n + n)*sin(Theta)*Pi*zn/Rho;
    Ne1n[1] = cos(Phi)*Tau*deriv/Rho;
    Ne1n[2] = -sin(Phi)*Pi*deriv/Rho;
  }  // end of EccentricMie::calcSpherHarm(...)


  //**********************************************************************************//
  // This function calculates the scattering coefficients required to calculate       //
  // both the near- and far-field parameters.                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it.             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   an, bn: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void EccentricMie::ScattCoeffs() {

    isScaCoeffsCalc_ = false;


    isScaCoeffsCalc_ = true;
  }  // end of EccentricMie::ScattCoeffs(...)


  complex<double> EccentricMie::calcTransCn(int n, int np, complex<double> Cnnpm1,
                                            complex<double> Cnm1np, complex<double> Cnnpp1) {

    complex<double> term1, term2, term3, factor;

    factor = sqrt(double(n + n + 3)/double(np + np + 1))/double(n + 1);
    term1  = double(np)*Cnnpm1*sqrt(double(n + n + 1)/double(np + np - 1));
    term2  = double(n)*Cnm1np*sqrt(double(np + np + 1)/double(n + n - 1));
    term3  = -double(np + 1)*Cnnpp1*sqrt(double(n + n + 1)/double(np + np + 3));

    return factor*(term1 + term2 + term3);
  }

  complex<double> EccentricMie::calcTransCm(int m, int n, int np,
                                            complex<double> zd, complex<double> Cnnp,
                                            complex<double> Cnnpp1, complex<double> Cnnpm1) {

    complex<double> term1, term2, term3, factor;

    factor = 1/sqrt((n-m+1)*(n+m)*(np+np+1));
    term1  = Cnnp*sqrt((np-m+1)*(np+m)*(np+np+1));
    term2  = -zd*Cnnpp1*sqrt((np-m+2)*(np-m+1)/(np+np+3));
    term3  = -zd*Cnnpm1*sqrt((np+m)*(np+m-1)/(np+np-1));

    return factor * (term1 + term2 + term3);
  }

  complex<double> EccentricMie::calcTransA(int m, int np,
                                           complex<double> zd, complex<double> Cnnp,
                                           complex<double> Cnnpp1, complex<double> Cnnpm1) {

    complex<double> term1, term2;

    term1  = -zd*Cnnpp1*sqrt((np - m + 1)*(np + m + 1)/double((np + np + 1)*(np + np + 3)))/double(np + 1);
    term2  = -zd*Cnnpm1*sqrt(double((np - m)*(np + m))/double((np + np + 1)*(np + np - 1)))/double(np);

    return  Cnnp + term1 + term2;
  }

  complex<double> EccentricMie::calcTransB(int m, int np,
                                           complex<double> zd, complex<double> Cnnp) {

    complex<double> c_i(0.0, 1.0);
    return  -c_i*double(m)*zd*Cnnp/double(np*(np + 1));
  }

  complex<double> EccentricMie::calcTU(complex<double> m1, complex<double> m2,
                                       complex<double> Q, complex<double> Trans,
                                       complex<double> Psi_z2, complex<double> D1_z2,
                                       complex<double> D3_x2,
                                       complex<double> Zeta1_z2, complex<double> D3_z2) {

    return Trans*(m1*D3_x2*(2.0*Psi_z2 + Zeta1_z2*Q) - m2*(2.0*Psi_z2*D1_z2 + Zeta1_z2*D3_z2*Q)) ;
  }

  complex<double> EccentricMie::calcQ(complex<double> m1, complex<double> m2,
                                      complex<double> Psi, complex<double> D1,
                                      complex<double> Zeta1, complex<double> D3) {

    return  2.0*Psi*D1*(m1/Zeta1 - m2)/(m2*D1 - m1*D3);
  }


  //**********************************************************************************//
  // This function calculates the actual scattering parameters and amplitudes         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax_: Maximum number of multipolar expansion terms to be used for the         //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void EccentricMie::RunMieCalculation() {

    isExpCoeffsCalc_ = false;
    isScaCoeffsCalc_ = false;
    isMieCalculated_ = false;

    complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    complex<double> m1, m2, z1, z2, z3, zd;
    double x2;

    m1 = refractive_index_inc_;
    m2 = refractive_index_host_;
    x2 = size_param_host_;
    z1 = m1*size_param_inc_;
    z2 = m2*size_param_host_;
    z3 = m2*size_param_inc_;
    zd = m2*shift_inc_;

    // Estimate maximum number of iterations
    calcNmax();

    vector<complex<double> > Psi_x2(nmax_ + 1, c_zero), D1_x2(nmax_ + 1, c_zero);
    vector<complex<double> > Zeta1_x2(nmax_ + 1, c_zero), D3_x2(nmax_ + 1, c_zero);

    vector<complex<double> > Psi_z2(nmax_ + 1, c_zero), D1_z2(nmax_ + 1, c_zero);
    vector<complex<double> > Zeta1_z2(nmax_ + 1, c_zero), D3_z2(nmax_ + 1, c_zero);

    vector<complex<double> > Psi_z3(nmax_ + 1, c_zero), D1_z3(nmax_ + 1, c_zero);
    vector<complex<double> > Zeta1_z3(nmax_ + 1, c_zero), D3_z3(nmax_ + 1, c_zero);
    
    vector<complex<double> > cnm(nmax_ + 1, c_zero), dnm(nmax_ + 1, c_zero);

    // For z=x2
    // First, calculate the logarithmic derivatives
    calcD1D3(x2, D1_x2, D3_x2);
    // Now, calculate Riccati-Bessel functions
    calcPsiZeta(x2, Psi_x2, Zeta1_x2);

    // For z=z2
    // First, calculate the logarithmic derivatives
    calcD1D3(z2, D1_z2, D3_z2);
    // Now, calculate Riccati-Bessel functions
    calcPsiZeta(z2, Psi_z2, Zeta1_z2);

    // For z=z3
    // First, calculate the logarithmic derivatives
    calcD1D3(z3, D1_z3, D3_z3);
    // Now, calculate Riccati-Bessel functions
    calcPsiZeta(z3, Psi_z3, Zeta1_z3);

//  .......................................................................
//  .     Now we must solve for the Quality factors by using equations    .
//  .     7.1 and 7.2.                                                      .
//  .......................................................................

    vector<complex<double> > Qr(nmax_ + 1, c_zero), Qs(nmax_ + 1, c_zero);
    int n;

    if (PEC_layer_position_ < 0)
      for (n = 0; n < nmax_; n++) {
        Qr[n] = calcQ(m1, m2, Psi_z3[n], D1_z3[n], Zeta1_z3[n], D3_z3[n]); //eq. 7.1
        Qs[n] = calcQ(m2, m1, Psi_z3[n], D1_z3[n], Zeta1_z3[n], D3_z3[n]); //eq. 7.2
      }
    else // The inclusion is a perfect conductor
      for (n = 0; n < nmax_; n++) {
        Qr[n] = -2.0*Psi_z3[n]/Zeta1_z3[n];
        Qs[n] = -2.0*Psi_z3[n]*D1_z3[n]/Zeta1_z3[n]/D3_z3[n];
      }


//  ........................................................................
//  .    Calculate translation coefficients when the separation of the     .
//  .    center to center radii is d.  We will first find the bessel       .
//  .    function which describes the separation, m2xd=zd, used in 10.1 and 10.2.
//  ........................................................................


    vector<complex<double> > Psi_d(nmax_ + 1, c_zero), D1_d(nmax_ + 1, c_zero);
    vector<complex<double> > Zeta1_d(nmax_ + 1, c_zero), D3_d(nmax_ + 1, c_zero);

    // For z=zd
    // First, calculate the logarithmic derivatives
    calcD1D3(zd, D1_d, D3_d);
    // Now, calculate Riccati-Bessel functions
    calcPsiZeta(zd, Psi_d, Zeta1_d);

//  ...............................................
//  .    using 10.1 and 10.2 first.                   .
//  .    C00 = c0,0(np) and                     .
//  .    C0m1 = c-1,0(np) are the starters      .
//  ...............................................

    vector<complex<double> > C00(nmax_), C0m1(nmax_);
    vector<vector<complex<double> > > Cmn(nmax_), Cmm1n(nmax_), TransAm(nmax_), TransBm(nmax_);

    int np;
    for (n = 0; n < nmax_; n++) {
      Cmn[n].resize(2*nmax_ + 1 - n);
      Cmm1n[n].resize(2*nmax_ + 1 - n);

      TransAm[n].resize(nmax_);
      TransBm[n].resize(nmax_);
    }

    for (np = 0; np < 2*nmax_ + 1; np++) {
      C00[np] = sqrt(np + np + 1.0)*Psi_d[np]/zd ;
      C0m1[np] = -C00[np];
      Cmn[0][np] = C00[np];
    }

//  .................................................
//  .     Now use equation 10.3, this is where we     .
//  .     perform the switcharoos and recursively   .
//  .     solve for the translation coeff. for m=0  .
//  .     The result is stored into the Cmn      .
//  .................................................

    for (n = 0; n < nmax_; n++) // [0..nmax_-1]
      for (np = 0; np < 2*nmax_ + 1 - n; np++)
        Cmn[n + 1][np] = calcTransCn(n, np, Cmn[n][np - 1], Cmn[n - 1][np], Cmn[n][np + 1]);

    // Calculate angular functions Pi and Tau
    vector<vector<double> > Pi(nmax_), Tau(nmax_);
    for (int i = 0; i < nmax_; i++) {
      Pi[i].resize(nmax_);
      Tau[i].resize(nmax_);
    }

    calcPiTau (cos(PI_*angle_inc_/180.), Pi, Tau);

//  .......................................................
//  . Now that we have the elements in the matrix         .
//  . Cmn[0][n][np], we can start the            .
//  . iterations in m to get the rest of the translation  .
//  . coeff. Then we'll have the complete set of Cmn[m][n][np] .
//  .......................................................

    for (int m = 0; m < nmax_; m++) { // main loop for m
//  ................................................................
//  . Now to find the rest of the matrix in C(n,m)np.  For each    .
//  . value of m that is greater than 0, say 1, we can find the    .
//  . matrix C(n,1)np and store that into Cmn(n,np); then       .
//  . the next iteration of m, i.e. m=2, we can find C(n,2)np      .
//  . by already having C(n,1)np stored in the matrix Cmm1n(n,np).
//  . Use equation 10.4                                              .
//  ................................................................

      if (m > 0) {
        //save the last Cmn in Cmm1n
        for (n = 0; n < nmax_; n++)
          for (np = 0; np < 2*nmax_ + 1; np++)
            Cmm1n[n][np] = Cmn[n][np];

        for (n = m; n < nmax_; n++)
          for (np = m; np < 2 * nmax_ - m; np++)
            Cmn[n][np] = calcTransCm(m, n, np, zd, Cmm1n[n][np], Cmm1n[n][np + 1], Cmm1n[n][np - 1]);

      } // if m>0

//    ..................................................
//    . Cmn[n][np] now has the most recent values    .
//    . for the most recent value of m.                .
//    ..................................................

//    ............................................................
//    . Now that we have the matrix which contains the C(n,m)np  .
//    . we can go ahead and find the translation coefficients    .
//    . A(n,m)np and B(n,m)np by using eqns 9.1 & 9.2. Note,       .
//    . we do need to store TransAm and TransBm into a matrix    .
//    . like we did for the C(n,m)np for later manipulation      .
//    ............................................................

      for (n = m; n < nmax_; n++)
        for (np = m; np < nmax_; np++) {
          TransAm[n][np] = calcTransA(m, np, zd, Cmn[n][np], Cmn[n][np + 1], Cmn[n][np - 1]);
          TransBm[n][np] = calcTransB(m, np, zd, Cmn[n][np]);
        }

//    ...............................................................
//    .     Now we will calculate the incident field coefficients   .
//    .     for the TE case and store them in cnm(i) and dnm(i)       .
//    ...............................................................

      complex<double> factor;
      for (n = m; n < nmax_; n++) {
        factor = 2.0*pow(c_i, n + 2)/double(n*n + n);
        cnm[n] = factor*Tau[m][n];
        dnm[n] = factor*Pi[m][n];
      }

//    Fill the matrix and load de solution vectors (cnm & dnm) in the last column

      vector<vector<complex<double> > > matrix;
      complex<double> gamma_n;
      // matriz [nmax_-m][nmax_-m+1]
      int numel, j;

      numel = nmax_ - m;

      int i = 0;
      matrix.resize(numel);
      for (n = m; n < nmax_; n++) {
        matrix[i].resize(numel + 1);
        j = 0;
        for (np = m; np < nmax_; np++) {
          matrix[i][j] = calcTU(1.0, m2, Qr[np], TransAm[n][np], Psi_z2[n], D1_z2[n], D3_x2[n], Zeta1_z2[n], D3_z2[n]); // eq. 16.2
          matrix[i][j + numel] = calcTU(1.0, m2, Qs[np], TransBm[n][np], Psi_z2[n], D1_z2[n], D3_x2[n], Zeta1_z2[n], D3_z2[n]); // eq. 16.4
          matrix[i + numel][j] = calcTU(m2, 1.0, Qr[np], TransBm[n][np], Psi_z2[n], D1_z2[n], D3_x2[n], Zeta1_z2[n], D3_z2[n]); // eq. 16.3
          matrix[i + numel][j + numel] = calcTU(m2, 1.0, Qs[np], TransAm[n][np], Psi_z2[n], D1_z2[n], D3_x2[n], Zeta1_z2[n], D3_z2[n]); // eq. 16.5
          j++;
        } //np
        gamma_n = m2*Psi_x2[n]*(D3_x2[n] - D1_x2[n]); //Equation 16.1
        matrix[i][j] = cnm[n]*gamma_n;
        matrix[i + numel][j] = dnm[n]*gamma_n;
        i++;
      } //n

//    ................................................................
//    . THE MATRIX IS FILLED!!  For TE case                          .
//    ................................................................

    } // loop m


/*    // # of iterations for the small sphere and loops
    nsm = round(abs(z3) + 4.0*pow(abs(z3), 1.0/3.0) + 2.0);





    // # of iterations for the big sphere and loops
    nmax_ = round(abs(x2) + 4.0*pow(abs(x2), 1.0/3.0) + 2.0);
    // # of iterations for the small sphere and loops
    nsm = round(abs(z3) + 4.0*pow(abs(z3), 1.0/3.0) + 2.0);

    int i, j, jj, jjj, n, m, mn, np, nmax_, nsm, ntot, nel;
    double Qext1, Qsca1, Qabs, gg, f_angle, angle_inc__0, tau, pie, dd;
    double alpha, d_ang, b_angle, angledd, theta;
    // For the refractive indices and k-vectors
    double x2;
    complex<double> m2, m1, z2, z3, z1, zd;
    // For the MEDIUM, which in our case is just plain old air.
    vector<complex<double> > Jb_0(nmax_ + 2), tmp(nmax_ + 2), psi_0(nmax_ + 2), dpsi_0(nmax_ + 2);
    vector<complex<double> > Hk_0(nmax_ + 2), ctmp(nmax_ + 2), Zeta_0(nmax_ + 2), dzeta_0(nmax_ + 2);
    // For the HOST SPHERE, all the arrays need be at least nmax_
    vector<complex<double> > Hk1_1(nmax_ + 2), Zeta1_z2(nmax_ + 2), D3_z2(nmax_ + 2), Hk2_1(nmax_ + 2), Zeta2_1(nmax_ + 2), D4_1(nmax_ + 2);
    // For the INCLUSION SPHERE
    vector<complex<double> > Jb_2(nmax_ + 2), Psi_z3(nmax_ + 2), D1_z3(nmax_ + 2), Hk1_2(nmax_ + 2), Zeta1_z3(nmax_ + 2), D3_z3(nmax_ + 2);
    vector<complex<double> > Hk2_2(nmax_ + 2), Zeta2_2(nmax_ + 2), D4_2(nmax_ + 2);
    // For the separation, and Jb_d need be at least 4*nmax_
    vector<complex<double> > Jb_d(4*nmax_ + 2);
    // For the Internal Coefficients
    vector<vector<complex<double> > > TintTE(nmax_ + 1), UintTE(nmax_ + 1), TintTM(nmax_ + 1), UintTM(nmax_ + 1);
    // For the SCATTERING
    vector<vector<complex<double> > > c_TE(nmax_ + 1), d_TE(nmax_ + 1), c_TM(nmax_ + 1), d_TM(nmax_ + 1);
    vector<complex<double> > Qr(nmax_ + 1), Qs(nmax_ + 1);
    vector<complex<double> > s1_TE(nmax_ + 1), s2_TM(nmax_ + 1), s3_TE(nmax_ + 1), s4_TM(nmax_ + 1);
    complex<double> s1, s2, s3, s4, sum1 , sum2, fact1, fact2, fact3, fact4;
    // For the NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
    vector<double > legpol, legpol1, legpol2;
    double angle, x;
    // For the TRANSLATION COEFFICIENTS
    vector<complex<double> > c00(4*nmax_ + 2), cm10(4*nmax_ + 2);
    vector<vector<complex<double> > > c00m(4*nmax_ + 2), cm10m(4*nmax_ + 2);
    complex<double> TransA, TransB;
    vector<vector<complex<double> > > TransAm, TransBm;
    complex<double> c, cc, ifact;
    // For the MATRICES IN THE LU_DEC0MPOSTION
    vector<int > indx(nmax_ + nmax_ + 2);
    vector<complex<double> > b, bd;
    vector<vector<complex<double> > > matrix, matrixlu;
    vector<complex<double> > a_TE(nmax_ + 1), b_TE(nmax_ + 1);
    // MUELLER MATRIX
    vector<double > mm(4);
    for (int n = 0; n < 4; n++) {
      mm[n].resize(4);
    }

    // Set and initialize arrays
    for (int n = 0; n <= nmax_; n++) {
      TintTE[n].resize(nmax_ + 1);
      UintTE[n].resize(nmax_ + 1);
      TintTM[n].resize(nmax_ + 1);
      UintTM[n].resize(nmax_ + 1);

      c_TE[n].resize(nmax_ + 1);
      d_TE[n].resize(nmax_ + 1);
      c_TM[n].resize(nmax_ + 1);
      d_TM[n].resize(nmax_ + 1);

      TransAm.resize(nmax_ + 1);
      TransBm.resize(nmax_ + 1);
    }

    for (int n = 0; n < (4*nmax_ + 2); n++) {
      c00m[n].resize(4*nmax_ + 2);
      cm10m[n].resize(4*nmax_ + 2);
    }

    for (int i = 0; i < (nmax_ + 2); i++) {
      Jb_0[i] = 0.0;
      tmp[i] = 0.0;
      psi_0[i] = 0.0;
      dpsi_0[i] = 0.0;
      Hk_0[i] = 0.0;
      ctmp[i] = 0.0;
      Zeta_0[i] = 0.0;
      dzeta_0[i] = 0.0;
      Hk1_1[i] = 0.0;
      Zeta1_z2[i] = 0.0;
      D3_z2[i] = 0.0;
      Hk2_1[i] = 0.0;
      Zeta2_1[i] = 0.0;
      D4_1[i] = 0.0;
      Jb_2[i] = 0.0;
      Psi_z3[i] = 0.0;
      D1_z3[i] = 0.0;
      Hk1_2[i] = 0.0;
      Zeta1_z3[i] = 0.0;
      D3_z3[i] = 0.0;
      Hk2_2[i] = 0.0;
      Zeta2_2[i] = 0.0;
      D4_2[i] = 0.0;
    }

    // ................................................................
    // . THE MATRIX IS FILLED!!  That was the heart of the problem    .
    // ................................................................
    // ...............................................................
    // . Now we can load up the solution vectors a_TE[n] and b_TE[n]     .
    // . and start the LU Decomposition to get t[n] and u[n]         .
    // . Note that a_TE[n] and b_TE[n] contain the incident angle        .
    // . Utilize equations 37 - 43.                                    .
    // ...............................................................
    // .............
    // For TE case .
    // .............
    for (int n = m; n <= nmax_; n++) {
      c = dzeta_0[n]*psi_0[n];
      c = m2*(c - dpsi_0[n]*Zeta_0[n]);
      b[n - m] = a_TE[n]*c;
      b[n - m + nel] = b_TE[n]*c;
      bd[n - m] = b[n - m];
      bd[n - m + nel] = b[n - m + nel];
    }
    //writematrix('orig.txt', matrixlu);
    dd = ludcmp(indx, matrixlu);
    //writematrix('final.txt', matrixlu);
    lubksb(indx, matrixlu, b);
    improve(indx, bd, matrix, matrixlu, b);
    for (int n = m; n <= nmax_; n++) {
      TintTE[m][n] = b[n - m];
      UintTE[m][n] = b[n - m + nel];
    }
    // ...................................................................
    // .  Now we must find the Scattering coefficients c_TE and d_TE,    .
    // .  and from these coeff we can determine the scattering matrices. .
    // .  The equations are  33 and 35 (or 34 and 36).                   .
    // ...................................................................
    for (int n = m; n <= nmax_; n++) {
      sum1 = 0.0;
      sum2 = 0.0;
      for (int np = m; np <= nmax_; np++) {
        fact1 = D4_1[n] + Qr[np]*D3_z2[n];
        fact1 = TransAm[np][n]*TintTE[m, np]*fact1;
        fact2 = D4_1[n] + Qs[np]*D3_z2[n];
        fact2 = TransBm[np][n]*UintTE[m, np]*fact2;
        sum1 = fact1 + fact2 + sum1;

        fact3 = Zeta2_1[n] + Qr[np]*Zeta1_z2[n];
        fact3 = TransBm[np][n]*TintTE[m, np]*fact3;
        fact4 = Zeta2_1[n] + Qs[np]*Zeta1_z2[n];
        fact4 = TransAm[np][n]*UintTE[m, np]*fact4;
        sum2 = sum2 + fact3 + fact4;
      }
      c_TE[m][n] = (sum1 - a_TE[n]*dpsi_0[n])/dzeta_0[n];
      d_TE[m][n] = (sum2 - b_TE[n]*psi_0[n])/Zeta_0[n];
    }
    // .................
    // FOR TM CASE!!! .
    // .................
    for (int n = m; n <= nmax_; n++) {
      c = dzeta_0[n]*psi_0[n];
      c  = m2*(c - dpsi_0[n]*Zeta_0[n]);
      b[n - m] = _i*b_TE[n]*c;
      b[n - m + nel] = _i*a_TE[n]*c;
      bd[n - m] = b[n - m];
      bd[n - m + nel] = b[n - m + nel];
    }
    lubksb(indx, matrixlu, b);
    improve(indx, bd, matrix, matrixlu, b);
    for (int n = m; n <= nmax_; n++) {
      TintTM[m][n] = b[n - m];
      UintTM[m][n] = b[n - m + nel];
    }
    // ...................................................................
    // .  Now we must find the Scattering coefficients c_TM and d_TM,    .
    // .  and from these coeff we can determine the scattering matrices. .
    // .  The equations are gotten from 33 and 35                        .
    // ...................................................................
    for (int n = m; n <= nmax_; n++) {
      sum1 = 0.0;
      sum2 = 0.0;
      for (int np = m; np <= nmax_; np++) {
        fact1 = D4_1[n] + Qr[np]*D3_z2[n];
        fact1 = TransAm[np][n]*TintTM[m, np]*fact1;
        fact2 = D4_1[n] + Qs[np]*D3_z2[n];
        fact2 = TransBm[np][n]*UintTM[m, np]*fact2;
        sum1 = fact1 + fact2 + sum1;

        fact3 = Zeta2_1[n] + Qr[np]*Zeta1_z2[n];
        fact3 = TransBm[np][n]*TintTM[m, np]*fact3;
        fact4 = Zeta2_1[n] + Qs[np]*Zeta1_z2[n];
        fact4 = TransAm[np][n]*UintTM[m, np]*fact4;
        sum2 = sum2 + fact3 + fact4;
      }
      c_TM[m][n] = dpsi_0[n]*_i*b_TE[n];
      c_TM[m][n] = (sum1 - c_TM[m][n])/dzeta_0[n];
      d_TM[m][n] = _i*a_TE[n]*psi_0[n];
      d_TM[m][n] = (sum2 - d_TM[m][n])/Zeta_0[n];
    }
    // .......................................................................
    // . We can now calculate the Efficiencies (Qext1, Qsca1, Qabs) by using   .
    // . the equations 58 - 60.                                              .
    // .......................................................................
    for (int n = 1; n <= nmax_; n++) {
      tau = c_TE[m][n].r*c_TE[m][n].r + c_TE[m][n].i*c_TE[m][n].i;
      tau = tau + d_TE[m][n].r*d_TE[m][n].r + d_TE[m][n].i*d_TE[m][n].i;
      pie = c_TM[m][n].r*c_TM[m][n].r + c_TM[m][n].i*c_TM[m][n].i;
      pie = pie + d_TM[m][n].r*d_TM[m][n].r + d_TM[m][n].i*d_TM[m][n].i;
      x = (n*n + n)*(tau + pie);

      s1 = c_TE[m][n]*conjg(a_TE[n]);
      s1 = s1 + d_TE[m][n]*conjg(b_TE[n]);
      s2 = c_TM[m][n]*conjg(_i*b_TE[n]);
      s2 = s2 + d_TM[m][n]*conjg(_i*a_TE[n]);
      s3 = (n*n + n)*(s2 + s1);

      if (m == 0) {
        x = 2.0*x;
        s3 = -2.0*s3;
      } else {
        x = 4.0*x;
        s3 = -4.0*s3;
      }

      Qsca1 = Qsca1 + x;
      Qext1 = Qext1 + s3.r;
    }
    // .......................................................................
    // .  That finishes one loop in m, now we move to the next m and repeat  .
    // .  the same procedure all over again, storing our internal results    .
    // .  into T_te, U_te, T_tm, U_tm, and storing the scattering results    .
    // .  into C_te, D_te, C_tm, D_tm.                                       .
    // .......................................................................
  }
  // ............................................
  // DONE!!  The sum over m is now complete!!! .
  // ............................................
  Qext = Qext1/x2/x2;
  Qsca = Qsca1/x2/x2;
  Qabs = Qext - Qsca;

  // COMMENT THIS BLOCK OUT IF YOU DON'T WANT TO FIND
  // THE ASYMMETRY PARAMETER.  THIS FIND_G SUBROUTINE
  // IS VERY TIME CONSUMING (Use equation 61)
  alpha = ai;
  gg = find_g(nmax_, isize1, alpha, c_TE, d_TE, c_TM, d_TM);
  gg = gg/x2/x2/Qsca;
  // .......................................................................
  // . Here's the meat of the problem.  We will calculate the amplitude    .
  // . scattering matrices for the fields that we've found.  The angle     .
  // . theta is our observation angle.  Num is the number of angles we     .
  // . want in determining the scattering intensity.  The rest of this     .
  // . should look almost like the Mie Code found in Bohren and Huffman    .
  // . We have here 3 nested loops.  The first one is to determine the     .
  // . angles we want to look at, the second one is the sum over m, and    .
  // . the third one sums over n.                                          .
  // . If we want to have the observation angle's resolution of 1 degree   .
  // . set noai = 360, if we want finer resolution, let noai be a bigger.
  // . number such as 720 (to have resolution of 0.5degrees)               .
  // .......................................................................
  angle_inc__0 = angle_inc_*pi/180.0;
  i = 0;
  // noai = 1;
  d_ang = 360.0/noai;
  b_angle = ai + 180.0;
  f_angle = ai;
  for (int j = 1; j <= noai + 1; j++) {
    angle = d_ang*(j - 1);
    angledd = angle;
    if (angle > 180.0) {
      angle = 360.0 - angle;
      angle_inc_ = pi + angle_inc__0;
    } else {
      angle_inc_ = angle_inc__0;
    }
    i = i + 1;
    s1_TE[i] = 0.0;
    s3_TE[i] = 0.0;
    s2_TM[i] = 0.0;
    s4_TM[i] = 0.0;
    theta = angle*pi/180.0;
    x = cos(theta);
    for (int m = 0; m <= nmax_; m++) {
      // ...........................................................
      // .  The next few lines will calc the Legendre polynomial   .
      // .  at the angle theta                                     .
      // ...........................................................
      legpol2 = nplgndr(m, nmax_ + 1, x);
      legpol1 = nplgndr(m + 1, nmax_ + 1, x);
      if (m = 0) {
        legpol.resize(legpol1.size() + 1);
        for (int jj = 0; jj <= legpol1.size(); jj++)
          legpol[jj] = -legpol1[jj];
      } else
        legpol = nplgndr(m - 1, nmax_ + 1, x);

      for (int n = m; n <= nmax_; n++) {
        // .........................................................
        // .  Calculate the corresponding pie and tau              .
        // .........................................................
        tau = sqrt((n + m + 1.0)*(n - m))*legpol1[n];
        pie = tau;
        tau = -(tau - sqrt((n + m)*(n - m + 1.0))*legpol[n])/2.0;
        if (abs(x) > 0.5) {
          pie = pie + sqrt((n + m)*(n - m + 1.0))*legpol[n];
          pie = pie*0.5/x;
        } else
          pie = m*legpol2[n]/sin(theta);
        // .......................................................
        // .Done with finding Tau and Pie!                       .
        // .Now to find the scattering amplitude S1, S2, S3, and S4 .
        // .......................................................
        ifact = (-_i)**n;
        if (m = 0) {
          fact1 = d_TE[0, n]*pie + c_TE[0, n]*tau;
          fact2 = c_TM[0, n]*pie + d_TM[0, n]*tau;
          fact3 = c_TE[0, n]*pie + d_TE[0, n]*tau;
          fact4 = d_TM[0, n]*pie + c_TM[0, n]*tau;
          s1 = ifact*fact1;
          s2 = ifact*fact2;
          s3 = ifact*fact3;
          s4 = ifact*fact4;
        } else {
          fact1 = d_TE[m][n]*pie + c_TE[m][n]*tau;
          fact1 = fact1*2.0*cos(m*angle_inc_);
          fact2 = c_TM[m][n]*pie + d_TM[m][n]*tau;
          fact2 = fact2*2.0*cos(m*angle_inc_);
          fact3 = c_TE[m][n]*pie + d_TE[m][n]*tau;
          fact3 = fact3*2.0*_i*sin(m*angle_inc_);
          fact4 = d_TM[m][n]*pie + c_TM[m][n]*tau;
          fact4 = fact4*2.0*_i*sin(m*angle_inc_);
          s1 = ifact*fact1;
          s2 = ifact*fact2;
          s3 = ifact*fact3;
          s4 = ifact*fact4;
        }
        s1_TE[i] = s1_TE[i] + s1;
        s2_TM[i] = s2_TM[i] + s2;
        s3_TE[i] = s3_TE[i] + s3;
        s4_TM[i] = s4_TM[i] + s4;
      }
      // ..............................................
      // . Finished summing over n, now need to       .
      // . sum over m in order to have the complete   .
      // . solution for the Scat Amp Matrix           .
      // ..............................................
    }
    s1_TE[i] = s1_TE[i];
    s2_TM[i] = -_i*s2_TM[i];
    s3_TE[i] = -_i*s3_TE[i];
    s4_TM[i] = s4_TM[i];
    // ..............................................
    // . Done summing over m. Now we have one       .
    // . complete set of scat amp S[i]              .
    // ..............................................
    // ...............................................................
    // .  Lastly, we can now use the scattering amplitude matrix     .
    // .  and solve for the Mueller Matrix elements.  Once we have   .
    // .  the MM, we have all the informations we need.              .
    // ...............................................................
    for (int jj = 1; jj <= 4; jj++) {
      mm[1][jj] = 0.0;
      mm[2][jj] = 0.0;
      mm[3][jj] = 0.0;
      mm[4][jj] = 0.0;
    }
    fact1 = s1_TE[i]*conjg(s1_TE[i])/2.0;
    fact2 = s2_TM[i]*conjg(s2_TM[i])/2.0;
    fact3 = s3_TE[i]*conjg(s3_TE[i])/2.0;
    fact4 = s4_TM[i]*conjg(s4_TM[i])/2.0;
    c = s2_TM[i]*conjg(s3_TE[i]);
    cc = s1_TE[i]*conjg(s4_TM[i]);
    mm[1][1] = (fact1 + fact2 + fact3 + fact4).r;
    mm[1][2] = (fact2 - fact1 - fact3 + fact4).r/mm[1][1];
    mm[1][3] = (c + cc).r/mm[1][1];
    mm[1][4] = (c - cc).i/mm[1][1];
    mm[2][1] = (fact2 - fact1 + fact3 - fact4).r;
    mm[2][2] = (fact2 + fact1 - fact3 - fact4).r/mm[1][1];
    mm[2][3] = (c - cc).r/mm[1][1];
    mm[2][4] = (c + cc).i/mm[1][1];
    fact1 = s2_TM[i]*conjg(s4_TM[i]);
    fact2 = s1_TE[i]*conjg(s3_TE[i]);
    fact3 = s1_TE[i]*conjg(s2_TM[i]);
    fact4 = s3_TE[i]*conjg(s4_TM[i]);
    c = s2_TM[i]*conjg(s1_TE[i]);
    cc = s4_TM[i]*conjg(s3_TE[i]);
    mm[3][1] = (fact1 + fact2).r;
    mm[3][2] = (fact1 - fact2).r/mm[1][1];
    mm[3][3] = (fact3 + fact4).r/mm[1][1];
    mm[3][4] = (c + cc).i/mm[1][1];
    fact1 = s4_TM[i]*conjg(s2_TM[i]);
    fact3 = s1_TE[i]*conjg(s2_TM[i]);
    mm[4][1] = (fact1 + fact2).i;
    mm[4][2] = (fact1 - fact2).i/mm[1][1];
    mm[4][3] = (fact3 - fact4).i/mm[1][1];
    mm[4][4] = (fact3 - fact4).r/mm[1][1];
  }
  S11 = mm[1][1];
  S12 = mm[1][2];
  S21 = mm[3][3];
  S22 = mm[3][4];
  // .........................................................
  // .  Done with calculating the amplitude coefficients     .
  // .........................................................
*/
    isMieCalculated_ = true;
  }


  //**********************************************************************************//
  // This function calculates the expansion coefficients inside the particle,         //
  // required to calculate the near-field parameters.                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it.             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   aln, bln, cln, dln: Complex scattering amplitudes inside the particle          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void EccentricMie::ExpanCoeffs() {
    if (!isScaCoeffsCalc_)
      throw invalid_argument("(ExpanCoeffs) You should calculate external coefficients first!");

    isExpCoeffsCalc_ = false;

/*    complex<double> c_one(1.0, 0.0), c_zero(0.0, 0.0);

    const int L = refractive_index_.size();

    aln_.resize(L + 1);
    bln_.resize(L + 1);
    cln_.resize(L + 1);
    dln_.resize(L + 1);
    for (int l = 0; l <= L; l++) {
      aln_[l].resize(nmax_);
      bln_[l].resize(nmax_);
      cln_[l].resize(nmax_);
      dln_[l].resize(nmax_);
    }

    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int n = 0; n < nmax_; n++) {
      aln_[L][n] = an_[n];
      bln_[L][n] = bn_[n];
      cln_[L][n] = c_one;
      dln_[L][n] = c_one;

      printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", L, n, real(aln_[L][n]), imag(aln_[L][n]), L, n, real(bln_[L][n]), imag(bln_[L][n]), L, n, real(cln_[L][n]), imag(cln_[L][n]), L, n, real(dln_[L][n]), imag(dln_[L][n]));
    }

    vector<complex<double> > D1z(nmax_ + 1), D1z1(nmax_ + 1), D3z(nmax_ + 1), D3z1(nmax_ + 1);
    vector<complex<double> > Psiz(nmax_ + 1), Psiz1(nmax_ + 1), Zetaz(nmax_ + 1), Zetaz1(nmax_ + 1);
    complex<double> denomZeta, denomPsi, T1, T2, T3, T4;

    auto& m = refractive_index_;
    vector< complex<double> > refractive_index_host_(L);

    for (int l = 0; l < L - 1; l++) refractive_index_host_[l] = m[l + 1];
    refractive_index_host_[L - 1] = complex<double> (1.0, 0.0);

    complex<double> z, z1;
    for (int l = L - 1; l >= 0; l--) {
      z = size_param_[l]*m[l];
      z1 = size_param_[l]*refractive_index_host_[l];

      calcD1D3(z, D1z, D3z);
      calcD1D3(z1, D1z1, D3z1);
      calcPsiZeta(z, Psiz, Zetaz);
      calcPsiZeta(z1, Psiz1, Zetaz1);

      for (int n = 0; n < nmax_; n++) {
        int n1 = n + 1;

        denomZeta = refractive_index_host_[l]*Zetaz[n1]*(D1z[n1] - D3z[n1]);
        denomPsi  =  refractive_index_host_[l]*Psiz[n1]*(D1z[n1] - D3z[n1]);

        T1 = aln_[l + 1][n]*Zetaz1[n1] - dln_[l + 1][n]*Psiz1[n1];
        T2 = bln_[l + 1][n]*Zetaz1[n1] - cln_[l + 1][n]*Psiz1[n1];

        T3 = D1z1[n1]*dln_[l + 1][n]*Psiz1[n1] - D3z1[n1]*aln_[l + 1][n]*Zetaz1[n1];
        T4 = D1z1[n1]*cln_[l + 1][n]*Psiz1[n1] - D3z1[n1]*bln_[l + 1][n]*Zetaz1[n1];

        // aln
        aln_[l][n] = (D1z[n1]*refractive_index_host_[l]*T1 + m[l]*T3)/denomZeta;
        // bln
        bln_[l][n] = (D1z[n1]*m[l]*T2 + refractive_index_host_[l]*T4)/denomZeta;
        // cln
        cln_[l][n] = (D3z[n1]*m[l]*T2 + refractive_index_host_[l]*T4)/denomPsi;
        // dln
        dln_[l][n] = (D3z[n1]*refractive_index_host_[l]*T1 + m[l]*T3)/denomPsi;

        printf("aln_[%02i, %02i] = %g,%g; bln_[%02i, %02i] = %g,%g; cln_[%02i, %02i] = %g,%g; dln_[%02i, %02i] = %g,%g\n", l, n, real(aln_[l][n]), imag(aln_[l][n]), l, n, real(bln_[l][n]), imag(bln_[l][n]), l, n, real(cln_[l][n]), imag(cln_[l][n]), l, n, real(dln_[l][n]), imag(dln_[l][n]));
      }  // end of all n
    }  // end of all l

    // Check the result and change  aln_[0][n] and aln_[0][n] for exact zero
    for (int n = 0; n < nmax_; ++n) {
//      printf("n=%d, aln_=%g,%g,   bln_=%g,%g \n", n, real(aln_[0][n]), imag(aln_[0][n]),
//	     real(bln_[0][n]), imag(bln_[0][n]));
      if (abs(aln_[0][n]) < 1e-10) aln_[0][n] = 0.0;
      else throw invalid_argument("Unstable calculation of aln_[0][n]!");
      if (abs(bln_[0][n]) < 1e-10) bln_[0][n] = 0.0;
      else throw invalid_argument("Unstable calculation of bln_[0][n]!");
    }*/

    isExpCoeffsCalc_ = true;
  }  // end of   void EccentricMie::ExpanCoeffs()


  //**********************************************************************************//
  // given nxn matrix a with physical dim NP, this routine replaces it                //
  // by the LU decomposition of a rowwise permutation of itself.  A and               //
  // n are input. A is output. INDX is an output vector which records the             //
  // row permutation effected by the partial pivoting; D is output                    //
  // as +-1 depending on whether the number of row interchanges                       //
  // was even or odd.  This routine is used in combination with Lubksb                //
  // to solve linear equations or invert a matrix.                                    //
  //**********************************************************************************//
  void EccentricMie::LUDecomp(const int n, const int np, vector<vector<complex<double> > >& a,
                              vector<int>& indx, int d) {

    complex<double> suma, cdum, aamax;
    vector<complex<double> > vv(nmax_);
    complex<double> mytiny(1.0e-80,1.0e-80);
    int imax;
    d = 1.0;

    for (int i = 0; i < n; i++) {
      aamax = complex<double>(0.0, 0.0);
      for (int j = 0; j < n; j++) {
        if (abs(a[i][j]) > abs(aamax)) aamax = a[i][j];
      }
      if (abs(aamax) == 0) throw invalid_argument("Singular matrix!!!");
      vv[i] = 1.0/aamax;
    }

    for (int j = 0; j < n; j++) {
      for (int i = 0; i < j - 1; i++) {
        suma = a[i][j];
        for (int k = 0; k < i - 1; k++)
          suma -= a[i][k]*a[k][j];

        a[i][j] = suma;
      }

      aamax = complex<double>(0.0, 0.0);
      for (int i = j - 1; i < n; i++) {
        suma = a[i][j];
        for (int k = 0; k < j - 1; k++)
          suma -= a[i][k]*a[k][j];

          a[i][j] = suma;
          cdum = vv[i]*suma;
          if (abs(cdum) >= abs(aamax)) {
            imax = i;
            aamax = cdum;
          }
      }

      if (j != imax) {
        for (int k = 0; k < n; k++) {
          cdum = a[imax][k];
          a[imax][k] = a[j][k];
          a[j][k] = cdum;
        }
        d = -d;
        vv[imax] = vv[j];
      }

      indx[j] = imax;
      if(abs(a[j][j]) == 0) a[j][j] = mytiny;
      if (j != n) {
        cdum = 1.0/a[j][j];
        for (int i = j; i < n; i++)
          a[i][j] = a[i][j]*cdum;
      }
    }
  }  // end of void EccentricMie::LUDecomp(...)


  //**********************************************************************************//
  // Solves the set of n linear equations A.X=B. Here A is input, not                 //
  // as the matrix A but rather as its LU decomp, from LUDCMP.  INDX                  //
  // is input as the permutation vector returned by LUDCMP. B is input                //
  // as the right-hand side vector B, and returns with the solution                   //
  // Vector X.  A,N,NP and INDX are not modified by this routine and                  //
  // can be left in place for successive calls with different right-hand              //
  // sides B.  This routine takes into account the possibility that B                 //
  // can begin with many zero elements, so it is efficient for use in                 //
  // matrix inversion.                                                                //
  //**********************************************************************************//
  void EccentricMie::LUSolve(const int n, const int np, vector<vector<complex<double> > >& a,
                             vector<int>& indx, vector<complex<double> > b) {

    complex<double> suma;
    int ii = 0, ll;

    for (int i = 0; i < n; i++) {
      ll = indx[i];
      suma = b[ll];
      b[ll] = b[i];

      if (ii != 0) {
        for (int j = ii - 1; j< i - 1; j++)
          suma -= a[i][j]*b[j];
      } else {
        if (abs(suma) != 0.0)  ii = i;
      }
      b[i] = suma;
    }

    for (int i = n - 1; i >= 0; i--) {
      suma = b[i];
      if (i < n) {
        for (int j = i; j< n; j++)
          suma -= a[i][j]*b[j];
      }
      b[i] = suma/a[i][i];
    }
  }  // end of void EccentricMie::LUSolve(...)


  //**********************************************************************************//
  // This function calculates the electric (E) and magnetic (H) fields inside and     //
  // around the particle.                                                             //
  //                                                                                  //
  // Input parameters (coordinates of the point):                                     //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic fields                                     //
  //**********************************************************************************//
  void EccentricMie::calcField(const double Rho, const double Theta, const double Phi,
                                vector<complex<double> >& E, vector<complex<double> >& H)  {

    /*complex<double> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    vector<complex<double> > ipow = {c_one, c_i, -c_one, -c_i} // Vector containing precomputed integer powers of i to avoid computation
    vector<complex<double> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    vector<complex<double> > M1o1n(3), M1e1n(3), N1o1n(3), N1e1n(3);
    vector<complex<double> > jn(nmax_ + 1), jnp(nmax_ + 1), h1n(nmax_ + 1), h1np(nmax_ + 1);
    vector<double> Pi(nmax_), Tau(nmax_);

    vector<complex<double> > Ei(3), Hi(3);

    int l = 0;  // Layer number
    complex<double> ml;

    // Initialize E and H
    for (int i = 0; i < 3; i++) {
      E[i] = c_zero;
      H[i] = c_zero;
    }

    if (Rho > size_param_.back()) {
      l = size_param_.size();
      ml = c_one;
    } else {
      for (int i = size_param_.size() - 1; i >= 0 ; i--) {
        if (Rho <= size_param_[i]) {
          l = i;
        }
      }
      ml = refractive_index_[l];
    }

    // Calculate spherical Bessel and Hankel functions and their derivatives
    sbesjh(Rho*ml, jn, jnp, h1n, h1np);

    // Calculate angular functions Pi and Tau
    calcPiTau(cos(Theta), Pi, Tau);*/
/*    printf("Thetd = %g, cos(Theta) = %g\n", Theta, cos(Theta));
    printf("Pi:\n");
    for (auto p : Pi) printf("%11.4e\n",p);
    printf("Tau:\n");
    for (auto p : Tau) printf("%11.4e\n",p);
*/
    /*for (int n = nmax_ - 2; n >= 0; n--) {
      int n1 = n + 1;
      double rn = static_cast<double>(n1);

      // using BH 4.12 and 4.50
      calcSpherHarm(Rho, Theta, Phi, jn[n1], jnp[n1], Pi[n], Tau[n], rn, M1o1n, M1e1n, N1o1n, N1e1n);
      calcSpherHarm(Rho, Theta, Phi, h1n[n1], h1np[n1], Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);

      // Total field in the lth layer: eqs. (1) and (2) in Yang, Appl. Opt., 42 (2003) 1710-1720
      complex<double> En = ipow[n1 % 4]*(rn + rn + 1.0)/(rn*rn + rn);
      for (int i = 0; i < 3; i++) {
        // electric field E [V m - 1] = EF*E0
        E[i] += En*(cln_[l][n]*M1o1n[i] - c_i*dln_[l][n]*N1e1n[i]
              + c_i*aln_[l][n]*N3e1n[i] -     bln_[l][n]*M3o1n[i]);

        H[i] += En*(-dln_[l][n]*M1e1n[i] - c_i*cln_[l][n]*N1o1n[i]
              +  c_i*bln_[l][n]*N3o1n[i] +     aln_[l][n]*M3e1n[i]);
        Ei[i] += En*(M1o1n[i] - c_i*N1e1n[i]);
        Hi[i] += En*(-M1e1n[i] - c_i*N1o1n[i]);

      }
    }  // end of for all n

    //printf("rho = %11.4e; angle_inc_ = %11.4eº; theta = %11.4eº; x[%i] = %11.4e; m[%i] = %11.4er%+10.5ei\n", Rho, Phi*180./PI_, Theta*180./PI_, l, size_param_[l], l, real(ml), imag(ml));
    // magnetic field
    double hffact = 1.0/(cc_*mu_);
    for (int i = 0; i < 3; i++) {
      H[i] = hffact*H[i];
      Hi[i] *= hffact;
      printf("E[%i] = %10.5er%+10.5ei; H[%i] = %10.5er%+10.5ei\n", i, real(E[i]), imag(E[i]), i, real(H[i]), imag(H[i]));
    }*/
   }  // end of EccentricMie::calcField(...)


  //**********************************************************************************//
  // This function calculates complex electric and magnetic field in the surroundings //
  // and inside the particle.                                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to 0 (zero) and the function will calculate it.       //
  //   ncoord: Number of coordinate points                                            //
  //   Coords: Array containing all coordinates where the complex electric and        //
  //           magnetic fields will be calculated                                     //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic field at the provided coordinates          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  void EccentricMie::RunFieldCalculation() {
    double Rho, Theta, Phi;

    // Calculate scattering coefficients an_ and bn_
    ScattCoeffs();

    // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
    ExpanCoeffs();
    // for (int i = 0; i < nmax_; ++i) {
    //   printf("cln_[%i] = %11.4er%+10.5ei;  dln_[%i] = %11.4er%+10.5ei\n", i, real(cln_[0][i]), imag(cln_[0][i]),
    // 	     i, real(dln_[0][i]), imag(dln_[0][i]));
    // }


    long total_points = coords_[0].size();
    E_.resize(total_points);
    H_.resize(total_points);
    for (auto& f : E_) f.resize(3);
    for (auto& f : H_) f.resize(3);

    for (int point = 0; point < total_points; point++) {
      const double& Xp = coords_[0][point];
      const double& Yp = coords_[1][point];
      const double& Zp = coords_[2][point];

      // Convert to spherical coordinates
      Rho = sqrt(pow2(Xp) + pow2(Yp) + pow2(Zp));

      // If Rho=0 Theta is undefined. Just set it to zero to avoid problems
      Theta = (Rho > 0.0) ? acos(Zp/Rho) : 0.0;

      // If Xp=Yp=0 Phi is undefined. Just set it to zero to avoid problems
      if (Xp == 0.0)
        Phi = (Yp != 0.0) ? asin(Yp/sqrt(pow2(Xp) + pow2(Yp))) : 0.0;
      else
        Phi = acos(Xp/sqrt(pow2(Xp) + pow2(Yp)));

      // Avoid convergence problems due to Rho too small
      if (Rho < 1e-5) Rho = 1e-5;

      //printf("X = %g; Y = %g; Z = %g; pho = %g; angle_inc_ = %g; theta = %g\n", Xp, Yp, Zp, Rho, Phi*180./PI_, Theta*180./PI_);

      //*******************************************************//
      // external scattering field = incident + scattered      //
      // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
      // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
      //*******************************************************//

      // This array contains the fields in spherical coordinates
      vector<complex<double> > Es(3), Hs(3);

      // Firstly the easiest case: the field outside the particle
      //      if (Rho >= GetSizeParameter()) {
      //        fieldExt(Rho, Theta, Phi, Es, Hs);
      // } else {
      calcField(Rho, Theta, Phi, Es, Hs);  //Should work fine both: inside and outside the particle
      //}
 { //Now, convert the fields back to cartesian coordinates
        E_[point][0] = sin(Theta)*cos(Phi)*Es[0] + cos(Theta)*cos(Phi)*Es[1] - sin(Phi)*Es[2];
        E_[point][1] = sin(Theta)*sin(Phi)*Es[0] + cos(Theta)*sin(Phi)*Es[1] + cos(Phi)*Es[2];
        E_[point][2] = cos(Theta)*Es[0] - sin(Theta)*Es[1];

        H_[point][0] = sin(Theta)*cos(Phi)*Hs[0] + cos(Theta)*cos(Phi)*Hs[1] - sin(Phi)*Hs[2];
        H_[point][1] = sin(Theta)*sin(Phi)*Hs[0] + cos(Theta)*sin(Phi)*Hs[1] + cos(Phi)*Hs[2];
        H_[point][2] = cos(Theta)*Hs[0] - sin(Theta)*Hs[1];
      }
    }  // end of for all field coordinates
  }  //  end of EccentricMie::RunFieldCalculation()
}  // end of namespace eccmie

