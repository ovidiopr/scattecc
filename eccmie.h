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

#define VERSION "0.3.1"
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace eccmie {
  int ScattCoeffs(double& xi, double& xh, const int pl, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const int nmax, std::vector<std::complex<double> > &an, std::vector<std::complex<double> > &bn);
  int eccMie(double& xi, double& xh, const int pl, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int eccMie(double& xi, double& xh, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int eccMie(double& xi, double& xh, const int pl, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int eccMie(double& xi, double& xh, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int eccField(double& xi, double& xh, const int pl, std::complex<double>& mi, std::complex<double>& mh, double& dx, double& alpha, const int nmax, const unsigned int ncoord, const std::vector<double>& Xp_vec, const std::vector<double>& Yp_vec, const std::vector<double>& Zp_vec, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H);

  class EccentricMie {
   public:
    // Run calculation
    void RunMieCalculation();
    void RunFieldCalculation();

    // Return calculation results
    double GetQext();
    double GetQsca();
    double GetQabs();
    double GetQbk();
    double GetQpr();
    double GetAsymmetryFactor();
    double GetAlbedo();
    std::vector<std::complex<double> > GetS1();
    std::vector<std::complex<double> > GetS2();

    std::vector<std::complex<double> > GetAn(){return an_;};
    std::vector<std::complex<double> > GetBn(){return bn_;};

    // Problem definition
    // Modify size of inclusion
    void SetIncSize(const double inc_size);
    // Modify size of host
    void SetHostSize(const double host_size);
    // Modify inclusion's shift
    void SetIncShift(const double inc_shift);
    // Modify angle of incidence
    void SetIncAngle(const double inc_angle);
    // Modify refractive index of inclusion
    void SetIncIndex(const std::complex<double> inc_index);
    // Modify refractive index of host
    void SetHostIndex(const std::complex<double> host_index);
    // Modify scattering (theta) angles
    void SetAngles(const std::vector<double>& angles);
    // Modify coordinates for field calculation
    void SetFieldCoords(const std::vector< std::vector<double> >& coords);
    // Modify PEC layer
    void SetPECLayer(int layer_position = -1);

    // Set a fixed value for the maximun number of terms
    void SetMaxTerms(int nmax);
    // Get maximun number of terms
    int GetMaxTerms() {return nmax_;};

    // Clear layer information
    void ClearLayers();

    // Applied units requests
    double GetSizeParameter();
    double GetLayerWidth(int layer_position = 0);
    std::vector<double> GetLayersSize();
    std::vector<std::complex<double> > GetLayersIndex();
    std::vector<std::array<double, 3> > GetFieldCoords();

    std::vector<std::vector< std::complex<double> > > GetFieldE(){return E_;};   // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<double> > > GetFieldH(){return H_;};
  private:
    void calcNstop();
    void calcNmax(unsigned int first_layer);

    std::complex<double> calc_an(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
                                 std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                 std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1);
    std::complex<double> calc_bn(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
                                 std::complex<double> PsiXL, std::complex<double> ZetaXL,
                                 std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1);

    std::complex<double> calc_S1(int n, std::complex<double> an, std::complex<double> bn,
                                 double Pi, double Tau);
    std::complex<double> calc_S2(int n, std::complex<double> an, std::complex<double> bn,
                                 double Pi, double Tau);
    void calcD1D3(std::complex<double> z,
                  std::vector<std::complex<double> >& D1,
                  std::vector<std::complex<double> >& D3);
    void calcPsiZeta(std::complex<double> x,
                     std::vector<std::complex<double> >& Psi,
                     std::vector<std::complex<double> >& Zeta);
    void calcPiTau(const double& costheta,
                   std::vector<double>& Pi, std::vector<double>& Tau);
    void calcSpherHarm(const double Rho, const double Theta, const double Phi,
                       const std::complex<double>& zn, const std::complex<double>& dzn,
                       const double& Pi, const double& Tau, const double& n,
                       std::vector<std::complex<double> >& Mo1n, std::vector<std::complex<double> >& Me1n, 
                       std::vector<std::complex<double> >& No1n, std::vector<std::complex<double> >& Ne1n);
    void ScattCoeffs();
    void ExpanCoeffs();

    void LUDecomp(const int n, const int np, std::vector<std::vector<std::complex<double> > >& a,
                  std::vector<std::complex<double> >& indx, int d);
    void LUSolve(const int n, const int np, std::vector<std::vector<std::complex<double> > >& a,
                 std::vector<std::complex<double> >& indx, std::vector<std::complex<double> > b);

    void calcField(const double Rho, const double Theta, const double Phi,
                   std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H);

    bool isExpCoeffsCalc_ = false;
    bool isScaCoeffsCalc_ = false;
    bool isMieCalculated_ = false;

    // Size parameter for inclusion and host, inclusion's shift and angle of incidence
    double size_param_inc_, size_param_host_, shift_inc_, angle_inc_;
    // Refractive index for inclusion and host
    std::complex<double> refractive_index_inc_, refractive_index_host_;
    // Scattering angles for scattering pattern in radians
    std::vector<double> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;

    // with calcNmax(int first_layer);
    int nmax_ = -1;
    int nmax_preset_ = -1;
    // Scattering coefficients
    std::vector<std::complex<double> > an_, bn_;
    std::vector< std::vector<double> > coords_;
    // TODO: check if l index is reversed will lead to performance
    // boost, if $a^(L+1)_n$ stored in aln_[n][0], $a^(L)_n$ in
    // aln_[n][1] and so on...
    // at the moment order is forward!
    std::vector< std::vector<std::complex<double> > > aln_, bln_, cln_, dln_;
    /// Store result
    double Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    std::vector<std::vector< std::complex<double> > > E_, H_;  // {X[], Y[], Z[]}
    std::vector<std::complex<double> > S1_, S2_;

    //Used constants
    const double PI_=3.14159265358979323846;
    // light speed [m s-1]
    double const cc_ = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    double const mu_ = 4.0*PI_*1.0e-7;

    //Temporary variables
    std::vector<std::complex<double> > PsiZeta_;


  };  // end of class EccentricMie

}  // end of namespace eccmie
