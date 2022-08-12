#ifndef DIFFUSIVETEMRS_HPP
#define DIFFUSIVETEMRS_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Diffusive terms
//
//-----------------------------------------------------------------------------------//

#include "variables.hpp"
#include "pressureAndEOS.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------//

struct DT_NONE{

  static double Psi(
  InteractionData<double, double> &I,
  InteractionData<double, double> &J,
  ConstantVariables &C,
  double drs)
  {return 0;}

};

//-----------------------------------------------------------------------------------//

struct DT_Molteni{

  static double Psi(
  InteractionData<double, double> &I,
  InteractionData<double, double> &J,
  ConstantVariables &C,
  double drs)
  {

    double D;
    D = 2. * C.c0 *C.h * C.delta * (J.rho - I.rho)/(drs*drs);
    return D;
  }

};

//-----------------------------------------------------------------------------------//

struct DT_Fourtakas{

  static double Psi(
  InteractionData<double, double> &I,
  InteractionData<double, double> &J,
  ConstantVariables &C,
  double drs)
  {

    double D;
    double dph = 1. + (I.r.z - J.r.z)* C.rho0 * 9.81 / C.cb;
    double drhoH = C.rho0 * pow(dph, 1./7) - C.rho0;
    D = 2. * C.c0 *C.h * C.delta * ((J.rho - I.rho) - drhoH)/(drs*drs);
    return D;

  }

};

//-----------------------------------------------------------------------------------//

struct DT_Antuono{

  static double Psi()
  {
    double D;
    return D;
  }

};

//-----------------------------------------------------------------------------------//

#endif
