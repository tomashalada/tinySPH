#ifndef DIFFUSIVETEMRS_HPP
#define DIFFUSIVETEMRS_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Diffusive terms
//
//-----------------------------------------------------------------------------------//

#include "variables.hpp"

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

  static double Psi()
  {
    double D;
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
