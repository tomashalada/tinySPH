#ifndef PRESSUREANDEOS_HPP
#define PRESSUREANDEOS_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Pressure - equations of state
//
//-----------------------------------------------------------------------------------//

#include <cmath>

//-----------------------------------------------------------------------------------//

double WCSPH_EOS_Tait(double rho, double c0, double rho0)
{

  double p;

  const double gamma = 7.; //Weakly compressible state eq. Poisson constant
  const double b_const = c0*c0*rho0/gamma; //Weakly compressible state eq. constant
  const double rho_r = rho/rho0; //Weakly compressible state eq. relative density

  //Weakly compressible Tait state equation
  p = b_const*( pow(rho_r, gamma) - 1 );

  return p;

}

//-----------------------------------------------------------------------------------//

double WCSPH_EOS_TaitLinearized(double rho, double c0, double rho0)
{

  double p;

  //Linearized weakly compressible Tait state equation
  p = c0*c0*( rho - rho0 );

  return p;

}

//-----------------------------------------------------------------------------------//

void DensityToPressure(Variables<double, double> &V, ConstantVariables &C)
{

  const double c0 = C.c0;
  const double rho0 = C.rho0;
  const unsigned int N = V.N;


  for(unsigned int i = 0; i < N; i++)
  {
    V.p[i] = WCSPH_EOS_Tait(V.rho[i], c0, rho0);
    //V.p[i] = WCSPH_EOS_TaitLinearized(V.rho[i], c0, rho0);
  }

}

//-----------------------------------------------------------------------------------//

#endif
