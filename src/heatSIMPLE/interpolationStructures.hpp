#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: heatSIMPLE - interaction structures
//
//----------------------------------------------------------------------------------//

#include <tuple>
#include "vector.hpp"
#include "variables.hpp"

#include "kernel.hpp"

//----------------------------------------------------------------------------------//


template<typename K, typename S>
struct InterpolationData
{

  //Data required for the interaction
  Vector3<double> r = {0., 0., 0.};
  double T = 0.;
  double rho = 0.;

  double gamma = 0;
  double Wsum = 0;

  void CopyDataIn(Variables<K,S> &vars,  unsigned int i);
  void CopyNodeDataIn(Variables<K,S> &vars,  unsigned int i);
  void CopyDataOut(Variables<K,S> &vars,  unsigned int i);

};

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InterpolationData<K, S>::CopyDataIn(Variables<K, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> rho = vars.rho[i];
  this -> T = vars.T[i];

}

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InterpolationData<K, S>::CopyNodeDataIn(Variables<K, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];

}

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InterpolationData<K, S>::CopyDataOut(Variables<K, S> &vars,  unsigned int i)
{

  if(gamma > 0.5){
    vars.T[i] = T/gamma;
    vars.rho[i] = rho/gamma;

    //vars.T[i] = p/Wsum;
    //vars.rho[i] = rho/Wsum;
  }
  else{
    vars.T[i] = 0.;
    vars.rho[i] = 0.;
  }

}

//----------------------------------------------------------------------------------//

struct heatSIMPLE_INTERPOLATION{

template<
typename K,
typename S,
typename KERNEL
>
static void SolidSolidInterpolation(InterpolationData<K, S> &I,
                                    InterpolationData<K, S> &J,
                                    ConstantVariables &C)
{

  double h = C.h; double m = C.m;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  double drs = dr.length();

  double W = KERNEL::W(drs, h);
  double V = m/J.rho;

  I.gamma += W*V;

  I.T += J.T*W*V;
  I.rho += W*V;

  //I.Wsum += W;
  //I.T += J.T*W;
  //I.rho += J.rho*W;

}

};

//-----------------------------------------------------------------------------------//

#endif
