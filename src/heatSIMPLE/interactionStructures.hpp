#ifndef INTERACTION_HPP
#define INTERACTION_HPP

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
struct InteractionData
{

  //Data required for the interaction
  Vector3<double> r = {0., 0., 0.};
  double T = 0.;
  double rho = 0.;

  //Results of the interaction
  double dT = 0;

  double gamma = 0; //BI

  void CopyDataIn(Variables<K,S> &vars,  unsigned int i);
  void CopyBoundaryDataIn(Variables<K,S> &vars,  unsigned int i);

  void CopyDataOut(Variables<K,S> &vars,  unsigned int i);
  void CopyBoundaryDataOut(Variables<K,S> &vars,  unsigned int i);

};

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InteractionData<K, S>::CopyDataIn(Variables<K, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> rho = vars.rho[i];
  this -> T = vars.T[i];

}

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InteractionData<K, S>::CopyBoundaryDataIn(Variables<K, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> rho = vars.rho[i];
  this -> T = vars.T[i];

}


//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InteractionData<K, S>::CopyDataOut(Variables<K, S> &vars,  unsigned int i)
{

  vars.dT[i] = dT;

}

//----------------------------------------------------------------------------------//

template<typename K, typename S>
void InteractionData<K, S>::CopyBoundaryDataOut(Variables<K, S> &vars,  unsigned int i)
{

}

//----------------------------------------------------------------------------------//

struct heatSIMPLE_SOLIDSOLID{

template<
typename K,
typename S,
typename KERNEL
>
static void SolidSolidInteraction(InteractionData<K, S> &I,
                                  InteractionData<K, S> &J,
                                  ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double lambda = C.lambda; double cp = C.cp; double dp = C.dp;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  double drs = dr.length();

  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //I.dT += (1/cp)*(2*lambda/(I.rho*J.rho))*(I.T - J.T)*dot(dr, gradW)*m/(drs*drs);
  //I.dT += 1./J.rho*(I.T - J.T)*dot(dr, gradW)*m/(drs*drs);
  I.dT += 2./(I.rho*J.rho)*(I.T - J.T)*dot(dr, gradW)*m/(drs*drs);
  //I.dT += (I.T - J.T)*dot(dr, gradW)/(drs) * (dp*dp);

}

};

//-----------------------------------------------------------------------------------//

struct heatSIMPLE_SOLIDBOUND{

template<
typename K,
typename S,
typename KERNEL
>
static void SolidBoundInteraction(InteractionData<K, S> &I,
                                  InteractionData<K, S> &J,
                                  ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double lambda = C.lambda; double cp = C.cp; double dp = C.dp;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  double drs = dr.length();

  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //I.dT += (1/cp)*(2*lambda/(I.rho*J.rho))*(I.T - J.T)*dot(dr, gradW)*m/(drs*drs);
  I.dT += 2./(I.rho*J.rho)*(I.T - J.T)*dot(dr, gradW)*m/(drs*drs);
  //I.dT += (I.T - J.T)*dot(dr, gradW)/(drs) * (dp*dp);

}

};

//-----------------------------------------------------------------------------------//

struct heatSIMPLE_BC{

template<
typename K,
typename S,
typename KERNEL
>
static void UpdateBoundaryData(InteractionData<K, S> &I,
                               InteractionData<K, S> &J,
                               ConstantVariables &C)
{

}

};

//-----------------------------------------------------------------------------------//

#endif
