#ifndef INTERACTION_HPP
#define INTERACTION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Interaction - WCSPH with boundary integrals
//
//----------------------------------------------------------------------------------//

#include <tuple>
#include "vector.hpp"
#include "variables.hpp"

#include "kernel.hpp"
#include "pressureAndEOS.hpp"

//----------------------------------------------------------------------------------//

template<typename T, typename S>
struct InteractionData
{

  //Data required for the interaction
  Vector3<double> r = {0., 0., 0.};
  Vector3<double> v = {0., 0., 0.};
  double p = 0.;
  double rho = 0.;

  //Boundary data
  Vector3<double> n = {0., 0., 0.};

  //Results of the interaction
  double drho = 0;
  Vector3<double> a = {0., 0., 0.};

  double gamma = 0; //BI

  void CopyDataIn(Variables<T,S> &vars, unsigned int i);
  void CopyBoundaryDataIn(BoundVars<T,S> &vars, unsigned int i);

  void CopyDataOut(Variables<T,S> &vars, unsigned int i);
  void CopyBoundaryDataOut(Variables<T,S> &vars, unsigned int i);

};

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyDataIn(Variables<T, S> &vars, unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];

}

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyBoundaryDataIn(BoundVars<T, S> &vars, unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];
  this -> n = vars.n[i];

}

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyDataOut(Variables<T, S> &vars, unsigned int i)
{

  vars.drho[i] = drho/gamma;
  vars.a[i] = a/gamma;

}

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyBoundaryDataOut(Variables<T, S> &vars, unsigned int i)
{

  if(gamma > 0.0001){
    vars.rho[i] = (drho/gamma > 1000) ? (drho/gamma) : (1000.);
  }
  else
    vars.rho[i] = 1000.;

}

//-----------------------------------------------------------------------------------//

struct WCSPHBI_FLUIDFLUID{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void FluidFluidInteraction(InteractionData<T, S> &I,
                                  InteractionData<T, S> &J,
                                  ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = 0.;
  double delta = C.delta;

  //Resolve kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  std::pair<double, double> WF = KERNEL::WF(drs, h); //first W, second F
  Vector3<double> gradW = dr*WF.second;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += dot(dv, gradW)*m - diffTerm;

  //momentum equation
  double p_term = (I.p + J.p)/(I.rho * J.rho);
  double visco = VISCO_TERM::ViscousTerm(dr, dv, I.rho, J.rho, drs, C);

  I.a -= (p_term + visco)*gradW*m;

  //Shepard renormalisation - computed only for fluid neighbours
  I.gamma += WF.first*m/J.rho;

}

};

//-----------------------------------------------------------------------------------//

struct WCSPHBI_FLUIDBOUND_BI{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void FluidBoundInteraction(InteractionData<T, S> &I,
                                  InteractionData<T, S> &J,
                                  ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = 0.;
  double delta = C.delta; double dp = C.dp;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  std::pair<double, double> WF = KERNEL::WF(drs, h); //first W, second F
  Vector3<double> gradW = dr*WF.second;
  double W = WF.first;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += (-1)*dot(dv, J.n)*W*J.rho*dp; //no DT for boundary: - diffTerm;

  //momentum equation
  double p_term = (I.p + J.p)/(I.rho * J.rho);
  double visco = VISCO_TERM::ViscousTerm(dr, dv, I.rho, J.rho, drs, C);

  I.a -= (-1)*(p_term + visco)*W*J.n*J.rho*dp;

  //Shepard renormalisation - computed only for fluid neighbours
  //I.gamma += W*m/J.rho;

}

};

//-----------------------------------------------------------------------------------//

struct WCSPH_BI{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void UpdateBoundaryData(InteractionData<T, S> &I,
                               InteractionData<T, S> &J,
                               ConstantVariables &C)
{

  //missing h, m, rho0, c0, eta,
  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  //double W = KERNEL::W(drs, h);
  std::pair<double, double> WF = KERNEL::WF(drs, h); //first W, second F
  Vector3<double> gradW = dr*WF.second;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += dot(dv, gradW)*m;// - diffTerm;

  //Shepard renormalisation
  I.gamma += WF.first*m/J.rho;

}

};

//-----------------------------------------------------------------------------------//

struct BOUNDARY_SIMPLE{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void UpdateBoundaryData(InteractionData<T, S> &I,
                               InteractionData<T, S> &J,
                               ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = C.eps;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  double W = KERNEL::W(drs, h);

  //continuuity equation
  I.drho += m*W;
  I.gamma += W*(m/J.rho);
  I.a = {0., 0., 0.}; //zmiznout

}

};

//-----------------------------------------------------------------------------------//

#endif
