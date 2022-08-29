#ifndef INTERACTION_HPP
#define INTERACTION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Interaction - pure WCSPH
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

  //Results of the interaction
  double drho = 0;
  Vector3<double> a = {0., 0., 0.};

  Vector3<double> pw_vecTerm = {0., 0., 0.};//wall pressure
  double pw = 0; //wall pressure
  double gamma = 0; //renormalisation

  void CopyDataIn(Variables<T,S> &vars,  unsigned int i);
  void CopyBoundaryDataIn(Variables<T,S> &vars,  unsigned int i);

  void CopyDataOut(Variables<T,S> &vars,  unsigned int i);
  void CopyBoundaryDataOut(Variables<T,S> &vars,  unsigned int i);

};

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyDataIn(Variables<T, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];

}

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyBoundaryDataIn(Variables<T, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];

}


//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyDataOut(Variables<T, S> &vars,  unsigned int i)
{

  vars.drho[i] = drho;
  vars.a[i] = a;

}

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionData<T, S>::CopyBoundaryDataOut(Variables<T, S> &vars,  unsigned int i)
{

  Vector3<double> g = {0., 0., -9.81};

  //vars.p[i] = (pw + dot((g-a), pw_vecTerm))/gamma;
  if(gamma > 0)
  {
    vars.p[i] = pw/gamma;
    //vars.v[i] = (-1)*v/gamma;
  }
  else
  {
    vars.p[i] = WCSPH_EOS_Tait(1000., 34.3, 1000.);
    //vars.v[i] = {0., 0., 0.};
  }
  //vars.v[i] = (-1)*v/gamma;
  vars.a[i] = {0., 0., 0.};

}

//----------------------------------------------------------------------------------//

struct WCSPH_FLUIDFLUID{

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
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += dot(dv, gradW)*m - diffTerm;

  //momentum equation
  double p_term = (I.p + J.p)/(I.rho * J.rho);
  double visco = VISCO_TERM::ViscousTerm(dr, dv, I.rho, J.rho, drs, C);

  I.a -= (p_term + visco)*gradW*m;

}

};

//-----------------------------------------------------------------------------------//

struct WCSPH_FLUIDBOUND_DBC{

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
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += dot(dv, gradW)*m - diffTerm;

  //momentum equation
  double p_term = (I.p + J.p)/(I.rho * J.rho);
  double visco = VISCO_TERM::ViscousTerm(dr, dv, I.rho, J.rho, drs, C);

  I.a -= (p_term + visco)*gradW*m;

}

};

//-----------------------------------------------------------------------------------//

struct WCSPH_DBC{

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
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //continuuity equation
  double psi = DIFF_TERM::Psi(I, J, C, drs);
  double diffTerm = psi*dot(dr, gradW)*m/J.rho;

  I.drho += dot(dv, gradW)*m - diffTerm;
  I.a = {0., 0., 0.}; //zmiznout

}

};

//-----------------------------------------------------------------------------------//

struct WCSPH_GWBC{

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
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;
  Vector3<double> g = {0., 0., -9.81};


  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  Vector3<double> dv = I.v - J.v;
  double drs = dr.length();

  double W = KERNEL::W(drs, h);
  double F = KERNEL::F(drs, h);
  Vector3<double> gradW = dr*F;

  //I.pw += J.p*W;
  //I.pw_vecTerm += J.rho*dr*W;
  //I.pw_vecTerm += dr*W;

  I.pw += J.p*W + dot(g, dr)*J.rho*W;

  double p_term = (I.p + J.p)/(I.rho * J.rho);
  double visco = VISCO_TERM::ViscousTerm(dr, dv, I.rho, J.rho, drs, C);

  //I.v += J.v*W;
  //I.a -= (p_term + visco)*gradW*m;
  I.a = {0.,0.,0.};


  I.gamma += W;
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
