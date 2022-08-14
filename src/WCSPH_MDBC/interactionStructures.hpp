#ifndef INTERACTION_HPP
#define INTERACTION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Interaction - pure WCSPH
//
//----------------------------------------------------------------------------------//

#include <tuple>
#include "vector.hpp"
#include "matrix.hpp"
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
  Vector3<double> gn = {0., 0., 0.};

  //Results of the interaction
  double drho = 0;
  Vector3<double> a = {0., 0., 0.};

  double gamma = 0; //BI

  void CopyDataIn(Variables<T,S> &vars,  unsigned int i);
  void CopyBoundaryDataIn(BoundVars<T,S> &vars,  unsigned int i);

  void CopyDataOut(Variables<T,S> &vars,  unsigned int i);
  void CopyBoundaryDataOut(Variables<T,S> &vars,  unsigned int i);

};

//----------------------------------------------------------------------------------//

template<typename T, typename S>
struct InteractionDataGhostNode : InteractionData<T, S>{

  Matrix3<double> A = {0., 0., 0.,
                       0., 0., 0.,
                       0., 0., 0.};

  Vector3<double> b = {0., 0., 0.};

  Vector3<double> v = {0., 0., 0.};

  void CopyBoundaryDataOut(BoundVars<T,S> &vars,  unsigned int i);

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
void InteractionData<T, S>::CopyBoundaryDataIn(BoundVars<T, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];

  this -> gn = vars.gn[i];

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

}

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InteractionDataGhostNode<T, S>::CopyBoundaryDataOut(BoundVars<T, S> &vars, unsigned int i)
{

  if (A.det() > 0.001 ){
    Vector3<double> gradRho = solve(A, b);
    Vector3<double> dr_bgn = vars.r[i] - vars.gn[i];
    vars.rho[i] = gradRho.x + gradRho.y*dr_bgn.x + gradRho.z*dr_bgn.z;
    vars.v[i] = (-1)*v/A.a11;
  }
  else if (A.a11 > 0.){
    vars.rho[i] = b.x/A.a11;
    vars.v[i] = (-1)*v/A.a11;
  }
  else{
    vars.rho[i] = 1000.;
    vars.v[i] = {0., 0., 0.};
  }

  /* Can we do something with this, please? */
  if (vars.rho[i] < 1000.){
    vars.rho[i] = 1000.;
  }

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

struct WCSPH_MDBC{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void UpdateBoundaryData(InteractionDataGhostNode<T, S> &I,
                               InteractionData<T, S> &J,
                               ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  Matrix3<double> A = {0., 0., 0.,
                       0., 0., 0.,
                       0., 0., 0.};

  //Resolv kernel function values:
  Vector3<double> dr = J.r - I.gn; // POS_FLUID - POS_GN
  double drs = dr.length();

  std::pair<double, double> WF = KERNEL::WF(drs, h);
  double W = WF.first;
  Vector3<double> gradW = dr*WF.second;

  //continuuity equation
  double V = m/J.rho;

  A.a11 = W*V;       A.a12 = dr.x*W*V;       A.a13 = dr.z*W*V;
  A.a21 = gradW.x*V; A.a22 = dr.x*gradW.x*V; A.a23 = dr.z*gradW.x*V;
  A.a31 = gradW.z*V; A.a32 = dr.x*gradW.z*V; A.a33 = dr.z*gradW.z*V;

  I.A += A;

  I.b.x += W*m; I.b.y += gradW.x*m; I.b.z += gradW.z*m;

}

};

//-----------------------------------------------------------------------------------//

struct WCSPH_MDBCvelocity{

template<
typename T,
typename S,
typename KERNEL,
typename DIFF_TERM,
typename VISCO_TERM
>
static void UpdateBoundaryData(InteractionDataGhostNode<T, S> &I,
                               InteractionData<T, S> &J,
                               ConstantVariables &C)
{

  double h = C.h; double m = C.m; double rho0 = C.rho0;
  double c0 = C.c0; double eta = C.eta; double eps = 0.001;
  double delta = C.delta;

  Matrix3<double> A = {0., 0., 0.,
                       0., 0., 0.,
                       0., 0., 0.};

  //Resolv kernel function values:
  Vector3<double> dr = J.r - I.gn; // POS_FLUID - POS_GN
  double drs = dr.length();

  std::pair<double, double> WF = KERNEL::WF(drs, h);
  double W = WF.first;
  Vector3<double> gradW = dr*WF.second;

  //continuuity equation
  double V = m/J.rho;

  A.a11 = W*V;       A.a12 = dr.x*W*V;       A.a13 = dr.z*W*V;
  A.a21 = gradW.x*V; A.a22 = dr.x*gradW.x*V; A.a23 = dr.z*gradW.x*V;
  A.a31 = gradW.z*V; A.a32 = dr.x*gradW.z*V; A.a33 = dr.z*gradW.z*V;

  I.A += A;

  I.b.x += W*m; I.b.y += gradW.x*m; I.b.z += gradW.z*m;

  I.v += J.v*W*V;

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
