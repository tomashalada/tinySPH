#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Interpolation - pure WCSPH
//
//----------------------------------------------------------------------------------//

#include <tuple>

//----------------------------------------------------------------------------------//

template<typename T, typename S>
struct InterpolationData
{

  //Data required for the interaction
  Vector3<double> r = {0., 0., 0.};
  Vector3<double> v = {0., 0., 0.};
  double p = 0.;
  double rho = 0.;

  double gamma = 0;
  double Wsum = 0;

  void CopyDataIn(Variables<T,S> &vars,  unsigned int i);
  void CopyNodeDataIn(Variables<T,S> &vars,  unsigned int i);
  void CopyDataOut(Variables<T,S> &vars,  unsigned int i);

};

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InterpolationData<T, S>::CopyDataIn(Variables<T, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];
  this -> v = vars.v[i];
  this -> p = vars.p[i];
  this -> rho = vars.rho[i];

}

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InterpolationData<T, S>::CopyNodeDataIn(Variables<T, S> &vars,  unsigned int i)
{

  this -> r = vars.r[i];

}

//----------------------------------------------------------------------------------//

template<typename T, typename S>
void InterpolationData<T, S>::CopyDataOut(Variables<T, S> &vars,  unsigned int i)
{

  if(gamma > 0.5){
    vars.rho[i] = rho/gamma;
    vars.p[i] = p/gamma;
    vars.v[i] = v/gamma;

    //vars.rho[i] = rho/Wsum;
    //vars.p[i] = p/Wsum;
    //vars.v[i] = v/Wsum;
  }
  else{
    vars.rho[i] = 0.;
    vars.p[i] = 0.;
    vars.v[i] = {0., 0., 0.};
  }

}

//----------------------------------------------------------------------------------//


struct WCSPH_INTERPOLATION{

template<
typename T,
typename S,
typename KERNEL
>
static void FluidFluidInterpolation(InterpolationData<T, S> &I,
                                    InterpolationData<T, S> &J,
                                    ConstantVariables &C)
{

  double h = C.h; double m = C.m;

  //Resolv kernel function values:
  Vector3<double> dr = I.r - J.r;
  double drs = dr.length();

  double W = KERNEL::W(drs, h);
  double V = m/J.rho;

  I.gamma += W*V;

  I.rho += W*m;
  I.p += J.p*W*V;
  I.v += J.v*W*V;

  //I.Wsum += W;
  //I.rho += J.rho*W;
  //I.p += J.p*W;
  //I.v += J.v*W;

}

};

//-----------------------------------------------------------------------------------//

#endif
