#ifndef VISOUCSTERMS_HPP
#define VISOUCSTERMS_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Viscous forces
//
//----------------------------------------------------------------------------------//

#include "variables.hpp"
#include "vector.hpp"
#include "generalTools.hpp"

//----------------------------------------------------------------------------------//

struct VISCOSITY_NONE{

  static double ViscousTerm(Vector3<double> dr,
                            Vector3<double> dv,
                            double Irho,
                            double Jrho,
                            double drs,
                            ConstantVariables &C)
  {return 0;}

};

//----------------------------------------------------------------------------------//

struct VISCOSITY_Artificial{

  static double ViscousTerm(Vector3<double> dr,
                            Vector3<double> dv,
                            double Irho,
                            double Jrho,
                            double drs,
                            ConstantVariables &C)
  {
    double visco_mu = C.h*dot(dr, dv)/(drs*drs + C.eps*C.h*C.h);
    double visco = (dot(dr, dv) < 0) ? (-C.eta*C.c0*visco_mu/(0.5*(Irho+Jrho))) : (0.);

    return visco;
  }

};

//----------------------------------------------------------------------------------//

struct VISCOSITY_ArtificialPlus{

  static double ViscousTerm(Vector3<double> dr,
                            Vector3<double> dv,
                            double Irho,
                            double Jrho,
                            double drs,
                            ConstantVariables &C)
  {
    double eta = MAX(dot(dr, dv), C.eta);
    double visco_mu = C.h*dot(dr, dv)/(drs*drs + C.eps*C.h*C.h);
    double visco = (dot(dr, dv) < 0) ? (-eta*C.c0*visco_mu/(0.5*(Irho+Jrho))) : (0.);

    return visco;
  }

};


//----------------------------------------------------------------------------------//

struct VISCOSITY_Physical{

  static double ViscousTerm(Vector3<double> dr,
                            Vector3<double> dv,
                            double Irho,
                            double Jrho,
                            double drs,
                            ConstantVariables &C)
  {
    double visco;
    return visco;
  }

};

//----------------------------------------------------------------------------------//

#endif
