#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: heatSIMPLE - integrators
//
//----------------------------------------------------------------------------------//

#include "variables.hpp"

//----------------------------------------------------------------------------------//

class Integrator{
public:

  unsigned int N;
  double dt;
  Variables<double, double> *MODEL;

  void Assign_dt(double _dt){ dt = _dt; }

  Integrator(Variables<double, double> *V, double _dt) : MODEL(V), dt(_dt)
  {
    this -> N = V -> N;
  }

};

//----------------------------------------------------------------------------------//

class EulerScheme : public Integrator {
public:

  void ComputeStep();

  EulerScheme(Variables<double, double> *V, double dt) : Integrator(V, dt)
  {

  }
  ~EulerScheme(){}

};

//----------------------------------------------------------------------------------//

void EulerScheme::ComputeStep(){

  #pragma omp parallel for schedule(static)
  for(int p = 0; p < N; p++)
  {
    MODEL->T[p] = MODEL->T[p] + MODEL->dT[p]*dt;
  }

}

//----------------------------------------------------------------------------------//

#endif
