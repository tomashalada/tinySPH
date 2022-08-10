#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Interaction - pure WCSPH
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

  Integrator(Variables<double, double> *V, double _dt) : MODEL(V), dt(_dt) {
    this -> N = V -> N;
  }


};

//----------------------------------------------------------------------------------//

class VerletScheme : public Integrator {
public:

  std::vector<Vector3<double>> v_o;
  std::vector<Vector3<double>> v_oo;

  std::vector<double> rho_o;
  std::vector<double> rho_oo;

  double test;

  void ComputeStep();
  void ComputeStepEulerForStability();

  VerletScheme(Variables<double, double> *V, double dt) : Integrator(V, dt) {

    this -> v_o.resize(N);
    this -> v_oo.resize(N);

    this -> rho_o.resize(N);
    this -> rho_oo.resize(N);

    for(int p = 0; p < N; p++)
    {

      v_oo[p] = MODEL->v[p];; //swap pointers only instead!
      v_o[p] = MODEL->v[p];

      rho_oo[p] = MODEL->rho[p];; //swap pointers only instead!
      rho_o[p] = MODEL->rho[p];

    }
  }

};

//----------------------------------------------------------------------------------//

void VerletScheme::ComputeStep(){

  #pragma omp parallel for schedule(static)
  for(int p = 0; p < N; p++)
  {

    MODEL->v[p] = v_oo[p] + MODEL->a[p]*dt*2.;
    MODEL->r[p] = MODEL->r[p] + v_o[p]*dt + MODEL->a[p]*dt*dt*0.5;
    MODEL->rho[p] = rho_oo[p] + MODEL->drho[p]*dt*2;

    v_oo[p] = v_o[p]; //swap pointers only instead!
    v_o[p] = MODEL->v[p];

    rho_oo[p] = rho_o[p]; //swap pointers only instead!
    rho_o[p] = MODEL->rho[p];

  }

}

//----------------------------------------------------------------------------------//

void VerletScheme::ComputeStepEulerForStability(){

  #pragma omp parallel for schedule(static)
  for(int p = 0; p < N; p++)
  {

    MODEL->v[p] = MODEL->v[p] + MODEL->a[p]*dt;
    MODEL->r[p] = MODEL->r[p] + v_o[p]*dt + MODEL->a[p]*dt*dt*0.5;
    MODEL->rho[p] = MODEL->rho[p] + MODEL->drho[p]*dt;

    v_oo[p] = v_o[p]; //swap pointers only instead!
    v_o[p] = MODEL->v[p];

    rho_oo[p] = rho_o[p]; //swap pointers only instead!
    rho_o[p] = MODEL->rho[p];

  }

}

//----------------------------------------------------------------------------------//

class SymplecticScheme : public Integrator {
public:

  std::vector<Vector3<double>> r_o;
  std::vector<Vector3<double>> v_o;
  std::vector<double> rho_o;

  double test;

  void ComputePredictor();
  void ComputeCorrector();

  SymplecticScheme(Variables<double, double> *V, double dt) : Integrator(V, dt) {

    this -> r_o.resize(N);
    this -> v_o.resize(N);
    this -> rho_o.resize(N);

    for(int p = 0; p < N; p++)
    {

      r_o[p] = MODEL->r[p];
      v_o[p] = MODEL->v[p];
      rho_o[p] = MODEL->rho[p];

    }
  }

};

//----------------------------------------------------------------------------------//

void SymplecticScheme::ComputePredictor(){

  #pragma omp parallel for schedule(static)
  for(int p = 0; p < N; p++)
  {

    r_o[p] = MODEL->r[p];
    v_o[p] = MODEL->v[p];
    rho_o[p] = MODEL->rho[p];

    MODEL->r[p] = MODEL->r[p] + v_o[p]*dt*0.5;
    MODEL->v[p] = MODEL->v[p] + MODEL->a[p]*dt*0.5;
    MODEL->rho[p] = MODEL->rho[p] + MODEL->drho[p]*dt*0.5;

  }

}

//----------------------------------------------------------------------------------//

void SymplecticScheme::ComputeCorrector(){

  #pragma omp parallel for schedule(static)
  for(int p = 0; p < N; p++)
  {

    double epsilon = -(MODEL->drho[p]/MODEL->rho[p])*dt;
    MODEL->rho[p] = rho_o[p] * ((2. - epsilon) / (2. + epsilon));

    MODEL->v[p] = v_o[p] + MODEL->a[p]*dt;
    MODEL->r[p] = r_o[p] + (MODEL->v[p] + v_o[p])*dt*0.5;

  }

}

//----------------------------------------------------------------------------------//

#endif
