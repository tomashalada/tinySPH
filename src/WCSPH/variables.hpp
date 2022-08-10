#ifndef VARIABLES_HPP
#define VARIABLES_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Boundary integrals formulation variables
//
//-----------------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <iterator>

#include "vector.hpp"

//-----------------------------------------------------------------------------------//

class ConstantVariables{
public:

  double h;
  double m;
  double c0;
  double eta;
  double rho0;
  double eps;
  double delta = 0.1;
  double cb;
  double dp;

  ConstantVariables(double _h, double _m, double _c0, double _eta, double _rho0, double _dp):
  	h(_h), m(_m), c0(_c0), eta(_eta), rho0(_rho0), dp(_dp)
  {
  	cb =  _c0 * _c0 * _rho0 / 7.;
  }

  ~ConstantVariables(){}

};

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
class Variables {
public:

  size_t N;

  std::vector<unsigned int> t;

  std::vector<Vector3d> r;
  std::vector<Vector3d> v;
  std::vector<Vector3d> a;

  std::vector<double> rho;
  std::vector<double> drho;

  std::vector<double> p;

  void allocate(size_t N);
  void initializeWithGeometryFile(std::string geoFileName);

  Variables(){}
  ~Variables(){}

};

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
class FluidVars: public Variables<T, S> {

};

//-----------------------------------------------------------------------------------//

template<typename T, typename S>
class BoundVars: public Variables<T, S> {

};

//-----------------------------------------------------------------------------------//

template <typename T, typename S>
void Variables<T, S>::allocate(size_t n){

  this-> N = n;

  this-> t.resize(n);
  this-> r.resize(n); this-> v.resize(n); this-> a.resize(n);
  this-> rho.resize(n); this-> drho.resize(n);
  this-> p.resize(n);

}

//-----------------------------------------------------------------------------------//

template <typename T, typename S>
void Variables<T, S>::initializeWithGeometryFile(std::string geoFileName){

  std::ifstream geoFile(geoFileName);
  std::string geoFileText;

  unsigned int line = 0;

  while (getline (geoFile, geoFileText)) {

    if(line == 0){
      this-> N = std::stoi(geoFileText);
      this-> Variables::allocate(N);
      line++;
      continue;
    }

    std::istringstream iss(geoFileText);
    std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};

    this-> r[line-1].x = std::stod(tokens[0]);
    this-> r[line-1].y = std::stod(tokens[1]);
    this-> r[line-1].z = std::stod(tokens[2]);

    this-> v[line-1].x = std::stod(tokens[3]);
    this-> v[line-1].y = std::stod(tokens[4]);
    this-> v[line-1].z = std::stod(tokens[5]);

    this-> rho[line-1] = std::stod(tokens[6]);
    this-> p[line-1] = std::stod(tokens[7]);

    line++;
    if(line == N + 1){break;}
  }

}

//-----------------------------------------------------------------------------------//

#endif
