#ifndef VARIABLES_HPP
#define VARIABLES_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: heatSIMPLE - variables
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
  double rho0;
  double lambda;
  double cp;
  double dp;

  ConstantVariables(double _h, double _m, double _rho0, double _lambda, double _cp, double _dp):
  	h(_h), m(_m), rho0(_rho0), lambda(_lambda), cp(_cp), dp(_dp)
  {

  }

  ~ConstantVariables(){}

};

//-----------------------------------------------------------------------------------//

template<typename K, typename S>
class Variables {
public:

  size_t N;

  std::vector<unsigned int> t;

  std::vector<Vector3d> r;

  std::vector<double> rho;
  std::vector<double> T;
  std::vector<double> dT;

  void allocate(size_t N);
  void initializeWithGeometryFile(std::string geoFileName);

  Variables(){}
  ~Variables(){}

};

//-----------------------------------------------------------------------------------//

template<typename K, typename S>
class SolidVars: public Variables<K, S> {

};

//-----------------------------------------------------------------------------------//

template<typename K, typename S>
class BoundVars: public Variables<K, S> {

};

//-----------------------------------------------------------------------------------//

template <typename K, typename S>
void Variables<K, S>::allocate(size_t n){

  this-> N = n;

  this-> t.resize(n);
  this-> r.resize(n);
  this-> rho.resize(n);
  this-> T.resize(n); this-> dT.resize(n);

}

//-----------------------------------------------------------------------------------//

template <typename K, typename S>
void Variables<K, S>::initializeWithGeometryFile(std::string geoFileName){

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

    this-> rho[line-1] = std::stod(tokens[3]);
    this-> T[line-1] = std::stod(tokens[4]);

    line++;
    if(line == N + 1){break;}
  }

}

//-----------------------------------------------------------------------------------//

#endif
