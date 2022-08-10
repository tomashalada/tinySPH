#ifndef GENERALTOOLS_HPP
#define GENERALTOOLS_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: general tools
//
//-----------------------------------------------------------------------------------//

#include <string>
#include "vector.hpp"
#include "matrix.hpp"

//-----------------------------------------------------------------------------------//

#define POS(x,y,X,Y) (x*Y + y)


std::string coutVector2(Vector3<double> X)
{
  std::string vectString;
  vectString = "[ " + std::to_string(X.x) + "," + std::to_string(X.z) + "]";
  return vectString;
}

std::string coutVector3(Vector3<double> X)
{
  std::string vectString;
  vectString = "[ " + std::to_string(X.x) + "," + std::to_string(X.y) + "," + std::to_string(X.z) + "]";
  return vectString;
}

//-----------------------------------------------------------------------------------//

std::string coutMatrix3(Matrix3<double> A)
{
  std::string matString;
  matString = "[ " + std::to_string(A.a11) + "," + std::to_string(A.a12) + "," + std::to_string(A.a13) + "\n" \
                   + std::to_string(A.a21) + "," + std::to_string(A.a22) + "," + std::to_string(A.a23) + "\n" \
                   + std::to_string(A.a31) + "," + std::to_string(A.a32) + "," + std::to_string(A.a33) + " ]";
  return matString;
}

//-----------------------------------------------------------------------------------//

std::string stepToName(std::string name, unsigned int step)
{
  std::string nameWithStep;
  nameWithStep = name + "_" + std::to_string(step);
  return nameWithStep;
}

std::string stepToNameWithPtcsExtension(std::string name, unsigned int step)
{
  std::string nameWithStep;
  nameWithStep = name + "_" + std::to_string(step) + ".ptcs";
  return nameWithStep;
}

std::string addToName(std::string name, std::string add)
{
  std::string nameString;
  nameString = name + "_" + add;
  return nameString;
}

//-----------------------------------------------------------------------------------//

#endif

