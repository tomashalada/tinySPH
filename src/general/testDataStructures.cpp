#include <iostream>

#include "point.hpp"
#include "vector.hpp"
#include "matrix.hpp"

#include "generalTools.hpp"

int main(){

  Matrix3<double> A(3, 8, 6, 4, 1, 2, 3.5, 7, 1);
  Matrix3<double> B(16, 2.2, 6, 3.2, 2, 9, 8, 7, 3);

  std::cout << "Matrix A times consant c = 2 ..... : " << std::endl;
  Matrix3<double> Ac2 = A*2;
  std::cout << coutMatrix3(Ac2) << std::endl;

  std::cout << "Matrix A times matrix B .......... : " << std::endl;
  Matrix3<double> AB = A * B;
  Matrix3<double> ABresult(121.6, 64.6, 108., 83.2, 24.8, 39., 86.4, 28.7, 87.0);
  std::cout << coutMatrix3(AB) << std::endl;

  Vector3<double> b(2, 3, 11.6);
  std::cout << coutVector3(b) << std::endl;

  std::cout << "Matrix M times vector b .......... : " << std::endl;
  Vector3<double> Ab = A*b;
  Vector3<double> Abresult(99.6, 34.2, 39.6);
  std::cout << coutVector3(Ab) << std::endl;

  std::cout << "Solve system Mx = b .............. : " << std::endl;
  Vector3<double> x = solve(A, b);
  //Vector3<double> Abresult(99.6, 34.2, 39.6);
  std::cout << coutVector3(x) << std::endl;

  std::cout << "Determinant of matrix A .......... : " << std::endl;
  double detA = A.det();
  double detAresult = 132;
  std::cout << detA << std::endl;


  return EXIT_SUCCESS;
}
