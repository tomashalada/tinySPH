#ifndef WRITEPARTICLEDATA_HPP
#define WRITEPARTICLEDATA_HPP

//----------------------------------------------------------------------------------//
//
// tinySPH: Write data to .ptcs
//
//----------------------------------------------------------------------------------//

#include "vector.hpp"
#include "variables.hpp"

//----------------------------------------------------------------------------------//

//template<typename F> //works
//void writeParticleData(F field, std::string fileName){ //work
template<typename T, typename S>
void writeParticleData(const Variables<T, S> &field, std::string fileName){

  std::ofstream outputFile;
  outputFile.open(fileName);

  for(unsigned int p = 0; p < field.N; p++)
    outputFile << \
    field.r[p].x << " " << \
    field.r[p].y << " " << \
    field.r[p].z << " " << \
    field.v[p].x << " " << \
    field.v[p].y << " " << \
    field.v[p].z << " " << \
    field.rho[p] << " " << \
    field.p[p] <<	std::endl;

}

//----------------------------------------------------------------------------------//

#endif
