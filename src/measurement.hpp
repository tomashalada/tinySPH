#ifndef MEASUREMENT_HPP
#define MEASUREMENT_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Measurement template
//
//-----------------------------------------------------------------------------------//

#include "vector.hpp"
#include "variables.hpp"

//-----------------------------------------------------------------------------------//

class MEASUREMENT
{
  public:

};

//-----------------------------------------------------------------------------------//

template< typename VARIABLES = Variables< double, double > >
class MEASUREMENT_pressure: public MEASUREMENT
{
public:

  float dt;
  std::vector< float > p1, p2, p3, p4;
  std::vector< Vector3d > positions;

  // This si way to go.
  //  std::vector< float > storage_p;
  //  std::vector< float > storage_rho;
  //  std::vector< Vector3d > storage_v;

  VARIABLES sensors;

  MEASUREMENT_pressure( std::vector< Vector3d > _positions, float _dt)
  : positions( _positions ), dt( _dt )
  {

  }
  ~MEASUREMENT_pressure(  ){  };

  void StoreSensorsData( )
  {
    p1.push_back(sensors.p[0]);
    p2.push_back(sensors.p[1]);
    p3.push_back(sensors.p[2]);
    p4.push_back(sensors.p[3]);
  }

  void WritePressureToFile( std::string fileName )
  {
    std::ofstream file;
    file.open( fileName );

    for( int i = 0; i < p1.size(); i++ )
      //for( int j = 0; j < sensors.N; j++ )
        file << dt*i << " " << p1[ i ] << " " << p2[ i ] << " " << p3[ i ] << " " << p4[ i ] << std::endl;

    file.close();
  }


};

//-----------------------------------------------------------------------------------//

#endif
