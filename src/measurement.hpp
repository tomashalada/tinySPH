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
class MEASUREMENT_WCSPH: public MEASUREMENT
{
public:

  unsigned int numberOfSensors;
  float dt;

  std::vector< std::vector< float > > storage_p;
  std::vector< std::vector< float > > storage_rho;
  std::vector< std::vector< Vector3d > > storage_v;

  VARIABLES sensors;

  MEASUREMENT_WCSPH( unsigned int _numerOfSensors, float _dt)
  : numberOfSensors( _numerOfSensors ), dt( _dt )
  {

    this-> storage_p.resize(_numerOfSensors);
    this-> storage_rho.resize(_numerOfSensors);
    this-> storage_v.resize(_numerOfSensors);

  }
  ~MEASUREMENT_WCSPH(  ){  };

  void StoreSensorsData( bool save_pressure = true,
                         bool save_density = true,
                         bool save_velocity = true)
  {

    for(int i = 0; i < numberOfSensors; i++)
    {
      if(save_pressure)
        storage_p[i].push_back(sensors.p[i]);

      //if(save_density)
      //  storage_p[i].push_back(sensors.p[i]);

      //if(save_velocity)
      //  storage_p[i].push_back(sensors.p[i]);
    }
  }

  void WritePressureToFile( std::string fileName )
  {
    std::ofstream file;
    file.open( fileName );

    for( int i = 0; i < storage_p[0].size(); i++ )
    {
      file << dt*i;
      for( int j = 0; j < sensors.N; j++ )
      {
        file << " " << storage_p[ j ][ i ];
      }
     file << "\n";
    }

    file.close();
  }

};

//-----------------------------------------------------------------------------------//

#endif
