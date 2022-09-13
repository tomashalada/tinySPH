//-----------------------------------------------------------------------------------//
//
// tinySPH: DamBreak-measurement: measure pressure
//
//-----------------------------------------------------------------------------------//

#include "measurement.hpp"

//-----------------------------------------------------------------------------------//

class MEASUREMENT_pressure
{
public:

  float dt;
  std::vector< float > p1, p2, p3, p4;
  std::vector< Vector3d > positions;

  MEASUREMENT_pressure( std::vector< Vector3d > _positions, float dt)
  : positions( _positions )
  {

  }

  void WritePressureToFile( std::string fileName )
  {
    std::ofstream file;
    file.open( fileName );

    for( int i = 0; i < p1.size(); i++ )
      file << dt*i << " " << p1[ i ] << p2[ i ] << p3[ i ] << p4[ i ] << std::endl;

    file.close();
  }

  ~MEASUREMENT_pressure();

};

//-----------------------------------------------------------------------------------//
