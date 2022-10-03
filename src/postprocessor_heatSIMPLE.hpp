#ifndef INTERPOLATIONHANDLER_HPP
#define INTERPOLATIONHANDLER_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Interaction handler
//
//-----------------------------------------------------------------------------------//

#include <CompactNSearch>

//-----------------------------------------------------------------------------------//

#include "vector.hpp"
#include "variables.hpp"
#include "interpolationStructures.hpp"
#include "kernel.hpp"

//-----------------------------------------------------------------------------------//

template<
typename INTERPOLATION,
typename KERNEL
>
class PostProcessingHandler{
public:

  void Interpolate(Variables<double, double> &HEAT_SIMPLEnodes,
                   Variables<double, double> &HEAT_SIMPLEsolid,
                   Variables<double, double> &HEAT_SIMPLEbound,
                   ConstantVariables HEAT_SIMPLEconstants);

  PostProcessingHandler(double h,
                        Variables<double, double> &HEAT_SIMPLEnodes,
                        Variables<double, double> &HEAT_SIMPLEsolid,
                        Variables<double, double> &HEAT_SIMPLEbound);

  ~PostProcessingHandler(){};

private:
  CompactNSearch::NeighborhoodSearch pp_nsearch;
  unsigned int positions_id;
  unsigned int positions_id_wall;

  unsigned int positions_id_nodes;

};

//-----------------------------------------------------------------------------------//

template<
typename INTERPOLATION,
typename KERNEL
>
PostProcessingHandler<
INTERPOLATION,
KERNEL>::PostProcessingHandler(double h,
                               Variables<double, double> &HEAT_SIMPLEnodes,
                               Variables<double, double> &HEAT_SIMPLEsolid,
                               Variables<double, double> &HEAT_SIMPLEbound)
                               : pp_nsearch(h, false)
{

  positions_id_nodes = pp_nsearch.add_point_set(&HEAT_SIMPLEnodes.r.front().x,
                                                HEAT_SIMPLEnodes.N);

  positions_id = pp_nsearch.add_point_set(&HEAT_SIMPLEsolid.r.front().x,
                                          HEAT_SIMPLEsolid.N);
  positions_id_wall = pp_nsearch.add_point_set(&HEAT_SIMPLEbound.r.front().x,
                                               HEAT_SIMPLEbound.N);

  pp_nsearch.set_active(positions_id_nodes, positions_id, true);
  pp_nsearch.set_active(positions_id_nodes, positions_id_wall, true);

}

//-----------------------------------------------------------------------------------//

template<
typename INTERPOLATION,
typename KERNEL
>
void PostProcessingHandler<
INTERPOLATION,
KERNEL>::Interpolate(Variables<double, double> &HEAT_SIMPLEnodes,
                     Variables<double, double> &HEAT_SIMPLEsolid,
                     Variables<double, double> &HEAT_SIMPLEbound,
                     ConstantVariables HEAT_SIMPLEconstants)
{


  pp_nsearch.find_neighbors();

  // Update fluid particles
  CompactNSearch::PointSet const& ps_1 = pp_nsearch.point_set(positions_id_nodes);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_1.n_points(); ++i)
  {
    InterpolationData<double, double> ptcI;
    InterpolationData<double, double> ptcJ;
    ptcI.CopyNodeDataIn(HEAT_SIMPLEnodes, i);

    //Get point set 1 neighbors of point set 1.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(HEAT_SIMPLEsolid, pid);
      INTERPOLATION::template SolidSolidInterpolation<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
    }

    /*
    //Get point set 1 neighbors of point set 2.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
      ptcJ.CopyBoundaryDataIn(HEAT_SIMPLEbound, pid);
      INTERPOLATION::template SolidSolidInterpolation<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
    }
    */

    ptcI.CopyDataOut(HEAT_SIMPLEnodes, i);
  }

}

//-----------------------------------------------------------------------------------//

#endif
