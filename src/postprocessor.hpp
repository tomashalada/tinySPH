#ifndef POSTPROCESSINGHANDLER_HPP
#define POSTPROCESSINGHANDLER_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: PostProcessor
//
//-----------------------------------------------------------------------------------//

#include <CompactNSearch>

//-----------------------------------------------------------------------------------//

#include "vector.hpp"
#include "variables.hpp"
#include "interpolationStructures.hpp"
#include "kernel.hpp"
#include "diffusiveTerms.hpp"
#include "generalTools.hpp"

//-----------------------------------------------------------------------------------//


template<
typename INTERPOLATION,
typename KERNEL
>
class PostProcessingHandler{
public:

  void Interpolate(Variables<double, double> &WCSPHnodes,
                   Variables<double, double> &WCSPHfluid,
                   Variables<double, double> &WCSPHwall,
                   ConstantVariables WCSPHconstants);

  void AddInterpolationPlane(Variables<double, double> &WCSPHinterplationNodes);

  PostProcessingHandler(double h,
                        Variables<double, double> &WCSPHnodes,
                        Variables<double, double> &WCSPHfluid,
                        Variables<double, double> &WCSPHwall);

  ~PostProcessingHandler(){};

private:
  CompactNSearch::NeighborhoodSearch pp_nsearch;
  unsigned int positions_id;
  unsigned int positions_id_wall;

  unsigned int positions_id_nodes;
  //std::vector<unsinged int> listOfPlanes;

};

//-----------------------------------------------------------------------------------//

template<
typename INTERPOLATION,
typename KERNEL
>
PostProcessingHandler<
INTERPOLATION,
KERNEL>::PostProcessingHandler(double h,
                               Variables<double, double> &WCSPHnodes,
                               Variables<double, double> &WCSPHfluid,
                               Variables<double, double> &WCSPHwall)
                               : pp_nsearch(h, false)
{

  positions_id_nodes = pp_nsearch.add_point_set(&WCSPHnodes.r.front().x, WCSPHnodes.N);
  positions_id = pp_nsearch.add_point_set(&WCSPHfluid.r.front().x, WCSPHfluid.N);
  positions_id_wall = pp_nsearch.add_point_set(&WCSPHwall.r.front().x, WCSPHwall.N);

  //pp_nsearch.set_active(positions_id, positions_id_nodes, false);
  //pp_nsearch.set_active(positions_id, positions_id_wall, false);
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
KERNEL>::AddInterpolationPlane(Variables<double, double> &WCSPHinterplationNodes)
{

  positions_id_nodes = pp_nsearch.add_point_set(&WCSPHinterplationNodes.r.front().x,
                                                WCSPHinterplationNodes.N);

  pp_nsearch.set_active(positions_id_nodes, positions_id, true);

}


//-----------------------------------------------------------------------------------//

template<
typename INTERPOLATION,
typename KERNEL
>
void PostProcessingHandler<
INTERPOLATION,
KERNEL>::Interpolate(Variables<double, double> &WCSPHnodes,
                     Variables<double, double> &WCSPHfluid,
                     Variables<double, double> &WCSPHbound,
                     ConstantVariables WCSPHconstants)
{

  //Compute pressure from density
  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  pp_nsearch.find_neighbors();

  // Update fluid particles
  CompactNSearch::PointSet const& ps_1 = pp_nsearch.point_set(positions_id_nodes);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_1.n_points(); ++i)
  {
    InterpolationData<double, double> ptcI;
    InterpolationData<double, double> ptcJ;
    ptcI.CopyNodeDataIn(WCSPHnodes, i);

    //Get point set 1 neighbors of point set 1.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(WCSPHfluid, pid);
      INTERPOLATION::template FluidFluidInterpolation<double, double, KERNEL>(ptcI, ptcJ, WCSPHconstants);
    }

    /*
    //Get point set 1 neighbors of point set 2.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
      ptcJ.CopyDataIn(WCSPHbound, pid);
      INTERPOLATION::template FluidFluidInterpolation<double, double, KERNEL>(ptcI, ptcJ, WCSPHconstants);
    }
    std::cout << std::endl;
    */

    ptcI.CopyDataOut(WCSPHnodes, i);
  }


}

//-----------------------------------------------------------------------------------//

#endif
