#ifndef INTERACTIONHANDLER_HPP
#define INTERACTIONHANDLER_HPP

//-----------------------------------------------------------------------------------//
//
// tinySPH: Interaction handler
//
//-----------------------------------------------------------------------------------//

#include <CompactNSearch>

//-----------------------------------------------------------------------------------//

#include "vector.hpp"
#include "variables.hpp"
#include "interactionStructures.hpp"
#include "kernel.hpp"
#include "diffusiveTerms.hpp"

//-----------------------------------------------------------------------------------//

#include "interpolationStructures.hpp"
#include "measurement.hpp"

//-----------------------------------------------------------------------------------//


template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
class InteractionHandler{
public:

  void Interact(Variables<double, double> &WCSPHfluid,
                BoundVars<double, double> &WCSPHbound,
                ConstantVariables WCSPHconstants,
                bool updateBoundary);

  void UpdateBoundary(Variables<double, double> &WCSPHfluid,
                      BoundVars<double, double> &WCSPHbound,
                      ConstantVariables WCSPHconstants);

  /* Measurement */
  template<
  typename VARIABLES = Variables< double, double>,
  typename INTERPOLATION
  >
  void Measure(MEASUREMENT_pressure< VARIABLES > &WCSPHpressure,
               Variables<double, double> &WCSPHfluid,
               BoundVars<double, double> &WCSPHbound,
               ConstantVariables WCSPHconstants);

  template <typename VARIABLES = Variables< double, double> >
  void AddSensors(MEASUREMENT_pressure< VARIABLES > &WCSPHpressure);

  InteractionHandler(double h,
                     Variables<double, double> &WCSPHfluid,
                     Variables<double, double> &WCSPHbound);

  ~InteractionHandler(){};

private:
  CompactNSearch::NeighborhoodSearch nsearch;
  unsigned int positions_id;
  unsigned int positions_id_wall;

  /* Measurement */
  unsigned int positions_id_sensors;

};

//-----------------------------------------------------------------------------------//

template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::InteractionHandler(double h,
                            Variables<double, double> &WCSPHfluid,
                            Variables<double, double> &WCSPHbound)
                            : nsearch(h, true)
{

  positions_id = nsearch.add_point_set(&WCSPHfluid.r.front().x, WCSPHfluid.N);
  positions_id_wall = nsearch.add_point_set(&WCSPHbound.r.front().x, WCSPHbound.N);

}

//-----------------------------------------------------------------------------------//

template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
template<
typename VARIABLES
>
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::AddSensors(MEASUREMENT_pressure< VARIABLES > &WCSPHpressure)
{

  positions_id_sensors = nsearch.add_point_set(&WCSPHpressure.positions.front().x,
                                                WCSPHpressure.positions.size());

  nsearch.set_active(positions_id, positions_id_sensors, false);
  nsearch.set_active(positions_id_wall, positions_id_sensors, false);
  nsearch.set_active(positions_id_sensors, positions_id_sensors, false);

}

//-----------------------------------------------------------------------------------//

template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::Interact(Variables<double, double> &WCSPHfluid,
                  BoundVars<double, double> &WCSPHbound,
                  ConstantVariables WCSPHconstants,
                  bool updateBoundary)
{

  //Compute pressure from density
  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  nsearch.find_neighbors();

  // Update fluid particles
  CompactNSearch::PointSet const& ps_1 = nsearch.point_set(positions_id);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_1.n_points(); ++i)
  {
    InteractionData<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyDataIn(WCSPHfluid, i);

    //Get point set 1 neighbors of point set 1.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(WCSPHfluid, pid);
      FLUID_FLUID::template FluidFluidInteraction<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    //Get point set 1 neighbors of point set 2.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
      ptcJ.CopyBoundaryDataIn(WCSPHbound, pid);
      FLUID_BOUND::template FluidBoundInteraction<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    ptcI.CopyDataOut(WCSPHfluid, i);
    WCSPHfluid.a[i].z -= 9.81; //apply external forces
  }

  // Update boundary particles
  if(updateBoundary)
  {
    CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_wall);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ps_2.n_points(); ++i)
    {

      InteractionData<double, double> ptcI;
      InteractionData<double, double> ptcJ;
      ptcI.CopyDataIn(WCSPHbound, i);
      for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
      {
        // Return the point id of the jth neighbor of the ith particle in the point_set_1.
        const unsigned int pid = ps_2.neighbor(positions_id, i, j);
        ptcJ.CopyDataIn(WCSPHfluid, pid);
        BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
      }

      ptcI.CopyBoundaryDataOut(WCSPHbound, i);

    }
  }

}

//-----------------------------------------------------------------------------------//

template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::UpdateBoundary(Variables<double, double> &WCSPHfluid,
                        BoundVars<double, double> &WCSPHbound,
                        ConstantVariables WCSPHconstants)
{

  nsearch.find_neighbors();

  // Update boundary particles
  CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_wall);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_2.n_points(); ++i)
  {

    InteractionData<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyDataIn(WCSPHbound, i);
    for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
    {
      // Return the point id of the jth neighbor of the ith particle in the point_set_1.
      const unsigned int pid = ps_2.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(WCSPHfluid, pid);
      BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    ptcI.CopyBoundaryDataOut(WCSPHbound, i);

  }

}

//-----------------------------------------------------------------------------------//

template<
typename FLUID_FLUID,
typename FLUID_BOUND,
typename BOUND_UPDATE,
typename DIFFUSIVE_TERM,
typename VISOUCS_TERM,
typename KERNEL
>
template<
typename VARIABLES,
typename INTERPOLATION
>
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::Measure(MEASUREMENT_pressure< VARIABLES > &WCSPHpressure,
                 Variables<double, double> &WCSPHfluid,
                 BoundVars<double, double> &WCSPHbound,
                 ConstantVariables WCSPHconstants)
{

   //Compute pressure from density
   DensityToPressure(WCSPHfluid, WCSPHconstants);
   DensityToPressure(WCSPHbound, WCSPHconstants);

   //pp_nsearch.find_neighbors(); //<--- do pointwise
                                  //<--- or find the nbs in main loop

   // Update fluid particles
   CompactNSearch::PointSet const& ps_3 = nsearch.point_set(positions_id_sensors);
   #pragma omp parallel for schedule(static)
   for (int i = 0; i < ps_3.n_points(); ++i)
   {
     InterpolationData<double, double> ptcI;
     InterpolationData<double, double> ptcJ;
     ptcI.CopyNodeDataIn(WCSPHpressure.sensors, i);

     //Get point set 1 neighbors of point set 1.
     for (size_t j = 0; j < ps_3.n_neighbors(positions_id, i); ++j)
     {
       const unsigned int pid = ps_3.neighbor(positions_id, i, j);
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

     ptcI.CopyDataOut(WCSPHpressure.sensors, i);
   }


}

//-----------------------------------------------------------------------------------//

#endif
