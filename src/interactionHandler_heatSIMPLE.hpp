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

//-----------------------------------------------------------------------------------//

template<
typename SOLID_SOLID,
typename SOLID_BOUND,
typename BOUND_UPDATE,
typename KERNEL
>
class InteractionHandler{
public:

  void Interact(Variables<double, double> &HEAT_SIMPLEsolid,
                BoundVars<double, double> &HEAT_SIMPLEbound,
                ConstantVariables HEAT_SIMPLEconstants,
                bool updateBoundary);

  /*
  void UpdateBoundary(Variables<double, double> &WCSPHfluid,
                      BoundVars<double, double> &WCSPHbound,
                      ConstantVariables WCSPHconstants);
  */

  InteractionHandler(double h,
                     Variables<double, double> &HEAT_SIMPLEsolid,
                     Variables<double, double> &HEAT_SIMPLEbound);

  ~InteractionHandler(){};

private:
  CompactNSearch::NeighborhoodSearch nsearch;
  unsigned int positions_id;
  unsigned int positions_id_wall;

};

//-----------------------------------------------------------------------------------//

template<
typename SOLID_SOLID,
typename SOLID_BOUND,
typename BOUND_UPDATE,
typename KERNEL
>
InteractionHandler<
SOLID_SOLID,
SOLID_BOUND,
BOUND_UPDATE,
KERNEL>::InteractionHandler(double h,
                            Variables<double, double> &HEAT_SIMPLEsolid,
                            Variables<double, double> &HEAT_SIMPLEbound)
                            : nsearch(h, true)
{

  positions_id = nsearch.add_point_set(&HEAT_SIMPLEsolid.r.front().x,
                                       HEAT_SIMPLEsolid.N);
  positions_id_wall = nsearch.add_point_set(&HEAT_SIMPLEbound.r.front().x,
                                            HEAT_SIMPLEbound.N);

}

//-----------------------------------------------------------------------------------//

template<
typename SOLID_SOLID,
typename SOLID_BOUND,
typename BOUND_UPDATE,
typename KERNEL
>
void InteractionHandler<
SOLID_SOLID,
SOLID_BOUND,
BOUND_UPDATE,
KERNEL>::Interact(Variables<double, double> &HEAT_SIMPLEsolid,
                  BoundVars<double, double> &HEAT_SIMPLEbound,
                  ConstantVariables HEAT_SIMPLEconstants,
                  bool updateBoundary)
{


  nsearch.find_neighbors();

  // Update fluid particles
  CompactNSearch::PointSet const& ps_1 = nsearch.point_set(positions_id);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_1.n_points(); ++i)
  {
    InteractionData<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyDataIn(HEAT_SIMPLEsolid, i);

    //Get point set 1 neighbors of point set 1.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(HEAT_SIMPLEsolid, pid);
      SOLID_SOLID::template SolidSolidInteraction<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
    }

    //Get point set 1 neighbors of point set 2.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
      ptcJ.CopyBoundaryDataIn(HEAT_SIMPLEbound, pid);
      SOLID_BOUND::template SolidBoundInteraction<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
    }

    ptcI.CopyDataOut(HEAT_SIMPLEsolid, i);
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
      ptcI.CopyDataIn(HEAT_SIMPLEbound, i);
      for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
      {
        // Return the point id of the jth neighbor of the ith particle in the point_set_1.
        const unsigned int pid = ps_2.neighbor(positions_id, i, j);
        ptcJ.CopyDataIn(HEAT_SIMPLEsolid, pid);
        BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
      }

      ptcI.CopyBoundaryDataOut(HEAT_SIMPLEbound, i);

    }
  }

}

//-----------------------------------------------------------------------------------//

/*
template<
typename SOLID_SOLID,
typename SOLID_BOUND,
typename BOUND_UPDATE,
typename KERNEL
>
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
KERNEL>::UpdateBoundary(Variables<double, double> &HEAT_SIMPLEsolid,
                        BoundVars<double, double> &HEAT_SIMPLEbound,
                        ConstantVariables HEAT_SIMPLEconstants)
{

  nsearch.find_neighbors();

  // Update boundary particles
  CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_wall);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_2.n_points(); ++i)
  {

    InteractionData<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyDataIn(HEAT_SIMPLEbound, i);
    for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
    {
      // Return the point id of the jth neighbor of the ith particle in the point_set_1.
      const unsigned int pid = ps_2.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(HEAT_SIMPLEsolid, pid);
      BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL>(ptcI, ptcJ, HEAT_SIMPLEconstants);
    }

    ptcI.CopyBoundaryDataOut(HEAT_SIMPLEbound, i);

  }

}
*/

//-----------------------------------------------------------------------------------//

#endif
