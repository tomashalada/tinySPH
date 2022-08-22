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
#include "generalTools.hpp"

//-----------------------------------------------------------------------------------//

#include "../cases/stillwater2D_withMDBC_periodic/periodicBC.hpp"

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
                PeriodicBC &WCSPHperiodic);

  void UpdateBoundary(Variables<double, double> &WCSPHfluid,
                      BoundVars<double, double> &WCSPHbound,
                      ConstantVariables WCSPHconstants);

  void FindNeighbors();

  InteractionHandler(double h,
                     Variables<double, double> &WCSPHfluid,
                     BoundVars<double, double> &WCSPHbound,
                     PeriodicBC &WCSPHperiodic);

  ~InteractionHandler(){};

  //temporary periodic BC
  /*
  unsigned int leftOverlap_FluidN = 0;
  std::vector<Vector3d> leftOverlap_r;
  std::vector<unsigned int> leftOverlap_idx;
  unsigned int rightOverlap_FluidN = 0;
  std::vector<Vector3d> rightOverlap_r;
  std::vector<unsigned int> rightOverlap_idx;
  */

private:
  CompactNSearch::NeighborhoodSearch nsearch;
  unsigned int positions_id;
  unsigned int positions_id_wall;
  unsigned int positions_id_ghostNodes;

  unsigned int positions_id_leftOverlap;
  unsigned int positions_id_rightOverlap;

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
                            BoundVars<double, double> &WCSPHbound,
                            PeriodicBC &WCSPHperiodic)
                            : nsearch(h, true)
{

  positions_id = nsearch.add_point_set(&WCSPHfluid.r.front().x, WCSPHfluid.N);
  positions_id_wall = nsearch.add_point_set(&WCSPHbound.r.front().x, WCSPHbound.N);
  positions_id_ghostNodes = nsearch.add_point_set(&WCSPHbound.gn.front().x, WCSPHbound.N);

  //periodic boundary condition
  positions_id_leftOverlap = nsearch.add_point_set(&WCSPHperiodic.leftOverlap_r.front().x, WCSPHperiodic.leftOverlap_r.size());
  positions_id_rightOverlap = nsearch.add_point_set(&WCSPHperiodic.rightOverlap_r.front().x, WCSPHperiodic.rightOverlap_r.size());

  nsearch.set_active(positions_id_wall, positions_id_ghostNodes, false);
  nsearch.set_active(positions_id_wall, positions_id_wall, false);
  nsearch.set_active(positions_id_wall, positions_id, false);

  nsearch.set_active(positions_id_ghostNodes, positions_id_ghostNodes, false);
  nsearch.set_active(positions_id_ghostNodes, positions_id_wall, false);

  nsearch.set_active(positions_id, positions_id_ghostNodes, false);

  //activation table for boundary conditions

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
                  PeriodicBC &WCSPHperiodic)
{

  //Compute pressure from density
  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  nsearch.find_neighbors();
  std::cout << "hit" << std::endl;

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

    //Resolve period bounary condition
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_leftOverlap, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_leftOverlap, i, j);
      if(pid < WCSPHperiodic.leftOverlap_FluidN)
        ptcJ.CopyDataIn(WCSPHfluid, WCSPHperiodic.leftOverlap_idx[pid]);
      else
        ptcJ.CopyDataIn(WCSPHbound, WCSPHperiodic.leftOverlap_idx[pid]);

      ptcJ.r = WCSPHperiodic.leftOverlap_r[pid];
      FLUID_FLUID::template FluidFluidInteraction<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    //Resolve period bounary condition
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_rightOverlap, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_rightOverlap, i, j);
      if(pid < WCSPHperiodic.rightOverlap_FluidN)
        ptcJ.CopyDataIn(WCSPHfluid, WCSPHperiodic.rightOverlap_idx[pid]);
      else
        ptcJ.CopyDataIn(WCSPHbound, WCSPHperiodic.rightOverlap_idx[pid]);

      ptcJ.r = WCSPHperiodic.rightOverlap_r[pid];
      FLUID_FLUID::template FluidFluidInteraction<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    //Get point set 1 neighbors of point set 2.
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
      ptcJ.CopyBoundaryDataIn(WCSPHbound, pid);
      FLUID_BOUND::template FluidBoundInteraction<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }


    ptcI.a.z -= 9.81; //apply external forces
    ptcI.CopyDataOut(WCSPHfluid, i);
  }

  // Update boundary particles
  CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_ghostNodes);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_2.n_points(); ++i)
  {

    InteractionDataGhostNode<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyBoundaryDataIn(WCSPHbound, i);

    for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
    {
      // Return the point id of the jth neighbor of the ith particle in the point_set_1.
      const unsigned int pid = ps_2.neighbor(positions_id, i, j);
      ptcJ.CopyDataIn(WCSPHfluid, pid);
      BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    //Resolve period bounary condition
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_leftOverlap, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_leftOverlap, i, j);
      if(pid < WCSPHperiodic.leftOverlap_FluidN)
        ptcJ.CopyDataIn(WCSPHfluid, WCSPHperiodic.leftOverlap_idx[pid]);
      else
        continue;

      ptcJ.r = WCSPHperiodic.leftOverlap_r[pid];
      BOUND_UPDATE::template UpdateBoundaryData<double, double, KERNEL, DIFFUSIVE_TERM, VISOUCS_TERM>(ptcI, ptcJ, WCSPHconstants);
    }

    //Resolve period bounary condition
    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_rightOverlap, i); ++j)
    {
      const unsigned int pid = ps_1.neighbor(positions_id_rightOverlap, i, j);
      if(pid < WCSPHperiodic.rightOverlap_FluidN)
        ptcJ.CopyDataIn(WCSPHfluid, WCSPHperiodic.rightOverlap_idx[pid]);
      else
        continue;

      ptcJ.r = WCSPHperiodic.rightOverlap_r[pid];
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
  CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_ghostNodes);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < ps_2.n_points(); ++i)
  {

    InteractionDataGhostNode<double, double> ptcI;
    InteractionData<double, double> ptcJ;
    ptcI.CopyBoundaryDataIn(WCSPHbound, i);
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
void InteractionHandler<
FLUID_FLUID,
FLUID_BOUND,
BOUND_UPDATE,
DIFFUSIVE_TERM,
VISOUCS_TERM,
KERNEL>::FindNeighbors()
{

  nsearch.find_neighbors();

}

//-----------------------------------------------------------------------------------//

#endif
