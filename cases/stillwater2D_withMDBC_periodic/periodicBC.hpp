#ifndef PERIODIC_HPP
#define PERIODIC_HPP

#pragma once

//-----------------------------------------------------------------------------------//

class PeriodicBC
{
public:

  //temporary periodic BC
  unsigned int leftOverlap_FluidN = 0;
  std::vector<Vector3d> leftOverlap_r;
  std::vector<unsigned int> leftOverlap_idx;
  unsigned int rightOverlap_FluidN = 0;
  std::vector<Vector3d> rightOverlap_r;
  std::vector<unsigned int> rightOverlap_idx;

  void
  ApplyPeriodicBoundaryCondition(Variables<double, double> &WCSPHfluid,
                                 BoundVars<double, double> &WCSPHbound,
                                 double leftTreshold,
                                 double rightTreshold,
                                 double leftBoundary,
                                 double rightBoundary);
};

//-----------------------------------------------------------------------------------//

void PeriodicBC::ApplyPeriodicBoundaryCondition(Variables<double, double> &WCSPHfluid,
                                                BoundVars<double, double> &WCSPHbound,
                                                double leftTreshold,
                                                double rightTreshold,
                                                double leftBoundary,
                                                double rightBoundary)
{

  //Reset
  leftOverlap_FluidN = 0;
  leftOverlap_r.clear();
  leftOverlap_idx.clear();

  rightOverlap_FluidN = 0;
  rightOverlap_r.clear();
  rightOverlap_idx.clear();

  //Copy boundary overlaps
  for(int i = 0; i < WCSPHfluid.N; i++)
  {
    if (WCSPHfluid.r[i].x <= leftTreshold){
      leftOverlap_r.push_back({rightBoundary + WCSPHfluid.r[i].x, 0., WCSPHfluid.r[i].z});
      leftOverlap_idx.push_back(i);
      leftOverlap_FluidN++;
    }
    else if(WCSPHfluid.r[i].x >= rightTreshold){
      rightOverlap_r.push_back({WCSPHfluid.r[i].x - rightBoundary, 0., WCSPHfluid.r[i].z});
      rightOverlap_idx.push_back(i);
      rightOverlap_FluidN++;
    }
  }

  for(int i = 0; i < WCSPHbound.N; i++)
  {
    if (WCSPHbound.r[i].x <= leftTreshold){
      leftOverlap_r.push_back({rightBoundary + WCSPHbound.r[i].x, 0., WCSPHbound.r[i].z});
      leftOverlap_idx.push_back(i);
    }
    else if(WCSPHbound.r[i].x >= rightTreshold){
      rightOverlap_r.push_back({WCSPHbound.r[i].x - rightBoundary, 0., WCSPHbound.r[i].z});
      rightOverlap_idx.push_back(i);
    }
  }

  /*
  std::cout << "Information about periodic boundary update: " << std::endl;
  std::cout << " - leftOverlap_FluidN ....... " << leftOverlap_FluidN << std::endl;
  std::cout << " - leftOverlap_r.size() ..... " << leftOverlap_r.size() << std::endl;
  std::cout << " - leftOverlap_idx.size() ... " << leftOverlap_idx.size() << std::endl;
  std::cout << " - rightOverlap_FluidN ...... " << rightOverlap_FluidN << std::endl;
  std::cout << " - rightOverlap_r.size() .... " << rightOverlap_r.size() << std::endl;
  std::cout << " - rightOverlap_idx.size() .. " << rightOverlap_idx.size() << std::endl;
  */

  /*
  if(step == 0)
  {
    positions_id_leftOverlap = nsearch.add_point_set(&leftOverlap_r.front().x, leftOverlap_r.size());
    positions_id_rightOverlap = nsearch.add_point_set(&rightOverlap_r.front().x, rightOverlap_r.size());
  }

  std::cout << " - positions_id_leftOverlap .. " << positions_id_leftOverlap << std::endl;
  std::cout << " - positions_id_rightOverlap . " << positions_id_rightOverlap << std::endl;

  std::cout << std::endl;
  std::cout << "Other pointsets for comparasement: " << std::endl;
  std::cout << " - positions_id .............. " << positions_id << std::endl;
  std::cout << " - positions_id_wall ......... " << positions_id_wall << std::endl;
  std::cout << " - positions_id_ghostNodes ... " << positions_id_ghostNodes << std::endl;
  */

  /*
  std::cout << std::endl;
  std::cout << "Print my overlap sets: " << std::endl;
  std::cout << " - leftone ...........  {";
  for(auto &p: leftOverlap_r) std::cout << " [" << p.x << "," << p.z << "] ";
  std::cout << "}" << std::endl;
  std::cout << " - rightone ..........  {";
  for(auto &p: rightOverlap_r) std::cout << " [" << p.x << "," << p.z << "] ";
  std::cout << "}" << std::endl;
  std::cout << std::endl;
  */

}

//-----------------------------------------------------------------------------------//

#endif
