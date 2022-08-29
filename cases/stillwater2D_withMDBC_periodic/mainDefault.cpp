#include <CompactNSearch>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

#include <omp.h>
#include <algorithm>


#include <string>
#include <string_view>
#include <cstdlib>

//-----------------------------------------------------------------------------------//

#include "interactionStructures.hpp"    // definition of particle-particle interactions
#include "kernel.hpp"                   // kernel functions
#include "diffusiveTerms.hpp"           // structures with diffusive terms
#include "visoucsTerms.hpp"             // structures with visoucs forces

//-----------------------------------------------------------------------------------//

#include "variables.hpp"                // variables of used model
//#include "interactionHandlerMDBC.hpp" // interaction manager
#include "interactionHandlerMDBC_withPeriodicBC.hpp"   // interaction manager with aux periodic BC
#include "integration.hpp"              // integrators (connected to the model)

//-----------------------------------------------------------------------------------//

#include "interpolationStructures.hpp"  // definition of particle-particle interactions
#include "postprocessor.hpp"            // interpolation and measuretools

//-----------------------------------------------------------------------------------//

#include "vector.hpp"                   // data structures
#include "generalTools.hpp"             // tools to work with string, names and paths
#include "writeParticleData.hpp"        // simple .ptcs writer

//-----------------------------------------------------------------------------------//

#include "parameters.hpp"               // case configuration parameters

//-----------------------------------------------------------------------------------//

#include "measureKineticEnergy.hpp"     // cutstom measurement tools

//-----------------------------------------------------------------------------------//

#include "periodicBC.hpp"

//-----------------------------------------------------------------------------------//

int main(){

  //Set the case
  ConstantVariables WCSPHconstants(
  smoothingLength, //h
  particleMass, //m
  numericalSpeedOfSound, //cs
  viscosityCoef, //eta
  initialDensity, //rho0
  initialParticleDistance //dp
  );

  //Read fluid data
  FluidVars<double, double> WCSPHfluid;
  WCSPHfluid.initializeWithGeometryFile(caseFolder + "/stillwater_fluid.ptcs");
  std::cout << "Number of fluid particles: " << WCSPHfluid.N << std::endl;

  //Read boundary data
  BoundVars<double, double> WCSPHbound;
  WCSPHbound.initializeWithGeometryFile(caseFolder + "/stillwater_wall.ptcs");
  std::cout << "Number of boundary particles: " << WCSPHbound.N << std::endl;

//-----------------------------------------------------------------------------------//

  // InitPeriodicBC:
  PeriodicBC WCSPHperiodic;
  WCSPHperiodic.ApplyPeriodicBoundaryCondition(WCSPHfluid,
                                               WCSPHbound,
                                               0. + (WCSPHconstants.h*2.01),
                                               1.0 - (WCSPHconstants.h*2.01),
                                               0.+ WCSPHconstants.dp,
                                               1.0 - WCSPHconstants.dp);

//-----------------------------------------------------------------------------------//

  InteractionHandler<
  WCSPH_FLUIDFLUID,                     // fluid-fluid interaction model
  WCSPH_FLUIDBOUND_DBC,                 // fluid-wall interaction model
  //WCSPH_MDBC,                         // wall particle update model
  WCSPH_MDBCvelocity,                   // wall particle update model
  //DT_Molteni,                         // diffusive term
  DT_Fourtakas,                         // diffusive term
  VISCOSITY_Artificial,                 // viscosity term
  WendlandKernel                        // SPH kernel function
  > WCSPH(WCSPHconstants.h*2, WCSPHfluid, WCSPHbound, WCSPHperiodic);

//-----------------------------------------------------------------------------------//

  //Interpolation nodes
  FluidVars<double, double> WCSPHinterpolationPlane;
  WCSPHinterpolationPlane.initializeWithGeometryFile(caseFolder + "/stillwater_interpolationPlane.ptcs");
  std::cout << "Number of nodes in interpolation plane: " << WCSPHinterpolationPlane.N << std::endl;

  PostProcessingHandler<
  WCSPH_INTERPOLATION,
  WendlandKernel
  > WCSPHmeasurement(WCSPHconstants.h*2, WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound);

  //Measure kineticEnergy and particle trajectories
  MEASUREMENT_TotalKineticEnergyOfSystem WCSPHEkinTot(stepEnd, initTimeStep);
  MEASUREMENT_TrackParticleMovement WCSPHtrackParticles(stepEnd, initTimeStep);

//-----------------------------------------------------------------------------------//

  PressureToDensity(WCSPHfluid, WCSPHconstants);
  WCSPH.UpdateBoundary(WCSPHfluid, WCSPHbound, WCSPHconstants);

//-----------------------------------------------------------------------------------//
// Symplectic integrator

SymplecticScheme WCSPHSymplecticFluid(&WCSPHfluid, initTimeStep);

for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  WCSPHperiodic.ApplyPeriodicBoundaryCondition(WCSPHfluid,
                                               WCSPHbound,
                                               0. + (WCSPHconstants.h*2.01),
                                               1.0 - (WCSPHconstants.h*2.01),
                                               0.+ WCSPHconstants.dp,
                                               1.0 - WCSPHconstants.dp);



  WCSPHSymplecticFluid.ComputePredictor();

  //WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants, false);
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants, WCSPHperiodic);

  WCSPHSymplecticFluid.ComputeCorrector();
  //WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants, true);
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants, WCSPHperiodic);

  if(step % saveOutput == 0)
  {
    writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(caseResults + "/FLUID/fluid", step));
    writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(caseResults + "/BOUND/bound", step));
    WCSPHmeasurement.Interpolate(WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound, WCSPHconstants);
    writeParticleData(WCSPHinterpolationPlane, stepToNameWithPtcsExtension(caseResults + "/INTERPOLATION/interpolation", step));
  }

  WCSPHEkinTot.ComputeKineticEnergy(WCSPHfluid, WCSPHconstants);
  WCSPHtrackParticles.TrackParticles(WCSPHfluid, LT, LM, LB, MT, MM, MB, RT, RM, RB);

}


//-----------------------------------------------------------------------------------//
// Verlet integrator

/*
VerletScheme WCSPHVerlet(&WCSPHfluid, initTimeStep);


for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  //WCSPH.ApplyPeriodicBoundaryCondition(WCSPHfluid, WCSPHbound, 0. + (WCSPHconstants.h*1.05), 0.8 - (WCSPHconstants.h*1.05), 0., 0.8, step);
  WCSPHperiodic.ApplyPeriodicBoundaryCondition(WCSPHfluid,
                                               WCSPHbound,
                                               0. + (WCSPHconstants.h*2.01),
                                               0.8 - (WCSPHconstants.h*2.01),
                                               0.+ WCSPHconstants.dp,
                                               0.8 - WCSPHconstants.dp);

  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants, WCSPHperiodic);

  if(step % 20 == 0){
  WCSPHVerlet.ComputeStepEulerForStability();
  }
  else{
  WCSPHVerlet.ComputeStep();
  }

  if(step % saveOutput == 0){
    writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/FLUID/fluid", step));
    writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/BOUND/bound", step));

  WCSPHmeasurement.Interpolate(WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound, WCSPHconstants);
  writeParticleData(WCSPHinterpolationPlane, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/INTERPOLATION/interpolation", step));
  }

  //Custom measuretools
  WCSPHEkinTot.ComputeKineticEnergy(WCSPHfluid, WCSPHconstants);
  WCSPHtrackParticles.TrackParticles(WCSPHfluid);

}
*/


//-----------------------------------------------------------------------------------//

  WCSPHEkinTot.WriteTotalKinetcEnergyToFile(caseResults + "/TotalKineticEnergy.dat");
  WCSPHtrackParticles.WriteParticleTrajectoryToFile(caseResults + "/ParticleTrajectory"); //no .dat here!
  WCSPHtrackParticles.WriteParticleVelocityToFile(caseResults + "/ParticleVelocity"); //no .dat here!


  std::cout << "Done..." << std::endl;

  return EXIT_SUCCESS;

}

