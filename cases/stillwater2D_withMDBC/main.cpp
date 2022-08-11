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
#include "interactionHandlerMDBC.hpp"   // interaction manager
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

int main(){

  //Set the case
  ConstantVariables WCSPHconstants(
  0.04, //h
  0.4, //m
  42.48, //cs
  0.02, //eta
  1000., //rho0
  0.02 //dp
  );

  //Read fluid data
  FluidVars<double, double> WCSPHfluid;
  WCSPHfluid.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC/stillwater_fluid.ptcs");
  			std::cout << "Number of fluid particles: " << WCSPHfluid.N << std::endl;

  //Read boundary data
  BoundVars<double, double> WCSPHbound;
  WCSPHbound.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC/stillwater_wall.ptcs");
  			std::cout << "Number of boundary particles: " << WCSPHbound.N << std::endl;

  InteractionHandler<
  WCSPH_FLUIDFLUID,                     // fluid-fluid interaction model
  WCSPH_FLUIDBOUND_DBC,                 // fluid-wall interaction model
  WCSPH_MDBC,                           // wall particle update model
  DT_Molteni,                           // diffusive term
  VISCOSITY_Artificial,                 // viscosity term
  WendlandKernel                        // SPH kernel function
  > WCSPH(WCSPHconstants.h*2, WCSPHfluid, WCSPHbound);

//-----------------------------------------------------------------------------------//

  //Interpolation nodes
  FluidVars<double, double> WCSPHinterpolationPlane;
  WCSPHinterpolationPlane.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC/stillwater_interpolationPlane.ptcs");
  std::cout << "Number of nodes in interpolation plane: " << WCSPHinterpolationPlane.N << std::endl;

  PostProcessingHandler<
  WCSPH_INTERPOLATION,
  WendlandKernel
  > WCSPHmeasurement(WCSPHconstants.h*2, WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound);


//-----------------------------------------------------------------------------------//

  PressureToDensity(WCSPHfluid, WCSPHconstants);
  WCSPH.UpdateBoundary(WCSPHfluid, WCSPHbound, WCSPHconstants);

//-----------------------------------------------------------------------------------//
// Symplectic integrator

/*
SymplecticScheme WCSPHSymplecticFluid(&WCSPHfluid, 0.000035);
SymplecticScheme WCSPHSymplecticBound(&WCSPHbound, 0.000035);

for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  WCSPHSymplecticFluid.ComputePredictor();
  WCSPHSymplecticBound.ComputePredictor();
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);

  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  WCSPHSymplecticFluid.ComputeCorrector();
  WCSPHSymplecticBound.ComputeCorrector();
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);

  if(step % saveOutput == 0){
    writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(casePath + "/OUTPUT/FLUID/fluid", step));
    writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(casePath + "/OUTPUT/BOUND/bound", step));
  }

}
*/

//-----------------------------------------------------------------------------------//
// Verlet integrator

VerletScheme WCSPHVerlet(&WCSPHfluid, 0.0001);


for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);

  if(step % 20 == 0){
  WCSPHVerlet.ComputeStepEulerForStability();
  }
  else{
  WCSPHVerlet.ComputeStep();
  }

  if(step % saveOutput == 0){
    writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(casePath + "/OUTPUT/FLUID/fluid", step));
    writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(casePath + "/OUTPUT/BOUND/bound", step));

  WCSPHmeasurement.Interpolate(WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound, WCSPHconstants);
  writeParticleData(WCSPHinterpolationPlane, stepToNameWithPtcsExtension(casePath + "/OUTPUT/INTERPOLATION/interpolation", step));
  }

}


//-----------------------------------------------------------------------------------//

  std::cout << "Done..." << std::endl;

  return EXIT_SUCCESS;

}
