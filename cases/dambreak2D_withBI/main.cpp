#include <CompactNSearch>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>

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
#include "interactionHandler.hpp"       // interaction manager
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
  0.0070711, //h
  0.025, //m
  34.3, //cs
  0.02, //eta
  1000., //rho0
  0.005 //dp
  );

  //Read fluid data
  FluidVars<double, double> WCSPHfluid;
  WCSPHfluid.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withBI/dambreak_fluid.ptcs");
  std::cout << "Number of fluid particles: " << WCSPHfluid.N << std::endl;

  //Read boundary data
  BoundVars<double, double> WCSPHbound;
  WCSPHbound.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withBI/dambreak_wall.ptcs");
  std::cout << "Number of boundary particles: " << WCSPHbound.N << std::endl;

  InteractionHandler<
  WCSPHBI_FLUIDFLUID, //fluid-fluid interaction scheme
  WCSPHBI_FLUIDBOUND_BI, //fluid-boundary interaction scheme
  BOUNDARY_SIMPLE, //way to update boundary particles
  DT_Molteni,
  VISCOSITY_Artificial,
  WendlandKernel//diffusive term
  > WCSPH(WCSPHconstants.h*2, WCSPHfluid, WCSPHbound);

//-----------------------------------------------------------------------------------//

  //Interpolation nodes
  FluidVars<double, double> WCSPHinterpolationPlane;
  WCSPHinterpolationPlane.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withBI/dambreak_interpolationPlane.ptcs");
  std::cout << "Number of nodes in interpolation plane: " << WCSPHinterpolationPlane.N << std::endl;

  PostProcessingHandler<
  WCSPH_INTERPOLATION,
  WendlandKernel
  > WCSPHmeasurement(WCSPHconstants.h*2, WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound);

//-----------------------------------------------------------------------------------//
// Symplectic integrator

//:SymplecticScheme WCSPHSymplectic(&WCSPHfluid, 0.00003);
//:
//:for(int step = 0; step < stepEnd + 1; step++)
//:{
//:
//:	std::cout << "STEP: " << step << std::endl;
//:	DensityToPressure(WCSPHfluid, WCSPHconstants);
//:	DensityToPressure(WCSPHbound, WCSPHconstants);
//:
//:	WCSPHSymplectic.ComputePredictor();
//:	WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);
//:
//:	DensityToPressure(WCSPHfluid, WCSPHconstants);
//:	DensityToPressure(WCSPHbound, WCSPHconstants);
//:
//:	WCSPHSymplectic.ComputeCorrector();
//:	WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);
//:
//:	if(step % saveOutput == 0){
//:	writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(casePath + "/OUTPUT/FLUID/fluid", step));
//:	writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(casePath + "/OUTPUT/BOUND/bound", step));
//:	}
//:
//:}

//-----------------------------------------------------------------------------------//
// Verlet integrator

VerletScheme WCSPHVerlet(&WCSPHfluid, 0.00003);

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

