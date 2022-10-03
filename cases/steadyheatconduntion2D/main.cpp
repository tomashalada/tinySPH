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

//-----------------------------------------------------------------------------------//

#include "variables.hpp"                       // variables of used model
#include "interactionHandler_heatSIMPLE.hpp"   // interaction manager
#include "integration.hpp"                     // integrators (connected to the model)

//-----------------------------------------------------------------------------------//

/*
#include "interpolationStructures.hpp"  // definition of particle-particle interactions
#include "postprocessor.hpp"            // interpolation and measuretools
*/

//-----------------------------------------------------------------------------------//

#include "vector.hpp"                   // data structures
#include "generalTools.hpp"             // tools to work with string, names and paths
#include "writeParticleData.hpp"        // simple .ptcs writer

//-----------------------------------------------------------------------------------//

#include "parameters.hpp"               // case configuration parameters

//-----------------------------------------------------------------------------------//

int main(){

  //Set the case
  ConstantVariables HEAT_EQNconstants(
  smoothingLength, //h
  particleMass, //m
  initialDensity, //rho0
  heatConductivity, //lambda
  specificHeat, //cp
  initialParticleDistance //dp
  );

  //Read fluid data
  SolidVars<double, double> HEAT_EQNsolid;
  HEAT_EQNsolid.initializeWithGeometryFile(caseFolder + "/heatconduction_solid.ptcs");
  std::cout << "Number of fluid particles: " << HEAT_EQNsolid.N << std::endl;

  //Read boundary data
  BoundVars<double, double> HEAT_EQNbound;
  HEAT_EQNbound.initializeWithGeometryFile(caseFolder + "/heatconduction_boundary.ptcs");
  std::cout << "Number of boundary particles: " << HEAT_EQNbound.N << std::endl;

  InteractionHandler<
  heatSIMPLE_SOLIDSOLID,                // solid-solid interaction model
  heatSIMPLE_SOLIDBOUND,                // solid-boundary interaction model
  heatSIMPLE_BC,                        // wall particle update model
  WendlandKernel                        // SPH kernel function
  > HEAT_EQN(HEAT_EQNconstants.h*2, HEAT_EQNsolid, HEAT_EQNbound);

//-----------------------------------------------------------------------------------//

  /*
  //Interpolation nodes
  FluidVars<double, double> WCSPHinterpolationPlane;
  WCSPHinterpolationPlane.initializeWithGeometryFile(caseFolder + "/stillwater_interpolationPlane.ptcs");
  std::cout << "Number of nodes in interpolation plane: " << WCSPHinterpolationPlane.N << std::endl;

  PostProcessingHandler<
  WCSPH_INTERPOLATION,
  WendlandKernel
  > WCSPHmeasurement(WCSPHconstants.h*2, WCSPHinterpolationPlane, WCSPHfluid, WCSPHbound);
  */

//-----------------------------------------------------------------------------------//
// Verlet integrator

EulerScheme HEAT_EQNEuler(&HEAT_EQNsolid, initTimeStep);

for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  HEAT_EQN.Interact(HEAT_EQNsolid, HEAT_EQNbound, HEAT_EQNconstants, false);

  HEAT_EQNEuler.ComputeStep();

  if(step % saveOutput == 0){
    writeParticleData(HEAT_EQNsolid, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/SOLID/solid", step));
    writeParticleData(HEAT_EQNbound, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/BOUND/bound", step));

  //writeParticleData(WCSPHinterpolationPlane, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/INTERPOLATION/interpolation", step));
  }

}

//-----------------------------------------------------------------------------------//

  std::cout << "Done..." << std::endl;

  return EXIT_SUCCESS;

}

