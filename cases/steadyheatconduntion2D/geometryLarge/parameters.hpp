//-----------------------------------------------------------------------------------//
//
// tinySPH: steady heat conduction 2D - parameters
//
//-----------------------------------------------------------------------------------//

 double endTime = 1.0;                                  // [seconds]
 double initTimeStep = 0.001;                           // [seconds]
 unsigned int stepEnd = ceil( endTime/initTimeStep );

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 100;                         // [steps]
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/steadyheatconduntion2D";
 std::string caseResults = "/home/tomas/Documents/temp/devel/heat";

//-----------------------------------------------------------------------------------//

 double smoothingLength = 0.04001;                      //h - smoothingLength
 double initialParticleDistance = 0.02;                 //dp - initial particle spacing
 double particleMass = 0.4;                             //m - particle mass
 double initialDensity = 1000.;                         //rho0 - referential density
 double heatConductivity = 300;                           //eta,alpha - viscosity value
 double specificHeat = 300;                           //eta,alpha - viscosity value

//-----------------------------------------------------------------------------------//

