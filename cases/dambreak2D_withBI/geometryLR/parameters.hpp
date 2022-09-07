//-----------------------------------------------------------------------------------//
//
// tinySPH: CASE_DAMBREAK2D - parameters
//
//-----------------------------------------------------------------------------------//

 double endTime = 0.51;                                 // [seconds]
 double initTimeStep = 0.00003;                         // [seconds]
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 1000;                        // [steps]
 std::string caseFolder = "pathTemplate";
 std::string caseResults = "pathTemplate";

//-----------------------------------------------------------------------------------//

 double smoothingLength = 0.0070711*1.3;                //h - smoothingLength
 double initialParticleDistance = 0.005;                //dp - initial particle spacing
 double particleMass = 0.025;                           //m - particle mass
 double numericalSpeedOfSound = 34.3;                   //cs - numerical speed of sound
 double initialDensity = 1000.;                         //rho0 - referential density
 double viscosityCoef = 0.02;                           //eta,alpha - viscosity value

//-----------------------------------------------------------------------------------//

