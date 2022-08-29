//-----------------------------------------------------------------------------------//
//
// tinySPH: still water case 2D (MDBC) - generate particles
//
//-----------------------------------------------------------------------------------//

 double endTime = 0.01; //seconds
 //double endTime = 0.0001; //seconds
 double initTimeStep = 0.0002; //seconds
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 10000; //steps
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC";
 std::string caseResults = "/home/tomas/Documents/temp/tinySPH_glob/MDBC";

//-----------------------------------------------------------------------------------//

 double smoothingLength = 0.04; //h
 double initialParticleDistance = 0.02; //dp
 double particleMass = 0.4; //m
 double numericalSpeedOfSound = 42.48; //cs
 double initialDensity = 1000.; //rho0
 double viscosityCoef = 0.01; //eta,alpha

//-----------------------------------------------------------------------------------//

 int LT = 24; int MT = 49; int RT = 74;
 int LM = 23; int MM = 48; int RM = 73;
 int LB = 22; int MB = 47; int RB = 72;

//-----------------------------------------------------------------------------------//
