//-----------------------------------------------------------------------------------//
//
// tinySPH: still water case 2D (MDBC) - generate particles
//
//-----------------------------------------------------------------------------------//

 double endTime = 60.01; //seconds
 //double endTime = 0.0001; //seconds
 double initTimeStep = 0.00005; //seconds
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 60000; //steps
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withMDBC";
 std::string caseResults = "/home/tomas/Documents/temp/tinySPH_glob/MDBC";

//-----------------------------------------------------------------------------------//

 double smoothingLength = 0.02; //h
 double initialParticleDistance = 0.01; //dp
 double particleMass = 0.1; //m
 double numericalSpeedOfSound = 42.48; //cs
 double initialDensity = 1000.; //rho0
 double viscosityCoef = 0.01; //eta,alpha

//-----------------------------------------------------------------------------------//

 int LT = 49; int MT = 99; int RT = 149;
 int LM = 48; int MM = 98; int RM = 148;
 int LB = 47; int MB = 97; int RB = 147;

//-----------------------------------------------------------------------------------//
