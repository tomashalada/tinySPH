//-----------------------------------------------------------------------------------//
//
// tinySPH: CASE_DAMBREAK2D - parameters
//
//-----------------------------------------------------------------------------------//

 double endTime = 0.51; //seconds
 double initTimeStep = 0.000035; //seconds
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 250; //steps
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withMDBC";
 std::string caseResults = "/home/tomas/Documents/temp/tinySPH_dev/dambreak2D_withMDBC";

//-----------------------------------------------------------------------------------//


