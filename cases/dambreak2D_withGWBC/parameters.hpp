//-----------------------------------------------------------------------------------//
//
// tinySPH: CASE_DAMBREAK2D - parameters
//
//-----------------------------------------------------------------------------------//

 double endTime = 0.31; //seconds
 double initTimeStep = 0.00003; //seconds
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 1000; //steps
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/dambreak2D_withGWBC";
 std::string caseResults = "/home/tomas/Documents/temp/tinySPH/dambreak2D_withGWBC";

//-----------------------------------------------------------------------------------//


