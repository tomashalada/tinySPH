//-----------------------------------------------------------------------------------//
//
// tinySPH: still water case 2D (MDBC) - generate particles
//
//-----------------------------------------------------------------------------------//

 double endTime = 4.01; //seconds
 //double endTime = 0.0001; //seconds
 double initTimeStep = 0.0002; //seconds
 unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

 unsigned int saveOutput = 10000; //steps
 std::string caseFolder = "/home/tomas/Documents/testovaci/cpp/TINYSPH_lib/cases/stillwater2D_withDBC";
 std::string caseResults = "/home/tomas/Documents/temp/tinySPH_dev/stillwater2D_withDBC";

//-----------------------------------------------------------------------------------//


