//-----------------------------------------------------------------------------------//
//
// tinySPH: CASE_DAMBREAK2D - parameters
//
//-----------------------------------------------------------------------------------//

double endTime = 0.51; //seconds
//double endTime = 0.0001; //seconds
double initTimeStep = 0.00005; //seconds
unsigned int stepEnd = ceil(endTime/initTimeStep);

//-----------------------------------------------------------------------------------//

unsigned int saveOutput = 1000; //steps
 std::string casePath = "/home/tomas/Documents/temp/tinySPH/dambreak2D_withMDBC";

//-----------------------------------------------------------------------------------//


