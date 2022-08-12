//-----------------------------------------------------------------------------------//
//
// tinySPH:
//
//-----------------------------------------------------------------------------------//

void MeasureTotalKineticEnergyOfSystem(Variables<double, double> &V,
                                       ConstantVariables &C,
                                       double time,
                                       std::string fileName)
{

  double EkinTot = 0.;
  double m = C.m;

  for(unsigned int i = 0; i < V.N; i++)
  {
    EkinTot += 0.5*dot(V.v[i], V.v[i])*m;
  }

  std::ofstream fileWL;
  fileWL.open(fileName, std::ios_base::app);
  fileWL << time << " " << EkinTot << std::endl;
  fileWL.close();

}

//-----------------------------------------------------------------------------------//
