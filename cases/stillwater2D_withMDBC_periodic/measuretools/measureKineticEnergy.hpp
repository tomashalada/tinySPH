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

class MEASUREMENT_TotalKineticEnergyOfSystem
{

public:

  double dt;
  std::vector<double> EkinTot;

  //constructor
  MEASUREMENT_TotalKineticEnergyOfSystem(long int TotalSteps,
                                    double _dt)
                                    : dt(_dt)
  {
    //EkinTot.resize(TotalSteps);
  };

  //destructor
  ~MEASUREMENT_TotalKineticEnergyOfSystem(){};

  //Compute kinetic energy
  void ComputeKineticEnergy(Variables<double, double> &V,
                       ConstantVariables &C)
  {
    double EkinTot_local = 0.;
    double m = C.m;

    for(unsigned int i = 0; i < V.N; i++)
    {
      EkinTot_local += 0.5*dot(V.v[i], V.v[i])*m;
    }
    EkinTot.push_back(EkinTot_local);
  }

  //Write EkinTot to file
  void WriteTotalKinetcEnergyToFile(std::string fileName)
  {
    std::ofstream fileWL;
    fileWL.open(fileName);

    for(int i = 0; i < EkinTot.size(); i++)
      fileWL << dt*i << " " << EkinTot[i] << std::endl;

    fileWL.close();
  }

};


//-----------------------------------------------------------------------------------//
