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

class MEASUREMENT_TrackParticleMovement
{
  public:

  /* For dp = 0.01:
  LT = 24, LM = 23, LB = 22, MT = 49, MM = 48, MB = 47, RT = 74, RM = 73, RB = 72 */

  double dt;
  std::vector<Vector3d> r_LT, r_LM, r_LB, r_MT, r_MM, r_MB, r_RT, r_RM, r_RB;
  std::vector<Vector3d> v_LT, v_LM, v_LB, v_MT, v_MM, v_MB, v_RT, v_RM, v_RB;

  //constructor
  MEASUREMENT_TrackParticleMovement(long int TotalSteps,
                                    double _dt)
                                    : dt(_dt)
  {
    //EkinTot.resize(TotalSteps);
  };

  //destructor
  ~MEASUREMENT_TrackParticleMovement(){};

  void TrackParticles(Variables<double, double> &V)
  {
    r_LT.push_back(V.r[24]); v_LT.push_back(V.v[24]);
    r_LM.push_back(V.r[23]); v_LM.push_back(V.v[23]);
    r_LB.push_back(V.r[22]); v_LB.push_back(V.v[22]);

    r_MT.push_back(V.r[49]); v_MT.push_back(V.v[49]);
    r_MM.push_back(V.r[48]); v_MM.push_back(V.v[48]);
    r_MB.push_back(V.r[47]); v_MB.push_back(V.v[47]);

    r_RT.push_back(V.r[74]); v_RT.push_back(V.v[74]);
    r_RM.push_back(V.r[73]); v_RM.push_back(V.v[73]);
    r_RB.push_back(V.r[72]); v_RB.push_back(V.v[72]);
  }

  void WriteParticleTrajectoryToFile(std::string fileName)
  {
    // Writte LT
    std::ofstream fileLT;
    std::string fileNameLT = fileName + "LT.dat";
    fileLT.open(fileNameLT);
    for(int i = 0; i < r_LT.size(); i++)
      fileLT << dt*i << " " << r_LT[i].x << " " << r_LT[i].z << std::endl;
    fileLT.close();

    // Writte LM
    std::ofstream fileLM;
    std::string fileNameLM = fileName + "LM.dat";
    fileLM.open(fileNameLM);
    for(int i = 0; i < r_LM.size(); i++)
      fileLM << dt*i << " " << r_LM[i].x << " " << r_LM[i].z << std::endl;
    fileLM.close();

    // Writte LB
    std::ofstream fileLB;
    std::string fileNameLB = fileName + "LB.dat";
    fileLB.open(fileNameLB);
    for(int i = 0; i < r_LB.size(); i++)
      fileLB << dt*i << " " << r_LB[i].x << " " << r_LB[i].z << std::endl;
    fileLB.close();

    // Writte MT
    std::ofstream fileMT;
    std::string fileNameMT = fileName + "MT.dat";
    fileMT.open(fileNameMT);
    for(int i = 0; i < r_MT.size(); i++)
      fileMT << dt*i << " " << r_MT[i].x << " " << r_MT[i].z << std::endl;
    fileMT.close();

    // Writte MM
    std::ofstream fileMM;
    std::string fileNameMM = fileName + "MM.dat";
    fileMM.open(fileNameMM);
    for(int i = 0; i < r_MM.size(); i++)
      fileMM << dt*i << " " << r_MM[i].x << " " << r_MM[i].z << std::endl;
    fileMM.close();

    // Writte MB
    std::ofstream fileMB;
    std::string fileNameMB = fileName + "MB.dat";
    fileMB.open(fileNameMB);
    for(int i = 0; i < r_MB.size(); i++)
      fileMB << dt*i << " " << r_MB[i].x << " " << r_MB[i].z << std::endl;
    fileMB.close();

    // Writte RT
    std::ofstream fileRT;
    std::string fileNameRT = fileName + "RT.dat";
    fileRT.open(fileNameRT);
    for(int i = 0; i < r_RT.size(); i++)
      fileRT << dt*i << " " << r_RT[i].x << " " << r_RT[i].z << std::endl;
    fileRT.close();

    // Writte RM
    std::ofstream fileRM;
    std::string fileNameRM = fileName + "RM.dat";
    fileRM.open(fileNameRM);
    for(int i = 0; i < r_RM.size(); i++)
      fileRM << dt*i << " " << r_RM[i].x << " " << r_RM[i].z << std::endl;
    fileRM.close();

    // Writte RB
    std::ofstream fileRB;
    std::string fileNameRB = fileName + "RB.dat";
    fileRB.open(fileNameRB);
    for(int i = 0; i < r_RB.size(); i++)
      fileRB << dt*i << " " << r_RB[i].x << " " << r_RB[i].z << std::endl;
    fileRB.close();
  }

  void WriteParticleVelocityToFile(std::string fileName)
  {
    // Writte LT
    std::ofstream fileLT;
    std::string fileNameLT = fileName + "LT.dat";
    fileLT.open(fileNameLT);
    for(int i = 0; i < v_LT.size(); i++)
      fileLT << dt*i << " " << v_LT[i].x << " " << v_LT[i].z << std::endl;
    fileLT.close();

    // Writte LM
    std::ofstream fileLM;
    std::string fileNameLM = fileName + "LM.dat";
    fileLM.open(fileNameLM);
    for(int i = 0; i < v_LM.size(); i++)
      fileLM << dt*i << " " << v_LM[i].x << " " << v_LM[i].z << std::endl;
    fileLM.close();

    // Writte LB
    std::ofstream fileLB;
    std::string fileNameLB = fileName + "LB.dat";
    fileLB.open(fileNameLB);
    for(int i = 0; i < v_LB.size(); i++)
      fileLB << dt*i << " " << v_LB[i].x << " " << v_LB[i].z << std::endl;
    fileLB.close();

    // Writte MT
    std::ofstream fileMT;
    std::string fileNameMT = fileName + "MT.dat";
    fileMT.open(fileNameMT);
    for(int i = 0; i < v_MT.size(); i++)
      fileMT << dt*i << " " << v_MT[i].x << " " << v_MT[i].z << std::endl;
    fileMT.close();

    // Writte MM
    std::ofstream fileMM;
    std::string fileNameMM = fileName + "MM.dat";
    fileMM.open(fileNameMM);
    for(int i = 0; i < v_MM.size(); i++)
      fileMM << dt*i << " " << v_MM[i].x << " " << v_MM[i].z << std::endl;
    fileMM.close();

    // Writte MB
    std::ofstream fileMB;
    std::string fileNameMB = fileName + "MB.dat";
    fileMB.open(fileNameMB);
    for(int i = 0; i < v_MB.size(); i++)
      fileMB << dt*i << " " << v_MB[i].x << " " << v_MB[i].z << std::endl;
    fileMB.close();

    // Writte RT
    std::ofstream fileRT;
    std::string fileNameRT = fileName + "RT.dat";
    fileRT.open(fileNameRT);
    for(int i = 0; i < v_RT.size(); i++)
      fileRT << dt*i << " " << v_RT[i].x << " " << v_RT[i].z << std::endl;
    fileRT.close();

    // Writte RM
    std::ofstream fileRM;
    std::string fileNameRM = fileName + "RM.dat";
    fileRM.open(fileNameRM);
    for(int i = 0; i < v_RM.size(); i++)
      fileRM << dt*i << " " << v_RM[i].x << " " << v_RM[i].z << std::endl;
    fileRM.close();

    // Writte RB
    std::ofstream fileRB;
    std::string fileNameRB = fileName + "RB.dat";
    fileRB.open(fileNameRB);
    for(int i = 0; i < v_RB.size(); i++)
      fileRB << dt*i << " " << v_RB[i].x << " " << v_RB[i].z << std::endl;
    fileRB.close();
  }

};

//-----------------------------------------------------------------------------------//
