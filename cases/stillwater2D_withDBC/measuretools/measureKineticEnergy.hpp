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
    EkinTot += 0.5*dot(V.v[i], V.v[i]);
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

  std::vector<Vector3d> r_init;
  std::vector<Vector3d> delta_r;

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
      //EkinTot_local += 0.5*dot(V.v[i], V.v[i]);
    }
    EkinTot.push_back(EkinTot_local);
    std::cout << EkinTot_local << std::endl;
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

  void TrackParticles(Variables<double, double> &V,
                      int LT, int LM, int LB,
                      int MT, int MM, int MB,
                      int RT, int RM, int RB)
  {
    r_LT.push_back(V.r[LT]); v_LT.push_back(V.v[LT]);
    r_LM.push_back(V.r[LM]); v_LM.push_back(V.v[LM]);
    r_LB.push_back(V.r[LB]); v_LB.push_back(V.v[LB]);

    r_MT.push_back(V.r[MT]); v_MT.push_back(V.v[MT]);
    r_MM.push_back(V.r[MM]); v_MM.push_back(V.v[MM]);
    r_MB.push_back(V.r[MB]); v_MB.push_back(V.v[MB]);

    r_RT.push_back(V.r[RT]); v_RT.push_back(V.v[RT]);
    r_RM.push_back(V.r[RM]); v_RM.push_back(V.v[RM]);
    r_RB.push_back(V.r[RB]); v_RB.push_back(V.v[RB]);
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
