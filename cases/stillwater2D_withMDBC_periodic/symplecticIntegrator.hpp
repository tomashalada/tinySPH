// Symplectic integrator

/*
SymplecticScheme WCSPHSymplecticFluid(&WCSPHfluid, 0.000035);
SymplecticScheme WCSPHSymplecticBound(&WCSPHbound, 0.000035);

for(int step = 0; step < stepEnd + 1; step++)
{

  std::cout << "STEP: " << step << std::endl;
  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  WCSPHSymplecticFluid.ComputePredictor();
  WCSPHSymplecticBound.ComputePredictor();
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);

  DensityToPressure(WCSPHfluid, WCSPHconstants);
  DensityToPressure(WCSPHbound, WCSPHconstants);

  WCSPHSymplecticFluid.ComputeCorrector();
  WCSPHSymplecticBound.ComputeCorrector();
  WCSPH.Interact(WCSPHfluid, WCSPHbound, WCSPHconstants);

  if(step % saveOutput == 0){
    writeParticleData(WCSPHfluid, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/FLUID/fluid", step));
    writeParticleData(WCSPHbound, stepToNameWithPtcsExtension(caseResults + "/OUTPUT/BOUND/bound", step));
  }

}
*/

