#include <CompactNSearch>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

#include <omp.h>
#include <algorithm>


#include <string>
#include <string_view>
#include <cstdlib>

#include "vector.hpp"
#include "variables.hpp"
#include "writeParticleData.hpp"

#include "interaction.hpp"
#include "kernel.hpp"
#include "integration.hpp"

#include "generalTools.hpp"



// ------------------------------------------------------------------------

int main(){

	//Set the case
	ConstantVariables WCSPHconstants(
			0.0145, //h
			0.1, //m
			40.0, //cs
			0.02, //eta
			1000. //rho0
			);

	//Read fluid data
	FluidVars<double, double> WCSPHfluid;
	WCSPHfluid.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH/cases/dambreak2D_withDBC/dambreak_fluid.ptcs");
				std::cout << "Number of fluid particles: " << WCSPHfluid.N << std::endl;

	//Read boundary data
	BoundVars<double, double> WCSPHbound;
	WCSPHbound.initializeWithGeometryFile("/home/tomas/Documents/testovaci/cpp/TINYSPH/cases/dambreak2D_withDBC/dambreak_wall.ptcs");
				std::cout << "Number of boundary particles: " << WCSPHbound.N << std::endl;

	double h = 0.0145;

	VerletScheme WCSPHVerlet(&WCSPHfluid, 0.00005);
	VerletScheme WCSPHVerletProv(&WCSPHbound, 0.00005);

	//Init CompactNSearch (neighbors searchgin)
	CompactNSearch::NeighborhoodSearch nsearch(h, true);
	//unsigned int positions_id = nsearch.add_point_set(positions.front().data(), positions.size());
	unsigned int positions_id = nsearch.add_point_set(&WCSPHfluid.r.front().x, WCSPHfluid.N);
	unsigned int positions_id_wall = nsearch.add_point_set(&WCSPHbound.r.front().x, WCSPHbound.N);


// ------------------------------------------------------------------------

int STEP_MAX = 2000;
for(int STEP = 0; STEP < STEP_MAX; STEP++)
{
	std::cout << "STEP: " << STEP << std::endl;

	nsearch.find_neighbors();
	//#pragma omp parallel
	//{
 //   printf("Hello World... from thread = %d\n",
 //          omp_get_thread_num());
	//}

	DensityToPressure(WCSPHfluid, WCSPHconstants);
	DensityToPressure(WCSPHbound, WCSPHconstants);

	CompactNSearch::PointSet const& ps_1 = nsearch.point_set(positions_id);
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < ps_1.n_points(); ++i)
	{
	    // Get point set 1 neighbors of point set 1.
					InteractionData<double, double> ptcI;
					InteractionData<double, double> ptcJ;
					ptcI.CopyDataIn(WCSPHfluid, i);
	    for (size_t j = 0; j < ps_1.n_neighbors(positions_id, i); ++j)
	    {
	        // Return the point id of the jth neighbor of the ith particle in the point_set_1.
	        const unsigned int pid = ps_1.neighbor(positions_id, i, j);
									ptcJ.CopyDataIn(WCSPHfluid, pid);
									WCSPH_FluidFluidInteraction<double, double, WendlandKernel>(ptcI, ptcJ, WCSPHconstants);
	    }

	    // Get point set 1 neighbors of point set 2, i.e. to boundary particles.
	    for (size_t j = 0; j < ps_1.n_neighbors(positions_id_wall, i); ++j)
	    {
	        // Return the point id of the jth neighbor of the ith particle in the point_set_1.
	        const unsigned int pid = ps_1.neighbor(positions_id_wall, i, j);
									ptcJ.CopyDataIn(WCSPHbound, pid);
									WCSPH_FluidBoundInteraction<double, double, WendlandKernel>(ptcI, ptcJ, WCSPHconstants);
									//std::cout << "Pair for particle i: " << i << " with particle j: " << pid << std::endl;

									//std::cout << "Pair for particle i: " << i << " with particle j: " << pid << " total #nbs:" << ps_1.n_neighbors(positions_id_wall, i) << std::endl;
									//std::cout << "Particle i (fluid) position: " << coutPOS2D(WCSPHfluid.r[i]) << " particle j (wall) position: " << coutPOS2D(WCSPHbound.r[pid]) << std::endl;
									//std::cout << std::endl;
	    }
					ptcI.a.z -= 9.81; //apply external forces
			 	ptcI.CopyDataOut(WCSPHfluid, i);
	}

	CompactNSearch::PointSet const& ps_2 = nsearch.point_set(positions_id_wall);
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < ps_2.n_points(); ++i)
	{

					InteractionData<double, double> ptcI;
					InteractionData<double, double> ptcJ;
					ptcI.CopyDataIn(WCSPHbound, i);
	    for (size_t j = 0; j < ps_2.n_neighbors(positions_id, i); ++j)
	    {
	        // Return the point id of the jth neighbor of the ith particle in the point_set_1.
	        const unsigned int pid = ps_2.neighbor(positions_id, i, j);
									ptcJ.CopyDataIn(WCSPHfluid, pid);
									WCSPH_UpdateBoundaryData<double, double, WendlandKernel>(ptcI, ptcJ, WCSPHconstants);
									//WCSPH_UpdateBoundaryDBCsimple<double, double, WendlandKernel>(ptcI, ptcJ, WCSPHconstants);
									//std::cout << "Pair for particle i: " << i << " with particle j: " << pid << std::endl;

									//std::cout << "Pair for particle i: " << i << " with particle j: " << pid << " total #nbs:" << ps_2.n_neighbors(positions_id, i) << std::endl;
									//std::cout << "Particle i (wall) position: " << coutPOS2D(WCSPHbound.r[i]) << " particle j (fluid) position: " << coutPOS2D(WCSPHfluid.r[pid]) << std::endl;
									//std::cout << std::endl;
	    }

	 //ptcI.CopyDataOutBound(WCSPHbound, i);
	 ptcI.CopyDataOut(WCSPHbound, i);

	}

	if(STEP % 20 == 0){
	WCSPHVerlet.ComputeStepEulerForStability();
	WCSPHVerletProv.ComputeStepEulerForStability();
	}
	else{
	WCSPHVerlet.ComputeStep();
	WCSPHVerletProv.ComputeStep();
	}

}

// ------------------------------------------------------------------------

	writeParticleData(WCSPHfluid, "__test-output-fluid.ptcs");
	writeParticleData(WCSPHbound, "__test-output-bound.ptcs");



	std::cout << "working..." << std::endl;

	return EXIT_SUCCESS;

}
