#include "simulation.h"
#include "tissue.h"

#include <string>

Simulation::Simulation(Tissue* T)
{
    T_ = T;
}

void Simulation::run(int n_timesteps, int ouput_interval)
{
    for (int t_step = 1; t_step <= n_timesteps; t_step++)
    {
        T_->runTimestep();
        if (t_step % ouput_interval == 0) 
        {
            T_->writeCellFile("cells" + std::to_string(t_step) + ".vtk");
            T_->writeDirectorFile("directors" + std::to_string(t_step) + ".vtk");
        }
        
    }
}
