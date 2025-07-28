#include "simulation.h"
#include "tissue.h"

#include <string>

Simulation::Simulation(Tissue* T)
{
    T_ = T;
}

void Simulation::Run(int n_timesteps, int ouput_interval)
{
    for (int t_step = 1; t_step <= n_timesteps; t_step++)
    {
        T_->RunTimestep();
        if (t_step % ouput_interval == 0) 
        {
            T_->WriteCellFile("cells" + std::to_string(t_step) + ".vtk");
            T_->WriteDirectorFile("directors" + std::to_string(t_step) + ".vtk");
        }
        
    }
}
