#ifndef SIMULATION_H
#define SIMULATION_H

class Tissue;

class Simulation
{
private:

    Tissue* T_;

public:

    Simulation(Tissue* T);
    void run(int n_timesteps, int output_interval);

};

#endif // SIMULATION_H
