/**
    @file discrete_double_gyre.cpp
    @brief Demo for FTLE calculation of discrete double-gyre model.

    This file contains a demo for performing FTLE calculation of a discrete double-gyre model. Both positive and negative FLTE are calculated here.
    */

#include <iostream>
#include "../../src/ftle.hpp"
#include "../../src/velocity_function.hpp"
#include "../../src/io.hpp"

using namespace LCS;

int main()
{
    // number of threads for OpenMP
    std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;

    Position<double,2> pos(100,50);
    pos.SetAll(0,2,0,1);

    ContinuousVelocity<double, VelocityFunction::DoubleGyreModel<double>,2> double_gyre_vel(100,50,pos);

    // write discrete data to file
    for (int t=0; t <= 20; ++t)
    {
        std::stringstream ss;
        ss << "double_gyre_" << t << ".txt";
        double_gyre_vel.UpdateTime(t);
        double_gyre_vel.SetAll();
        double_gyre_vel.WriteToFile(ss.str());
    }
    std::cout << "Discrete data written to files" << std::endl;

    DiscreteFlowField<double,2> double_gyre(1000,500,100,50);
    
    double_gyre.DataPosition().SetAll(0,2,0,1);
    double_gyre.InitialPosition().SetAll(0,2,0,1);

    double_gyre.SetVelocityFileNamePrefix("double_gyre_");

    double_gyre.SetDataDelta(1);
    double_gyre.SetDataTimeRange(0,20);
    double_gyre.SetDelta(.1);
    double_gyre.SetStep(200);
    double_gyre.Run();

    // positive FTLE
    FTLE<double,2> ftle(double_gyre);
    ftle.Calculate();
    ftle.WriteToFile("double_gyre_ftle_pos.txt");

    double_gyre.SetDirection(Direction::Backward);
    double_gyre.SetInitialTime(20);
    double_gyre.Run();

    // negative FTLE
    ftle.Calculate();
    ftle.WriteToFile("double_gyre_ftle_neg.txt");

    return 0;
}
