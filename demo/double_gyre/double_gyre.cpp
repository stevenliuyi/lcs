#include <iostream>
#include "../../ftle.hpp"
#include "../../velocity_function.hpp"
#include "../../io.hpp"

using namespace LCS;

int main()
{
    ContinuousFlowField<double,VelocityFunction::DoubleGyreModel<double>,2> double_gyre(1000,500);
    double_gyre.InitialPosition().SetAll(0,2,0,1);

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


