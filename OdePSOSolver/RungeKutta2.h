#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <list>
#include <stdexcept>

#include "Common.h"
#include "SolverIF.h"

namespace odepso
{
    class RungeKutta2 : public SolverIF
    {
    private:

        //Helper Solver Vectors
        Eigen::VectorXd k1, k2;

    public:

		//Default constructor
		RungeKutta2() = default;

		//Default destructor
		virtual ~RungeKutta2() = default;

		//Copy constructor
		RungeKutta2(const RungeKutta2& rungeKutta2In);

		//Assignment constructor
		RungeKutta2& operator=(const RungeKutta2& rungeKutta2In);

		//Initalize our vectors
		virtual void initalizeSolverVectors(const size_t& vecSizeIn) override;

		//Update the method forward in time a fixed number of steps
		virtual void updateTimeStep(
			const Eigen::VectorXd& initalConditionsIn,
			Eigen::VectorXd& updateStateIn,
			const double& dtIn,
			const unsigned int& numberOfStepsIn,
			const double& currentTimeIn) override;

		//Calculate one step propagating the derivative forward in time by a set dt
		virtual void updateTimeStep(
			const Eigen::VectorXd& initalConditionsIn,
			Eigen::VectorXd& updateStateIn,
			const double& dtIn,
			const double& currentTimeIn) override;

		//Get the theoretical convergence of the method
		virtual const double getTheoreticalConvergence() const override;

		//Deep copy of this pointer
		virtual SolverIF* deepClone() const override;
    };
}