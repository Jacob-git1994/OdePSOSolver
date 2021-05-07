#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <list>
#include <stdexcept>

#include "Common.h"
#include "SolverIF.h"

namespace odepso
{
	//Class to impliment RK4
	class RungeKutta4 : public SolverIF
	{
	private:

		//Helper Solver Vectors
		Eigen::VectorXd k1, k2, k3, k4;

	public:

		//Default constructor
		RungeKutta4() = default;

		//Default destructor
		virtual ~RungeKutta4() = default;

		//Delete the copy constructor
		RungeKutta4(const RungeKutta4& rungeKutta4In);

		//Delete the assign constructor
		RungeKutta4& operator=(const RungeKutta4& rungeKutta4In);

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