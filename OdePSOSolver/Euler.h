#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <list>
#include <stdexcept>

#include "SolverIF.h"

namespace odepso
{
	//Class to impliment Eulers method
	class Euler : public SolverIF
	{
	private:

		//Helper Solver Vectors
		Eigen::VectorXd k1;

	public:

		//Default constructor
		Euler() = default;

		//Default destructor
		virtual ~Euler() = default;

		//Delete the copy constructor
		Euler(const Euler& eulerIn);

		//Delete the assign constructor
		Euler& operator=(const Euler& eulerIn);

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