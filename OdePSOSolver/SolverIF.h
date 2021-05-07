#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <list>
#include <stdexcept>
#include <utility>

#include "ProblemWrapperIF.h"

namespace odepso
{
	//Class to hold common Solvers. Designed to only hold the method.
	class SolverIF
	{
	protected:

		//The wrapper to the problem we are trying to solve
		std::shared_ptr<ProblemWrapperIF> problem;

		//Generalize and initalize what solver vectors are needed
		void initalizeVectors(std::list<Eigen::VectorXd*>& solverHelperVectorsIn, const size_t& vecSizeIn);

	public:

		//Default constructor
		SolverIF() = default;

		//Copy Constructor
		SolverIF(const SolverIF& solverIFIn);

		//Assignment Constructor
		SolverIF& operator=(const SolverIF& solverIFIn);

		//Default destructor
		virtual ~SolverIF() = default;

		//Save the problem we are trying to solve
		void initalizeProblem(std::shared_ptr<ProblemWrapperIF> problemIn);

		//This will call initalizeVectors
		virtual void initalizeSolverVectors(const size_t& vecSizeIn) = 0;

		//Update the method forward in time a fixed number of steps
		virtual void updateTimeStep(
			const Eigen::VectorXd& initalConditionsIn,
			Eigen::VectorXd& updateStateIn,
			const double& dtIn,
			const unsigned int& numberOfStepsIn,
			const double& currentTimeIn) = 0;

		//Calculate one step propagating the derivative forward in time by a set dt
		virtual void updateTimeStep(
			const Eigen::VectorXd& initalConditionsIn,
			Eigen::VectorXd& updateStateIn,
			const double& dtIn,
			const double& currentTimeIn) = 0;


		//Get the theoretical convergence of the method
		virtual const double getTheoreticalConvergence() const = 0;

		//Deep copy of this pointer
		virtual SolverIF* deepClone() const = 0;
	};
}


