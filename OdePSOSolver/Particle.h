#pragma once

#include "OdeSolverParameters.h"
#include "Richardson.h"
#include "SolverIF.h"

#include <cmath>
#include <chrono>
#include <Eigen\Dense>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

namespace odepso
{
	class Particle
	{
	private:

		//Ode Parameters
		OdeSolverParameters particleParameters;

		//Richardson Parameters
		Richardson particleTables;

		//Best Result
		Eigen::VectorXd bestResult;

		//Cost values
		double errorCost; //Find error in range
		double timeCost; //Find smallest runtime for the error in range

		//Store the error
		std::vector<double> errorVector;

		//Parameter Velocities
		double dtVelocity;
		double richardsonVelocity;

	public:

		//Constructor
		Particle();

		//Our constructor
		Particle(const OdeSolverParameters& odeSolverParametersIn, const Richardson& richardsonIn);

		//Copy Constructor
		Particle(const Particle& particleIn);

		//Copy assisgn constructor
		Particle& operator=(const Particle& particleIn);

		//Destructor
		virtual ~Particle() = default;

		//Initalize all variables
		void initalize(std::minstd_rand& randEngineIn);

		//Update time step using new information
		void updateTimeStep(std::unique_ptr<SolverIF>& solverTypeIn, OdeSolverParameters& odeSolverParametersIn, const Eigen::VectorXd& currentStateIn, const double& finalTimeIn);

		//Update time step using new information
		void updateTimeStep(std::unique_ptr<SolverIF>& solverTypeIn, const Eigen::VectorXd& currentStateIn, const double& finalTimeIn);

		//Get the best result
		const Eigen::VectorXd& getBestResult() const;

		//Get the particle parameters
		const OdeSolverParameters& getParameters() const;

		//Get the error cost
		const double& getErrorCost() const;

		//Get the time cost
		const double& getTimeCost() const;

		//Estimate global error
		const double estimateGlobalError(const double& finalTimeIn) const;

		//Set the parameters
		void setParticleParameters(const OdeSolverParameters& odeSolverParametersIn);

		//Get dtVeclocity
		const double& getDtVelocity() const;

		//Set dtVelocity
		void setDtVelocity(const double& dtVelocityIn);

		//Get Richardson Velocity
		const double& getRichardsonVelocity() const;

		//Set Richardson Velocity
		void setRichardsonVelocity(const double& richardsonVelocityIn);

		//Clear particle history (essentially resets the particle in prep to be run again)
		void clearParticleHistory(const OdeSolverParameters& paramsIn);

		void setBestResult(const Eigen::VectorXd& bestStateIn);
	};
}
