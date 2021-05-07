#pragma once

#include "OdeSolverParameters.h"
#include "Particle.h"
#include "Richardson.h"
#include "SolverIF.h"

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

namespace odepso
{
	class ParticleProcessor
	{
	private:

		//Best parameters
		OdeSolverParameters bestParameters;

		//Current Particles and their results
		std::vector<Particle> particles;

		//Best Particles states of particles
		std::vector<Particle> bestLocalParticles;

		//Current way points
		double wayPointBegin;
		double wayPointEnd;

		//Best particle
		Particle bestParticle;

		//Run one particle until end time
		void runParticle(std::unique_ptr<SolverIF>& solverTypeIn, Particle& particleIn, const Eigen::VectorXd& currentStateIn);

		//Compare particles
		bool compParticle(const Particle& lhs, const Particle& rhs) const;

		//Check if error is in range
		bool globalErrorInRange(const Particle& lhs) const;

		//Check if error is less than desired tol
		bool isLessThanMaxError(const Particle& particleIn) const;

		//Generate and set particle
		Particle generateAndSetParticle(std::unique_ptr<SolverIF>& solverTypeIn, const Particle& particleIn, Particle& bestLocalParticleIn, std::minstd_rand& randIn, const Eigen::VectorXd& currentStateIn);

		//Ensure all parameters are valid
		void ensureAllParametersAreValid(OdeSolverParameters& paramsIn);

		//Particle Statistics
		double expectedDt;
		double expectedRichardson;
		double varDt;
		double varRich;

		//Calculate Stats
		void calcStatistics();

		double finalEndTime;

	public:

		//Constructor
		ParticleProcessor();

		//Delete Copy
		ParticleProcessor(const ParticleProcessor& pararticleProcessorIn) = delete;

		//Delete Assign copy
		ParticleProcessor& operator=(const ParticleProcessor& pararticleProcessorIn) = delete;

		virtual ~ParticleProcessor() = default;

		//Set best parameters
		void setParamsAndWayPoint(const OdeSolverParameters& odeSolverParametersIn, const double& wayPointBeginIn, const double& wayPointEndIn, const double& finalEndTimeIn);

		//Get best parameters
		const OdeSolverParameters& getParams() const;

		//Get the best particle
		const Particle& getBestParticle() const;

		//PSO Optimization
		void PSO(std::unique_ptr<SolverIF>& solverTypeIn, std::minstd_rand& randIn, const Eigen::VectorXd& currentStateIn);

		//Get the expected dt
		const double& getExpectedDt() const;

		//Get expected rich table size
		const double& getExpectedRich() const;

		//Get dt var
		const double& getVarDt() const;

		//Get rich var
		const double& getVarRich() const;
	};
}