#pragma once

#include <algorithm>
#include <cmath>
#include <Eigen\Dense>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>
#include <utility>

#include "MethodWrapper.h"
#include "SolverIF.h"
#include "OdeSolverParameters.h"
#include "ParticleProcessor.h"
#include "Particle.h"
#include "ProblemWrapperIF.h"

namespace odepso
{
	class OdePSOEngine
	{
	private:

		MethodWrapper methodWrapper;

		double wayPointDt;
		double wayPointStartTime;
		double wayPointEndTime;

		size_t wayPoints;

		void runMethod(std::unique_ptr<SolverIF>& solverIn, 
			OdeSolverParameters& paramsIn, 
			ParticleProcessor& particleProcessorIn, 
			std::vector<std::pair<OdeSolverParameters,
			Eigen::VectorXd>>& resultVector, 
			const Eigen::VectorXd& initialConditions,
			const bool& bypassPSOIn);

	public:

		OdePSOEngine();

		OdePSOEngine(const OdePSOEngine& odePSOEngineIn) = delete;

		OdePSOEngine& operator=(const OdePSOEngine& odePSOEngineIn) = delete;

		~OdePSOEngine() = default;

		//Run
		void run(std::shared_ptr<ProblemWrapperIF> problemIn, const Eigen::VectorXd& initialConditions, const OdeSolverParameters& paramsIn, const double& leftTime, const double& rightTime, const size_t& numWayPoints, const bool& bypassPSO);

		//Get results of the best method
		const std::vector < std::pair<OdeSolverParameters, Eigen::VectorXd>>& getResults() const;

		//Get results for a paticular method
		const std::vector< std::pair<OdeSolverParameters, Eigen::VectorXd>>& getResults(Common::SOLVER_TYPES solverTypeIn) const;

	};
}