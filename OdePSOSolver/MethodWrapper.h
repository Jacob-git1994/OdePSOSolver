#pragma once

#include "Common.h"
#include "Euler.h"
#include "OdeSolverParameters.h"
#include "Particle.h"
#include "ParticleProcessor.h"
#include "RungeKutta2.h"
#include "RungeKutta4.h"
#include "SolverIF.h"

#include <Eigen/Dense>
#include <exception>
#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <utility>

namespace odepso
{
	class MethodWrapper
	{

	private:

		//Hold our methods
		std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers> methodsMap;

		//Hold our parameters
		std::map<Common::SOLVER_TYPES, OdeSolverParameters, Common::CompSolvers> paramsMap;

		//Hold our results
		std::map<Common::SOLVER_TYPES, std::vector<std::pair<OdeSolverParameters, Eigen::VectorXd>>, Common::CompSolvers> resultsMap;

		//Our particle processor
		std::map<Common::SOLVER_TYPES, ParticleProcessor, Common::CompSolvers> psoProcessorMap;

		//Store some of our random variables
		struct RandomStruct
		{
			//Our inital seed
			long long seed;

			//Our random engine
			std::minstd_rand randomEngine;
		};

		RandomStruct randomStruct;

	public:

		MethodWrapper();

		MethodWrapper(const MethodWrapper& methodWrapperIn) = default;

		MethodWrapper& operator=(const MethodWrapper& methodWrapperIn) = default;

		virtual ~MethodWrapper() = default;

		//Initalize everything. Return false if anything fails
		bool initialize(const OdeSolverParameters& paramsIn);

		//Get the methods map
		std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>& getMethodsMap();

		//Get the params map
		std::map<Common::SOLVER_TYPES, OdeSolverParameters, Common::CompSolvers>& getParamsMap();

		//Get the results map
		std::map<Common::SOLVER_TYPES, std::vector<std::pair<OdeSolverParameters, Eigen::VectorXd>>, Common::CompSolvers>& getResultsMap();

		//Get the PSO Processor Map
		std::map<Common::SOLVER_TYPES, ParticleProcessor, Common::CompSolvers>& getPSOMap();

		//Get randomGenerator
		RandomStruct& getRandStruct();

		//Get the methods map
		const std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>& getMethodsMap() const;

		//Get the params map
		const std::map<Common::SOLVER_TYPES, OdeSolverParameters, Common::CompSolvers>& getParamsMap() const ;

		//Get the results map
		const std::map<Common::SOLVER_TYPES, std::vector<std::pair<OdeSolverParameters, Eigen::VectorXd>>, Common::CompSolvers>& getResultsMap() const;

		//Get the PSO Processor Map
		const std::map<Common::SOLVER_TYPES, ParticleProcessor, Common::CompSolvers>& getPSOMap() const;

		//Get randomGenerator
		const RandomStruct& getRandStruct() const;
	};
}
