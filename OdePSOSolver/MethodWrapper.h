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
		std::map<Common::SOLVER_TYPES, std::vector<std::pair<double, Eigen::VectorXd>>, Common::CompSolvers> resultsMap;

		//Our particle processor
		ParticleProcessor psoProcessor;

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
	};
}
