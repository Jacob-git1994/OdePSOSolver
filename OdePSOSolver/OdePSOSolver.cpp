// OdePSOSolver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <utility>
#include <vector>
#include <Eigen\Dense>
#include <iomanip>
#include <chrono>
#include <memory>
#include <random>

#include "Common.h"
#include "Euler.h"
#include "RungeKutta4.h"
#include "TestExp.h"
#include "RungeKutta2.h"
#include "Richardson.h"
#include "Particle.h"
#include "ParticleProcessor.h"

using namespace odepso;

int main()
{
#if 0
	Eigen::VectorXd sol;
	Eigen::VectorXd ic;
	ic.resize(1); ic(0) = 1.;
	RungeKutta2 solver;
	solver.initalizeProblem(std::shared_ptr<ProblemWrapperIF>(new TestExp));
	solver.initalizeSolverVectors(1);
	solver.updateTimeStep(ic, sol, .00001, 100000, 0.);
	auto solv2 = solver;
	solv2.updateTimeStep(ic, sol, .00001, 100000, 0.);
	std::cout << sol;
#elif 0

	Richardson richardsonTables;
	Eigen::VectorXd sol;
	Eigen::VectorXd ic;
	ic.resize(1); ic(0) = 1.;
	sol = ic;
	RungeKutta4 solver;
	solver.initalizeProblem(std::shared_ptr<ProblemWrapperIF>(new TestExp));
	solver.initalizeSolverVectors(1);
	
	//Matrix to store data
	Eigen::MatrixXd resultMat;

	//dt ranges and rich levels
	const std::vector<double> dtRanges = { 1., .5,.1, .05, .01 };
	const std::vector<unsigned int> richLevels = { 2,3,4,5,6,7,8,9,10 };

	resultMat.resize(dtRanges.size(), richLevels.size());
	resultMat.fill(0.0);

	unsigned int k = 0, j = 0;

	for (const auto& dtEl : dtRanges)
	{
		const double& dt = dtEl;
		j = 0;

		for (const auto& richEl : richLevels)
		{
			auto t1 = std::chrono::high_resolution_clock::now();

			//Initalize richardson
			richardsonTables.reset(richEl, ic.size(), 2, dt, solver.getTheoreticalConvergence());

			//Time interval Splitting
			double currentTime = 0;
			double endTime = 1;
			double tempError = 0.0;
			while (currentTime < endTime)
			{
				Eigen::VectorXd temp = ic;
				//Initalize richardson
				richardsonTables.reset(richEl, ic.size(), 2, dt, solver.getTheoreticalConvergence());
				for (size_t i = 0; i < richardsonTables.getTableSize(); ++i)
				{
					solver.updateTimeStep(temp, sol, dt / std::pow(2., i), ((unsigned int)((currentTime + dt) / dt) * std::pow(2, i)), currentTime);
					richardsonTables.append(i, sol);
				}

				richardsonTables.calculateRichardsonExtrapolation();
				//std::cout << std::setprecision(16) << richardsonTables.getBestResult() << "\n" << std::setprecision(16) << richardsonTables.getError() << "\n" << richardsonTables.getC() << "\n";

				currentTime += dt;
				temp = richardsonTables.getBestResult();
				tempError += richardsonTables.getError();
			}
			auto t2 = std::chrono::high_resolution_clock::now();
			double runTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
			resultMat(k, j++) = tempError ;//static_cast<double>((tempError >= 1e-8) && (tempError <= 1e-6))* tempError;//-std::log(((tempError >= 1e-15) && (tempError <= 1e-13)) * tempError); // *runTime;
		}
		++k;
	}

	std::cout << resultMat <<"\n\n" << "\n";
#elif 0
	std::unique_ptr<SolverIF> solver(new RungeKutta4());
	Richardson richardsonTables;
	Eigen::VectorXd sol;
	Eigen::VectorXd ic;
	ic.resize(1); ic(0) = 1.;
	sol = ic;

	solver->initalizeProblem(std::shared_ptr<ProblemWrapperIF>(new TestExp));
	solver->initalizeSolverVectors(1);

	ParticleParameters partParams(10, 5, 1, 1, 1, .9, 1);
	OdeSolverParameters params(1e-6, 1e-4, 5, 6, .01, .1, partParams, 2.);

	Particle myParticle (params, richardsonTables, solver);
	std::minstd_rand randThing(555);
	myParticle.initalize(randThing);
	Eigen::VectorXd currentState = ic;
	while (myParticle.getParticleParameters().getCurrentTime() < 1.)
	{
		myParticle.updateTimeStep(params, currentState, 1.);
		std::cout << std::setprecision(16) << 
			myParticle.getBestResult() << "\n" << 
			myParticle.getParticleParameters().getTotalError() << "\n" << 
			myParticle.getErrorCost() << "\n" <<
			myParticle.getTimeCost() << "\n" <<
			myParticle.getParticleParameters().getCurrentTime() << "\n" << 
			myParticle.getParticleParameters().getDt() << "\n" << 
			myParticle.getParticleParameters().getRichardsonLevels() << "\n" <<
			myParticle.estimateGlobalError(1.) << "\n\n";

		currentState = myParticle.getBestResult();
	}

#else 
	std::unique_ptr<SolverIF> solver(new Euler());
	Richardson richardsonTables;
	Eigen::VectorXd sol;
	Eigen::VectorXd ic;
	ic.resize(1); ic(0) = 1.;
	sol = ic; int z = 5;

	solver->initalizeProblem(std::shared_ptr<ProblemWrapperIF>(new TestExp));
	solver->initalizeSolverVectors(1);

	ParticleParameters partParams(10, 100, .5, 0., .5, 1, 1000, 1e-10);
	OdeSolverParameters params(1e-10, 1e-4, 2, 10, .01, 1, partParams, 4);

	std::minstd_rand randThing(654321);
	Eigen::VectorXd currentState = ic;

	ParticleProcessor processor;

	processor.setParamsAndWayPoint(params, 0, 1, 1.);

	processor.PSO(solver, randThing, currentState);

	std::cout << std::setprecision(16) << processor.getBestParticle().getParameters().getTotalError() << "\n" <<
		processor.getBestParticle().getBestResult() << "\n" <<
		processor.getBestParticle().getParameters().getDt() << "\n" <<
		processor.getBestParticle().getParameters().getRichardsonLevels() << "\n" <<
		processor.getVarDt() << "\n" <<
		processor.getExpectedDt() << "\n" <<
		processor.getVarRich() << "\n" <<
		processor.getParams().getC() << "\n" <<
		processor.getBestParticle().estimateGlobalError(1.) << "\n" <<
		processor.getBestParticle().getParameters().getTotalRunTime() << "\n";

#endif

}
