#include "OdePSOEngine.h"

void odepso::OdePSOEngine::runMethod(std::unique_ptr<SolverIF>& solverIn,
	OdeSolverParameters& paramsIn,
	ParticleProcessor& particleProcessorIn,
	std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>& resultVector,
	const Eigen::VectorXd& initialConditions, 
	const bool& bypassPSOIn)
{
	//Build up the way points
	double curTime = wayPointStartTime;
	double nextTime = 0.;

	//Should there be an optimization?
	bool bypassPSO = bypassPSOIn;

	//Bypass counter if PSO is never converging don't rerun multiple times in a row
	unsigned int bypassCounter = 0;

	//Current State
	Eigen::VectorXd tempState = initialConditions;

	//Current Particle
	odepso::Particle bestParticle = odepso::Particle();

	//Build a vector for the different way points
	std::vector<double> timePoints;
	timePoints.push_back(curTime);
	while (curTime + wayPointDt < wayPointEndTime - wayPointDt)
	{
		timePoints.push_back(curTime += wayPointDt);
	}
	timePoints.push_back(wayPointEndTime);

	//Time point iterator
	std::vector<double>::const_iterator timePointItr = timePoints.cbegin();

	//Our the first result (initial condition)
	std::pair<odepso::OdeSolverParameters, Eigen::VectorXd> tempResult(paramsIn, tempState);
	resultVector.push_back(tempResult);

	//If bypass pso set the fields the max dt and min richardson tables
	if (bypassPSOIn)
	{
		//Set our fields
		paramsIn.setDt(paramsIn.getMaxDt());
		paramsIn.setRichardsonLevels(paramsIn.getMinRichardsonTables());

		//Set our particle
		bestParticle.setParticleParameters(paramsIn);
	}

	//Go through our waypoints
	while (timePointItr != std::prev(timePoints.cend(), 1))
	{
		curTime = *timePointItr;
		nextTime = *(std::next(timePointItr));
		
		std::cout << curTime << "\t" << nextTime << "\n";

		//Perform PSO
		if (!bypassPSO)
		{
			//Set our parameters
			particleProcessorIn.setParamsAndWayPoint(paramsIn, curTime, nextTime, wayPointEndTime);

			//Perform PSO optimization on our problem
			particleProcessorIn.PSO(solverIn, methodWrapper.getRandStruct().randomEngine, tempState);

			//Update the best particle
			bestParticle = particleProcessorIn.getBestParticle();

			//Update the bypass counter
			bypassCounter++;
		}
		else
		{
			//Set our parameters
			particleProcessorIn.setParamsAndWayPoint(paramsIn, curTime, nextTime, wayPointEndTime);

			//Bypass the full PSO
			particleProcessorIn.runParticle(solverIn, bestParticle, tempState);

			//Reset the bypass counter if we sucessfully performed pso - otherwise we give up as the run time gets too large
			if (bypassCounter++ != 3 && bypassPSO)
			{
				bypassCounter = 0;
			}
		}

		//Check if PSO needs to be performed another time if the error is not less than max allowed error.
		if ((bestParticle.estimateGlobalError(wayPointEndTime) <= paramsIn.getMaxAllowedError() &&
			bestParticle.estimateGlobalError(wayPointEndTime) >= paramsIn.getMinAllowedError()) ||
			bypassPSOIn) 
		{
			bypassPSO = true;

			//Save off our best particle
			paramsIn = bestParticle.getParameters();
			tempState = bestParticle.getBestResult();

			//Add our results
			std::pair<odepso::OdeSolverParameters, Eigen::VectorXd> tempResult(paramsIn, tempState);
			resultVector.push_back(tempResult);

			//Update the current time
			timePointItr++;
		}
		else if (bestParticle.estimateGlobalError(wayPointEndTime) <= paramsIn.getMaxAllowedError() && bypassCounter >= 3)
		{
			//We are finding a non-preferred error, but a valid result so no need to run PSO
			bypassPSO = true;

			//Save off our best particle
			paramsIn = bestParticle.getParameters();
			tempState = bestParticle.getBestResult();

			//Add our current results
			std::pair<odepso::OdeSolverParameters, Eigen::VectorXd> tempResult(paramsIn, tempState);
			resultVector.push_back(tempResult);

			//Update the current time
			timePointItr++;
		}
		else if (!bypassPSO)
		{
			//We are finding a non-preferred error, but a valid result so no need to run PSO
			bypassPSO = true;

			//Save off our best particle
			paramsIn = bestParticle.getParameters();
			tempState = bestParticle.getBestResult();

			//Add our results
			std::pair<odepso::OdeSolverParameters, Eigen::VectorXd> tempResult(paramsIn, tempState);
			resultVector.push_back(tempResult);

			//Update the current time
			timePointItr++;
		}
		else
		{
			//We are failing to find a valid pso solution or we want to bypass
			bypassPSO = (false || bypassPSOIn) || bypassCounter >= 3;
		}

		//std::cout << curTime << "\t" << paramsIn.getCurrentTime() << "\n";
	}
}

odepso::OdePSOEngine::OdePSOEngine() :
	methodWrapper(odepso::MethodWrapper()),
	wayPointDt(0.),
	wayPointStartTime(0.),
	wayPointEndTime(0.),
	wayPoints(1)
{
	//Nothing else to do here
}

void odepso::OdePSOEngine::run(std::shared_ptr<ProblemWrapperIF> problemIn, const Eigen::VectorXd& initialConditions, const OdeSolverParameters& paramsIn, const double& leftTime, const double& rightTime, const size_t& numWayPoints, const bool& bypassPSO)
{
	if (!methodWrapper.initialize(paramsIn) || leftTime > rightTime)
	{
		throw std::runtime_error("Failed to initialize methods");
	}

	//Initialize the methods and set the initial time
	for (auto& methodItr : methodWrapper.getMethodsMap())
	{
		methodItr.second->initalizeProblem(problemIn);
		methodItr.second->initalizeSolverVectors(initialConditions.size());
		methodWrapper.getParamsMap()[methodItr.first].setCurrentTime(leftTime);
	}

	//Set our way point dt
	wayPointDt = (rightTime - leftTime) / static_cast<double>(numWayPoints);
	wayPointStartTime = leftTime;
	wayPointEndTime = rightTime;

	//Our thread vector
	std::vector<std::thread> threadVec;

	//Go over all our methods
	for (auto& methodItr : methodWrapper.getMethodsMap())
	{
		//Get our approprate params and solvers
		std::unique_ptr<SolverIF>& curMethod = methodItr.second;
		odepso::OdeSolverParameters& curParams = methodWrapper.getParamsMap()[methodItr.first];
		odepso::ParticleProcessor& curPartProcessor = methodWrapper.getPSOMap()[methodItr.first];
		std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>& curResult = methodWrapper.getResultsMap()[methodItr.first];

		//Push back our results
		//threadVec.push_back(std::thread(&OdePSOEngine::runMethod, this, std::ref(curMethod), std::ref(curParams), std::ref(curPartProcessor), std::ref(curResult), std::cref(initialConditions), std::cref(bypassPSO)));
		runMethod(curMethod, curParams, curPartProcessor, curResult, initialConditions, bypassPSO);
	}

	//Join our threads
	for (auto& thItr : threadVec)
	{
		thItr.join();
	}
}

const std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>& odepso::OdePSOEngine::getResults() const
{
	const std::map<Common::SOLVER_TYPES, std::vector<std::pair<OdeSolverParameters, Eigen::VectorXd>>, Common::CompSolvers>& resultsMap = methodWrapper.getResultsMap();
	const std::vector<std::pair<OdeSolverParameters, Eigen::VectorXd>>* bestResult = nullptr;

	//Go through our results
	double error = 0.;
	double minError = std::numeric_limits<double>::infinity();
	for (const auto& resultItr : resultsMap)
	{
		error = resultItr.second.back().first.getTotalError();
		if (error < minError)
		{
			bestResult = &resultItr.second;
			minError = error;
		}
	}

	return *bestResult;
}
