#include "MethodWrapper.h"

odepso::MethodWrapper::MethodWrapper() :
	methodsMap(std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>()),
	paramsMap(std::map<Common::SOLVER_TYPES, OdeSolverParameters, Common::CompSolvers>()),
	resultsMap(std::map<Common::SOLVER_TYPES, std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>, Common::CompSolvers>()),
	psoProcessorMap(std::map<Common::SOLVER_TYPES, ParticleProcessor, Common::CompSolvers>()),
	randomStruct(RandomStruct())
{
	//Nothing else to do here
}

bool odepso::MethodWrapper::initialize(const OdeSolverParameters& paramsIn)
{
	//Make sure everything is valid
	if (!paramsIn.paramsValid())
	{
		return false;
	}

	//If valid clear our lists
	methodsMap.clear();
	paramsMap.clear();
	psoProcessorMap.clear();
	resultsMap.clear();
	randomStruct = RandomStruct();

	try
	{
		//Build Euler
		if (paramsIn.getAllowedEuler() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::EULER, std::move(std::unique_ptr<odepso::SolverIF>(new Euler())));
		}

		//Build RK2
		if (paramsIn.getAllowedRK2() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::RK2, std::move(std::unique_ptr<odepso::SolverIF>(new RungeKutta2())));
		}

		//Build RK4
		if (paramsIn.getAllowedRK4() && !paramsIn.getIsSystemStiff())
		{
			methodsMap.emplace(Common::SOLVER_TYPES::RK4, std::move(std::unique_ptr<odepso::SolverIF>(new RungeKutta4())));
		}

		//Build our remaining maps
		for (std::map<Common::SOLVER_TYPES, std::unique_ptr<SolverIF>, Common::CompSolvers>::const_iterator methodItr = methodsMap.cbegin();
			methodItr != methodsMap.cend(); ++methodItr)
		{
			//Build parameter map
			paramsMap.emplace(methodItr->first, paramsIn);

			//Build result map
			resultsMap.emplace(methodItr->first, std::move(std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>()));

			//Build our particle processor map
			psoProcessorMap.emplace(methodItr->first, std::move(odepso::ParticleProcessor()));
		}

		//Set up our random generators
		randomStruct.seed = paramsIn.getParticleParameters().getSeed();
		randomStruct.randomEngine = std::minstd_rand(static_cast<unsigned int>(randomStruct.seed));
	}
	catch (std::exception& e)
	{
		std::cerr << e.what();
		return false;
	}

	//If everything passed return true
	return true;
}

std::map<odepso::Common::SOLVER_TYPES, std::unique_ptr<odepso::SolverIF>, odepso::Common::CompSolvers>& odepso::MethodWrapper::getMethodsMap()
{
	return methodsMap;
}

std::map<odepso::Common::SOLVER_TYPES, odepso::OdeSolverParameters, odepso::Common::CompSolvers>& odepso::MethodWrapper::getParamsMap()
{
	return paramsMap;
}

std::map<odepso::Common::SOLVER_TYPES, std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>, odepso::Common::CompSolvers>& odepso::MethodWrapper::getResultsMap()
{
	return resultsMap;
}

std::map<odepso::Common::SOLVER_TYPES, odepso::ParticleProcessor, odepso::Common::CompSolvers>& odepso::MethodWrapper::getPSOMap()
{
	return psoProcessorMap;
}

odepso::MethodWrapper::RandomStruct& odepso::MethodWrapper::getRandStruct()
{
	return randomStruct;
}

const std::map<odepso::Common::SOLVER_TYPES, std::unique_ptr<odepso::SolverIF>, odepso::Common::CompSolvers>& odepso::MethodWrapper::getMethodsMap() const
{
	return methodsMap;
}

const std::map<odepso::Common::SOLVER_TYPES, odepso::OdeSolverParameters, odepso::Common::CompSolvers>& odepso::MethodWrapper::getParamsMap() const
{
	return paramsMap;
}

const std::map<odepso::Common::SOLVER_TYPES, std::vector<std::pair<odepso::OdeSolverParameters, Eigen::VectorXd>>, odepso::Common::CompSolvers>& odepso::MethodWrapper::getResultsMap() const
{
	return resultsMap;
}

const std::map<odepso::Common::SOLVER_TYPES, odepso::ParticleProcessor, odepso::Common::CompSolvers>& odepso::MethodWrapper::getPSOMap() const
{
	return psoProcessorMap;
}

const odepso::MethodWrapper::RandomStruct& odepso::MethodWrapper::getRandStruct() const
{
	return randomStruct;
}
