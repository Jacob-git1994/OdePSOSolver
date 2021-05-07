#pragma once

namespace odepso
{
	class ParticleParameters
	{
	private:

		//Seed
		long long seed;

		//Number of particles allows
		size_t numOfParticles;

		//Control the velocity
		double omega;

		//Phi1
		double phi1;

		//Phi2
		double phi2;

		//Learning rate
		double learningRate;

		//PSO Convergence Criteria
		size_t maxIterations;
		double psoTolerance;

	public:

		ParticleParameters();

		//Constructor
		ParticleParameters(
			const long long& seedIn,
			const size_t& numOfParticlesIn,
			const double& omegaIn,
			const double& phi1In,
			const double& phi2In,
			const double& learningRateIn,
			const size_t& maxIterIn,
			const double& psoTolIn);

		//Default Copy Constructor
		ParticleParameters(const ParticleParameters& particleParametersIn) = default;

		//Default Assign copy constructor
		ParticleParameters& operator=(const ParticleParameters& particleParametersIn) = default;

		//Destructor
		virtual ~ParticleParameters() = default;

		//Get the seed
        const long long& getSeed() const;

		//Set the seed
        void setSeed(const long long& seedIn);

		//Get the number of allowed particles
        const size_t& getNumOfParticles() const;

		//Set the number of allowed particles
        void setNumOfParticles(const size_t& numOfParticlesIn);

		//Get omega
        const double& getOmega() const;

		//Set omega
        void setOmega(const double& omegaIn);

		//Get phi one
        const double& getPhi1() const;

		//Set phi one
        void setPhi1(const double& phi1In);

		//Get phi two
        const double& getPhi2() const;

		//Set phi two
        void setPhi2(const double& phi2In);

		//Get the learning rate
        const double& getLearningRate() const;

		//Set the learning rate
        void setLearningRate(const double& learningRateIn);

		//Get the max number of iterations
		const size_t& getMaxPSOIterations() const;

		//Set the max number of iterations
		void setMaxPSOITerations(const size_t& maxIterationsIn);

		//Get the pso tolerance
		const double& getPSOTol() const;

		//Set the pso tolerance
		void setPSOTol(const double& psoTolIn);
	};
}