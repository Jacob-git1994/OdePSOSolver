#include "ParticleParameters.h"

odepso::ParticleParameters::ParticleParameters() :
    seed(0),
    numOfParticles(0),
    omega(0.0),
    phi1(0.0),
    phi2(0.0),
    learningRate(0.0),
    maxIterations(0),
    psoTolerance(0.0)
{
    //Nothing else to do here
}

odepso::ParticleParameters::ParticleParameters(
    const long long& seedIn, 
    const size_t& numOfParticlesIn,
    const double& omegaIn, 
    const double& phi1In, 
    const double& phi2In, 
    const double& learningRateIn, 
    const size_t& maxIterIn,
    const double& psoTolIn) : 
    seed(seedIn),
    numOfParticles(numOfParticlesIn),
    omega(omegaIn),
    phi1(phi1In),
    phi2(phi2In),
    learningRate(learningRateIn),
    maxIterations(maxIterIn),
    psoTolerance(psoTolIn)
{
    //Nothing else to do here
}

const long long& odepso::ParticleParameters::getSeed() const
{
    return seed;
}

void odepso::ParticleParameters::setSeed(const long long& seedIn)
{
    this->seed = seedIn;
}

const size_t& odepso::ParticleParameters::getNumOfParticles() const
{
    return numOfParticles;
}

void odepso::ParticleParameters::setNumOfParticles(const size_t& numOfParticlesIn)
{
    this->numOfParticles = numOfParticlesIn;
}

const double& odepso::ParticleParameters::getOmega() const
{
    return omega;
}

void odepso::ParticleParameters::setOmega(const double& omegaIn)
{
    this->omega = omegaIn;
}

const double& odepso::ParticleParameters::getPhi1() const
{
    return phi1;
}

void odepso::ParticleParameters::setPhi1(const double& phi1In)
{
    this->phi1 = phi1In;
}

const double& odepso::ParticleParameters::getPhi2() const
{
    return phi2;
}

void odepso::ParticleParameters::setPhi2(const double& phi2In)
{
    this->phi2 = phi2In;
}

const double& odepso::ParticleParameters::getLearningRate() const
{
    return learningRate;
}

void odepso::ParticleParameters::setLearningRate(const double& learningRateIn)
{
    this->learningRate = learningRateIn;
}

const size_t& odepso::ParticleParameters::getMaxPSOIterations() const
{
    return maxIterations;
}

void odepso::ParticleParameters::setMaxPSOITerations(const size_t& maxIterationsIn)
{
    maxIterations = maxIterationsIn;
}

const double& odepso::ParticleParameters::getPSOTol() const
{
    return psoTolerance;
}

void odepso::ParticleParameters::setPSOTol(const double& psoTolIn)
{
    psoTolerance = psoTolIn;
}


