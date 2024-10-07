using System;
using MGroup.Stochastic.Interfaces;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Distributions.Multivariate;

namespace MGroup.Stochastic.Structural.StochasticRealizers
{
    public class RandomVariable : IProbabilityDistributionSampler
    {
        private readonly ISampleableDistribution<double[]> _distribution;

        public RandomVariable(ISampleableDistribution<double[]> distribution)
        {
            _distribution = distribution;
        }

        public double[,] GenerateSamples(int numSamples)
        {
			var samples = new double[numSamples, _distribution.Generate().Length];
			for (int i = 0; i < numSamples; i++)
			{
				var currentSample = _distribution.Generate();
				for (int j = 0; j < currentSample.Length; j++)
				{
					samples[i, j] = currentSample[j];
				}
			}
			return(samples);
        }

    }
}
