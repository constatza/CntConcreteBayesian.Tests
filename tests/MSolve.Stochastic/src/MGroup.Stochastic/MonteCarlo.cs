using MGroup.Stochastic.Interfaces;
using System;
using Accord.Math.Random;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using DotNumerics.Optimization;

namespace MGroup.Stochastic
{
    public class MonteCarlo
    {
		public int NumModelParameters { get; }
		public Func<double[], double[]> Model { get; }
		public ISampleableDistribution<double[]> Distribution { get; }
		public double[,] MonteCarloRealizationOutputs;
		public double[] MonteCarloMeanValue;
        public double[] MonteCarloStandardDeviation;

		/// <summary>Initializes a new instance of the <see cref="MonteCarlo"/> class.</summary>
		/// <param name="numberOfSamples">The no of samples.</param>
		/// <param name="model">The model.</param>
		/// <param name="distribution">The sampling distribution.</param>
		public MonteCarlo(int numModelParameters, Func<double[], double[]> model, ISampleableDistribution<double[]> distribution)
        {
			this.NumModelParameters = numModelParameters;
			this.Model = model;
			this.Distribution = distribution;
        }

        /// <summary>Evaluates 1st and 2nd statistical moments of the predesignated response.</summary>
        public double[,] GenerateSamples(int numSamples)
        {
			var sample = Distribution.Generate();
			var modelValue = Model(sample);
			MonteCarloRealizationOutputs = new double[numSamples, modelValue.Length];
			MonteCarloMeanValue = new double[modelValue.Length];
			MonteCarloStandardDeviation = new double[modelValue.Length];
			for (int j = 0; j < modelValue.Length; j++)
			{
				MonteCarloRealizationOutputs[0, j] = modelValue[j];
				MonteCarloMeanValue[j] += modelValue[j];
			}
			for (int i = 1; i < numSamples; i++)
			{
				sample = Distribution.Generate();
				modelValue = Model(sample);
				for (int j = 0; j < modelValue.Length; j++)
				{
					MonteCarloRealizationOutputs[i, j] = modelValue[j];
					MonteCarloMeanValue[j] += modelValue[j];
				}
			}
			for (int i = 0; i < modelValue.Length; i++)
			{
				MonteCarloMeanValue[i] /= numSamples;
				for (int j = 0; j < numSamples; j++)
				{
					MonteCarloStandardDeviation[i] += Math.Pow(MonteCarloRealizationOutputs[j, i] - MonteCarloMeanValue[i], 2);
				}
				MonteCarloStandardDeviation[i] = Math.Sqrt(MonteCarloStandardDeviation[i] / (numSamples - 1));
			}
			return MonteCarloRealizationOutputs;
		}
    }
}
