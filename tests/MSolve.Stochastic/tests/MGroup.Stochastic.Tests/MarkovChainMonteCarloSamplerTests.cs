using System;
using System.Collections.Generic;
using System.Linq;

using Accord;
using Accord.Statistics;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Xunit;

namespace MGroup.Stochastic.Tests
{

	public static class MarkovChainMonteCarloSamplerTests
	{
		[Fact]
		public static void MetropolisHastingsTest()
		{
			var proposal = new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 1 } });
			Func<double[], double> prior = x => new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 1 } }).LogProbabilityDensityFunction(x[0]);
			Func<double[], double> likelihood = x => new MultivariateNormalDistribution(new double[] { 2 }, new double[,] { { 0.5 * 0.5 } }).LogProbabilityDensityFunction(x[0]);
			Func<double[], double> model = x => likelihood(x) + prior(x);
			var sampler = new MetropolisHastings(1, model, proposal);
			var samples = sampler.GenerateSamples(10000);
			var mean = samples.Mean();
			var std = samples.StandardDeviation();
			Assert.True(Math.Abs(mean - 1.6) < 0.05);
			Assert.True(Math.Abs(std[0] - 0.45) < 0.1);
		}

		[Fact]
		public static void TransitionalMarkovChainMonteCarloTest()
		{
			var prior = new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 1 } });
			Func<double[]> priorSampler = () => prior.Generate();
			Func<double[], double> model = x => new MultivariateNormalDistribution(new double[] { 2 }, new double[,] { { 0.5 * 0.5 } }).LogProbabilityDensityFunction(x[0]);
			var sampler = new TransitionalMarkovChainMonteCarlo(1, model, priorSampler, scalingFactor: 0.2);
			var samples = sampler.GenerateSamples(10000);
			var mean = samples.Mean();
			var std = samples.StandardDeviation();
			Assert.True(Math.Abs(mean - 1.6) < 0.05);
			Assert.True(Math.Abs(std[0] - 0.45) < 0.1);
		}

		[Fact]
		public static void HamiltonianMonteCarloTest()
		{
			var proposal = new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 1 } });
			Func<double[], double> prior = x => new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 1 } }).LogProbabilityDensityFunction(x[0]);
			Func<double[], double> likelihood = x => new MultivariateNormalDistribution(new double[] { 2 }, new double[,] { { 0.5 * 0.5 } }).LogProbabilityDensityFunction(x[0]);
			Func<double[], double> model = x => likelihood(x) + prior(x);
			Func<double[], double> gradPrior = x => -x[0];
			Func<double[], double> gradLikelihood = x => -2 * (2 * x[0] - 4);
			Func<double[], int, double> gradModel = (x,y) => gradPrior(x) + gradLikelihood(x);
			var sampler = new HamiltonianMonteCarlo(1, model, gradModel, proposal);
			var samples = sampler.GenerateSamples(10000);
			var mean = samples.Mean();
			var std = samples.StandardDeviation();
			Assert.True(Math.Abs(mean - 1.6) < 0.05);
			Assert.True(Math.Abs(std[0] - 0.45) < 0.1);
		}

		[Fact]
		public static void HierarchicalTransitionalMarkovChainMonteCarloTest()
		{
			var mean_hyperprior = new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 0.1 } });
			var var_hyperprior = new MultivariateNormalDistribution(new double[] { 1 }, new double[,] { { 0.01 } });
			//var lower_bound_hyperprior = new MultivariateUniformDistribution(new double[] { 0 }, new double[] { 1 });
			//var upper_bound_hyperprior = new MultivariateUniformDistribution(new double[] { 1 }, new double[] { 2 });
			var hyperpriorList = new Dictionary<string, ISampleableDistribution<double[]>>();
			hyperpriorList.Add("mean", mean_hyperprior);
			hyperpriorList.Add("var", var_hyperprior);
			var prior = new MultivariateNormalDistribution(new double[] { 0 }, new double[,] { { 0.1 } });
			Func<double[]> priorSampler = () => HyperpriorGenerate(prior, hyperpriorList);
			Func<double[], double> model = x => new MultivariateNormalDistribution(new double[] { 5 }, new double[,] { { 0.5 * 0.5 } }).LogProbabilityDensityFunction(x[0]);
			var sampler = new TransitionalMarkovChainMonteCarlo(1, model, priorSampler, scalingFactor: 0.2);
			var samples = sampler.GenerateSamples(10000);
			var mean = samples.Mean(0);
			var std = samples.StandardDeviation();
			//Assert.True(Math.Abs(mean - 1.6) < 0.05);
			//Assert.True(Math.Abs(std[0] - 0.45) < 0.1);

			static double[] HyperpriorGenerate(ISampleableDistribution<double[]> prior, Dictionary<string, ISampleableDistribution<double[]>> hyperpriorList)
			{
				var length = prior.Generate().Length;
				foreach (var hyperpriorDistribution in hyperpriorList.Values)
				{
					length += hyperpriorDistribution.Generate().Length;
				}
				var totalPriorSample = new double[length];
				if(prior is MultivariateNormalDistribution)
				{
					var prior_temp = (MultivariateNormalDistribution)prior;
					if (hyperpriorList.ContainsKey("mean") == true && hyperpriorList.ContainsKey("var") == true)
					{
						var mean_hyperprior_sample = hyperpriorList["mean"].Generate();
						var var_hyperprior_sample = hyperpriorList["var"].Generate();
						var priorSampleDistribution = new MultivariateNormalDistribution(mean_hyperprior_sample, CreateDiagonalArrayFromVector(var_hyperprior_sample));
						var priorSample = priorSampleDistribution.Generate();
						totalPriorSample = priorSample.Concat(mean_hyperprior_sample).Concat(var_hyperprior_sample).ToArray();
					}
					if (hyperpriorList.ContainsKey("mean") == true && hyperpriorList.ContainsKey("var") == false)
					{
						var mean_hyperprior_sample = hyperpriorList["mean"].Generate();
						var priorSampleDistribution = new MultivariateNormalDistribution(hyperpriorList["mean"].Generate(), prior_temp.Covariance);
						var priorSample = priorSampleDistribution.Generate();
						totalPriorSample = (double[])priorSample.Concat(mean_hyperprior_sample).ToArray();
					}
					if (hyperpriorList.ContainsKey("mean") == false && hyperpriorList.ContainsKey("var") == true)
					{
						var var_hyperprior_sample = hyperpriorList["var"].Generate();
						var priorSampleDistribution = new MultivariateNormalDistribution(prior_temp.Mean, CreateDiagonalArrayFromVector(var_hyperprior_sample));
						var priorSample = priorSampleDistribution.Generate();
						totalPriorSample = (double[])priorSample.Concat(var_hyperprior_sample).ToArray();
					}

				}
				//var prior = new MultivariateNormalDistribution(mean_hyperprior_sample , new double[][] { var_hyperprior_sample } );
				double[,] CreateDiagonalArrayFromVector(double[] vector)
				{
					var diagonalArray = new double[vector.Length, vector.Length];
					for (int i = 0; i < vector.Length; i++)
					{
						diagonalArray[i, i] = vector[i];
					}
					return diagonalArray;
				}
				return totalPriorSample;
			}
		}
	}
}
