using System;
using Accord.Statistics;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using MGroup.Stochastic.Bayesian;
using MGroup.Stochastic.Tests.SupportiveClasses;

using Xunit;

namespace MGroup.Stochastic.Tests
{

	public static class BayesianUpdateTests
	{
		[Fact]
		public static void FiniteElementModelTest()
		{
			var proposal = new MultivariateNormalDistribution(new double[] { 1353000 }, new double[,] { { 100000d * 100000d } });
			var prior = new MultivariateUniformDistribution(new double[] { 353000 }, new double[] { 2353000 });
			var measurementValues = new double[] { -0.1541071095408108 };
			var measurementError = new double[] { 0.01 * 0.01 };
			var model = new FiniteElementModel();
			var bayesianInstance = new BayesianUpdate(model.CreateModel, prior, measurementValues, measurementError);
			var sampler = new MetropolisHastings(1, bayesianInstance.PosteriorModel, proposal);
			var samples = sampler.GenerateSamples(1000);
			var mean = samples.Mean();
			var std = samples.StandardDeviation();
			Assert.True(Math.Abs(mean - 1353000) < 50000);
			Assert.True(Math.Abs(std[0] - 100000d) < 20000);
		}

		[Fact]
		public static void FiniteElementModelTestTMCMC()
		{
			//var proposal = new MultivariateNormalDistribution(new double[] { 1353000 }, new double[,] { { 100000d * 100000d } });
			var prior = new MultivariateUniformDistribution(new double[] { 353000 }, new double[] { 2353000 });
			Func<double[]> priorSampler = () => prior.Generate();
			var measurementValues = new double[] { -0.1541071095408108 };
			var measurementError = new double[] { 0.01 * 0.01 };
			var model = new FiniteElementModel();
			var bayesianInstance = new BayesianUpdate(model.CreateModel, prior, measurementValues, measurementError);
			var sampler = new TransitionalMarkovChainMonteCarlo(1, bayesianInstance.LikelihoodModel, priorSampler, scalingFactor: 0.2);
			var samples = sampler.GenerateSamples(10000);
			var mean = samples.Mean();
			var std = samples.StandardDeviation();
			Assert.True(Math.Abs(mean - 1353000) < 50000);
			Assert.True(Math.Abs(std[0] - 100000d) < 20000);
		}

		//[Fact]
		//public static void FiniteElementModelTestHMC()
		//{
		//	var proposal = new MultivariateNormalDistribution(new double[] { 1353000 }, new double[,] { { 100000d * 100000d } });
		//	var prior = new MultivariateUniformDistribution(new double[] { 353000 }, new double[] { 2353000 });
		//	var measurementValues = new double[] { -0.1541071095408108 };
		//	var measurementError = new double[] { 0.01 * 0.01 };
		//	var model = new FiniteElementModel();
		//	var bayesianInstance = new BayesianUpdate(model.CreateModel, prior, measurementValues, measurementError);
		//	var sampler = new MetropolisHastings(1, bayesianInstance.PosteriorModel, proposal);
		//	var samples = sampler.GenerateSamples(1000);
		//	var mean = samples.Mean();
		//	var std = samples.StandardDeviation();
		//	Assert.True(Math.Abs(mean - 1353000) < 50000);
		//	Assert.True(Math.Abs(std[0] - 100000d) < 20000);
		//}
	}
}
