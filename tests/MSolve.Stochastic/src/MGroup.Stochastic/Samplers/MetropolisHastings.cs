using Accord.Math.Random;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using MGroup.Stochastic.Interfaces;
using System;
using System.Linq;

namespace MGroup.Stochastic
{
	public class MetropolisHastings : IProbabilityDistributionSampler
	{
		Func<double[], double> model;
		MultivariateNormalDistribution proposalDistribution;

		public double[] initialSample;
		private double[] currentSample;
		private double[] candidateSample;

		private int burnIn;
		private int rejectionRate;
		private int numModelParameters;
		private Random randomSource;
		private int totalDimensions;

		public MetropolisHastings(int numModelParameters, Func<double[], double> model, MultivariateNormalDistribution proposal, int burnIn = 0, int rejectionRate = 0)
		{
			this.model = model;
			this.proposalDistribution = proposal;
			this.burnIn = burnIn;
			this.rejectionRate = rejectionRate;
			this.numModelParameters = numModelParameters;
			this.totalDimensions = proposal.Dimension;
			this.currentSample = new double[totalDimensions];
			this.candidateSample = new double[totalDimensions];
			this.randomSource = Generator.Random;
		}

		public double[,] GenerateSamples(int numSamples)
		{
			var samples = new double[numSamples, totalDimensions];
			var acceptedSamples = 0;
			var currentEvaluation = 0d;
			if (initialSample == null)
				currentEvaluation = model(proposalDistribution.Mean.Take(numModelParameters).ToArray());
			else
				currentEvaluation = model(initialSample.Take(numModelParameters).ToArray());
			while (acceptedSamples < numSamples + burnIn)
			{
				candidateSample = proposalDistribution.Generate();
				var candidateEvaluation = model(candidateSample.Take(numModelParameters).ToArray());
				double Ratio = candidateEvaluation - currentEvaluation;
				if (Math.Log(randomSource.NextDouble()) < Ratio)
				{
					acceptedSamples++;
					if (acceptedSamples > burnIn)
					{
						for (int i = 0; i < totalDimensions; i++)
						{
							samples[acceptedSamples - burnIn - 1, i] = candidateSample[i];
							currentSample = candidateSample;
							currentEvaluation = candidateEvaluation;
						}
					}
					proposalDistribution = new MultivariateNormalDistribution(candidateSample, proposalDistribution.Covariance);
				}
			}
			return samples;
		}
	}
}
