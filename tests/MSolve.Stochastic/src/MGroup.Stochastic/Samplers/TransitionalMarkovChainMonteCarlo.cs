using System;
using System.Collections.Generic;
using System.Text;
using MGroup.Stochastic.Interfaces;
using Accord.Statistics;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions.Univariate;
using Accord.Math.Random;
using Accord.Math.Decompositions;
using Accord.Math;
using System.Linq;

namespace MGroup.Stochastic
{
    public class TransitionalMarkovChainMonteCarlo : IProbabilityDistributionSampler
    {
        Func<double[], double> model;
		Func<double[]> initialSampler;

        public double[] initialSample;
        private double[] currentSample;
        private double[] candidateSample;

        private int numModelParameters;
		private int totalDimensions;
        private double coefficientOfVariation;
        private double scalingFactor;
        private Random randomSource;

        public TransitionalMarkovChainMonteCarlo(int numModelParameters, Func<double[], double> model, Func<double[]> initialSampler, double coefficientOfVariation = 1, double scalingFactor = 0.2)
        {
            this.model = model;
            this.initialSampler = initialSampler;
			this.totalDimensions = initialSampler().Length;
			this.numModelParameters = numModelParameters;
            this.coefficientOfVariation = coefficientOfVariation;
            this.scalingFactor = scalingFactor;
            this.currentSample = new double[numModelParameters];
            this.candidateSample = new double[numModelParameters];
            this.randomSource = Generator.Random;
        }

        public double[,] GenerateSamples(int numSamples)
        {
            var samples = new double[numSamples, totalDimensions];
            var modelEvaluations = new double[numSamples];
			var likelihoodEvaluations = new double[numSamples];
			for (int i = 0; i < numSamples; i++)
            {
                var priorSample = initialSampler();
                for (int j = 0; j < totalDimensions; j++)
                {
                    samples[i, j] = priorSample[j];
                }
				likelihoodEvaluations[i] = model(priorSample);
				//likelihoodEvaluations[i] = likelihood.LogProbabilityDensityFunction(modelEvaluations[i]);
            }
            var weights = new double[numSamples];
            var p_current = 0d;
            var b = 0.2d;
            var S = 1d;
            while (p_current < 1)
            {
                var p_temp = 0d;
                var p_lower = p_current;
                var p_upper = 1d;
                while (p_upper - p_lower > 0.0001)
                {
                    p_temp = (p_lower + p_upper) / 2;
                    for (int i = 0; i < numSamples; i++)
                    {
                        weights[i] = Math.Exp((p_temp - p_current) * likelihoodEvaluations[i]);
                    }
                    var currentCoV = Measures.StandardDeviation(weights) / Measures.Mean(weights);
                    if (currentCoV > coefficientOfVariation)
                    {
                        p_upper = p_temp;
                    }
                    else
                    {
                        p_lower = p_temp;
                    }
                }
                var sumWeights = weights.Sum();
                p_current = p_temp;
                if (1 - p_current < 0.001)
                {
                    p_current = 1;
                }
                var normalizedWeights = new double[numSamples];
                for (int i = 0; i < numSamples; i++)
                {
                    normalizedWeights[i] = weights[i] / sumWeights;
                }
                var meanProposal = new double[totalDimensions];
                var covarProposal = new double[totalDimensions, totalDimensions];
                for (int j = 0; j < totalDimensions; j++)
                {
                    for (int i = 0; i < numSamples; i++)
                    {
                        meanProposal[j] += normalizedWeights[i] * samples[i,j];
                    }
                }
                for (int ii = 0; ii < totalDimensions; ii++)
                    for (int jj = 0; jj < totalDimensions; jj++)
                        for (int i = 0; i < numSamples; i++)
                        {
                            covarProposal[ii, jj] += b * b * normalizedWeights[i] * (samples[i, ii] - meanProposal[ii]) * (samples[i, jj] - meanProposal[jj]);
                        }
                int[] randomSampleInd = GeneralDiscreteDistribution.Random(normalizedWeights, numSamples);
                var currentSamples = samples.Copy();
                var currentmodelEvaluations = modelEvaluations.Copy();
				var currentLikehoodEvaluations = likelihoodEvaluations.Copy();
                for (int i = 0; i < numSamples; i++)
                {
                    double[][] candidateSample = MultivariateNormalDistribution.Generate(samples: 1, mean: currentSamples.GetRow<double>(randomSampleInd[i]), covariance: covarProposal);
                    var candidateLikelihoodEvaluation = model(candidateSample[0]);
					//var candidateLikelihoodEvaluation = likelihood.LogProbabilityDensityFunction(candidateModelEvaluation);
                    var likelihoodEvaluationRatio = Math.Exp(p_current * (candidateLikelihoodEvaluation - likelihoodEvaluations[i]));
                    if (likelihoodEvaluationRatio > Generator.Random.NextDouble())
                    {
                        for (int j = 0; j < totalDimensions; j++)
                        {
                            samples[i,j] = candidateSample[0][j];
                            currentSamples[randomSampleInd[i], j] = candidateSample[0][j];
                        }
						likelihoodEvaluations[i] = candidateLikelihoodEvaluation;
						currentLikehoodEvaluations[randomSampleInd[i]] = candidateLikelihoodEvaluation;
                    }
                    else
                    {
                        for (int j = 0; j < totalDimensions; j++)
                        {
                            samples[i, j] = currentSamples[randomSampleInd[i],j];
                        }
						likelihoodEvaluations[i] = currentLikehoodEvaluations[randomSampleInd[i]];
                    }
                }
            }
            return samples;
        }
    }
}
