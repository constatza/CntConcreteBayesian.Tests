using System;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Math.Random;
using MGroup.LinearAlgebra.Matrices;
using System.Linq;

namespace MGroup.Stochastic
{
    public class HamiltonianMonteCarlo : IProbabilityDistributionSampler
    {
        Func<double[], double> model;
        Func<double[], int, double> gradModel;
        Func<double[]> initialModel;
        MultivariateNormalDistribution proposalDistribution;

        public double[] initialSample;
        private double[] currentSample;
        private double[] candidateSample;

        private int numModelParameters;
		private int totalDimensions;
        private int numSteps;
        private double stepSize;
        private Random randomSource;

        public HamiltonianMonteCarlo(int numModelParameters, Func<double[], double> model, Func<double[], int, double> gradModel, MultivariateNormalDistribution proposal, int numSteps = 25, double stepSize = 0.25)
        {
            this.model = model;
            this.gradModel = gradModel;
            this.proposalDistribution = proposal;
            this.numModelParameters = numModelParameters;
			this.totalDimensions = proposal.Dimension;
            this.numSteps = numSteps;
            this.stepSize = stepSize;
            this.currentSample = new double[totalDimensions];
            this.candidateSample = new double[totalDimensions];
            this.randomSource = Generator.Random;
        }

        public double[,] GenerateSamples(int numSamples)
        {
            var samples = new double[numSamples, totalDimensions];
            var acceptedSamples = 0;
            var current_q = new double[totalDimensions];
            var current_gradModel = 0d;
            var current_gradProposal = 0d;
            var current_U = -model(current_q.Take(numModelParameters).ToArray());
            var current_K = 0d;
            var proposed_U = 0d;
            var proposed_K = 0d;

            while (acceptedSamples < numSamples)
            {
                var q = (double[])current_q.Clone();
                var p = proposalDistribution.Generate();
				var current_p = new double[p.Length];
				Array.Copy(p, current_p, p.Length);
                for (int k = 0; k < totalDimensions; k++)
                {
                    current_gradModel = -gradModel(q.Take(numModelParameters).ToArray(), k);
                    p[k] = p[k] - stepSize * current_gradModel / 2;
                    for (int j = 0; j < numSteps; j++)
                    {
                        current_gradProposal = CalculateGradMinusLogProposal(p, k);
                        q[k] = q[k] + stepSize * current_gradProposal;
                        if (j + 1 != numSteps)
                        {
                            current_gradModel = -gradModel(q.Take(numModelParameters).ToArray(), k);
                            p[k] = p[k] - stepSize * current_gradModel;
                        }
                    }
                    current_gradModel = -gradModel(q.Take(numModelParameters).ToArray(), k);
                    p[k] = p[k] - stepSize * current_gradModel / 2;
                    p[k] = -p[k];
                }
                current_K = CalculateMinusLogProposal(current_p);
                proposed_U = -model(q.Take(numModelParameters).ToArray());
                proposed_K = CalculateMinusLogProposal(p);

                if (Math.Log(randomSource.NextDouble()) < Math.Exp (current_U - proposed_U + current_K - proposed_K))
                {
                    acceptedSamples++;
                    for (int i = 0; i < totalDimensions; i++)
                    {
                        samples[acceptedSamples - 1, i] = q[i];
                        current_q[i] = q[i];
                        current_U = -model(current_q.Take(numModelParameters).ToArray());
                    }
                }
            }
            return samples;
        }

        private double CalculateMinusLogProposal(double[] input)
        {
            var proposalEvaluation = -proposalDistribution.LogProbabilityDensityFunction(input);
            return proposalEvaluation;
        }
        
        private double CalculateGradMinusLogProposal(double[] input, int dimension)
        {
            var invCovarianceMatrix = Matrix.CreateFromArray(proposalDistribution.Covariance).Invert();
            var invCovariance = invCovarianceMatrix.CopyToArray2D();
            var gradTerm = 0d;
            for (int j = 0; j < input.Length; j++)
            {
                gradTerm += (input[dimension] - proposalDistribution.Mean[dimension]) * invCovariance[dimension, j];
            }
            var gradProposalEvaluation = gradTerm;

            return gradProposalEvaluation;
        }
    }
}
