using System;

using Accord;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;

namespace MGroup.Stochastic.Bayesian
{
	public class BayesianUpdate
	{
		private ISampleableDistribution<double[]> priorDistribution;
		private MultivariateContinuousDistribution likelihoodFunction;
		private Func<double[], double[]> model;

		double[] measurementValues;
		double[] measurementError;

		public BayesianUpdate(Func<double[], double[]> model, ISampleableDistribution<double[]> priorDistribution, double[] measurementValues, double[] measurementError)
		{
			this.priorDistribution = priorDistribution;
			this.measurementValues = measurementValues;
			this.measurementError = measurementError;
			this.model = model;
			likelihoodFunction = CreateLikelihoodFunction();
		}

		public double PosteriorModel(double[] sample) => likelihoodFunction.LogProbabilityDensityFunction(model(sample)) + priorDistribution.LogProbabilityFunction(sample);

		public double LikelihoodModel(double[] sample) => likelihoodFunction.LogProbabilityDensityFunction(model(sample));

		public ISampleableDistribution<double[]> PriorDistribution { get => priorDistribution; }

		public double[] PriorDistributionSampler { get => priorDistribution.Generate(); }

		public ISampleableDistribution<double[]> LikelihoodFuction { get => likelihoodFunction; }

		private MultivariateContinuousDistribution CreateLikelihoodFunction()
		{
			var measurementErrors = new double[measurementValues.Length, measurementValues.Length];
			for (int i = 0; i < measurementErrors.GetLength(0); i++)
			{
				measurementErrors[i, i] = measurementError[i];
			}

			return new MultivariateNormalDistribution(measurementValues, measurementErrors);
		}

		public double PosteriorFunctionEvaluator(double[] input)
		{
			var modelEvaluation = model(input);
			var likelihoodEvaluation = likelihoodFunction.ProbabilityDensityFunction(modelEvaluation);
			var priorEvaluation = priorDistribution.LogProbabilityFunction(input);
			return likelihoodEvaluation * priorEvaluation;
		}

		public double LogPosteriorFunctionEvaluator(double[] input)
		{
			var modelEvaluation = model(input);
			var likelihoodEvaluation = likelihoodFunction.LogProbabilityDensityFunction(modelEvaluation);
			var priorEvaluation = priorDistribution.LogProbabilityFunction(input);
			return likelihoodEvaluation + priorEvaluation;
		}

		public double LikelihoodFunctionEvaluator(double[] input)
		{
			var modelEvaluation = model(input);
			var likelihoodEvaluation = likelihoodFunction.ProbabilityDensityFunction(modelEvaluation);
			return likelihoodEvaluation;
		}

		public double LogLikelihoodFunctionEvaluator(double[] input)
		{
			var modelEvaluation = model(input);
			var likelihoodEvaluation = likelihoodFunction.LogProbabilityDensityFunction(modelEvaluation);
			return likelihoodEvaluation;
		}

		public double GradLogModelEvaluator(double[] input, int dimension)
		{
			var increment = 1e-8;
			var modelEvaluation = model(input);
			var inputNext = new double[input.Length];
			Array.Copy(input, inputNext, input.Length);
			inputNext[dimension] += increment;
			var nextModelEvaluation = model(inputNext);
			var priorEvaluation = priorDistribution.LogProbabilityFunction(input);
			var nextPriorEvaluation = priorDistribution.LogProbabilityFunction(inputNext);
			var likelihoodEvaluation = LogLikelihoodFunctionEvaluator(modelEvaluation);
			var nextLikelihoodEvaluation = LogLikelihoodFunctionEvaluator(nextModelEvaluation);
			var gradLikelihood = (nextLikelihoodEvaluation - likelihoodEvaluation) / increment;
			var gradPrior = (nextPriorEvaluation - priorEvaluation);
			var gradLogModel = priorEvaluation * gradLikelihood + gradPrior * likelihoodEvaluation;
			//var invLikelihoodCovarianceMatrix = Matrix.CreateFromArray(likelihoodFunction.Covariance).Invert();
			//var invLikelihoodCovariance = invLikelihoodCovarianceMatrix.CopyToArray2D();
			//var invPriorCovarianceMatrix = Matrix.CreateFromArray(priorDistribution.Covariance).Invert();
			//var invPriorCovariance = invPriorCovarianceMatrix.CopyToArray2D();
			//var gradLikelihoodTerm = 0d;
			//var gradPriorTerm = 0d;
			//for (int j = 0; j < modelEvaluation.Length; j++)
			//         {
			//             gradLikelihoodTerm += (modelEvaluation[dimension] - likelihoodFunction.Mean[dimension]) * invLikelihoodCovariance[dimension, j];
			//             gradPriorTerm += (input[dimension] - priorDistribution.Mean[dimension]) * invPriorCovariance[dimension, j];
			//         }
			//         var gradLogLikelihoodEvaluation = gradLikelihoodTerm + gradPriorTerm; // -Math.Log(gradModel * gradTerm * proposalEvaluation * priorEvaluation)

			return gradLogModel;
		}

		public double[] PriorDistributionGenerator()
		{
			var priorSample = priorDistribution.Generate();
			return priorSample;
		}

		//public double[,] Solve(int numSamples)
		//{
		//	if (Sampler == null)
		//		throw new ArgumentException("Need to assign a sampler first");
		//	var samples = Sampler.GenerateSamples(numSamples);
		//	return samples;
		//}
	}
}
