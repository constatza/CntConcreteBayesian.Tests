using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

using Accord;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;

namespace MGroup.Stochastic.Bayesian
{
	public class HierarchicalBayesianUpdate
	{
		private ISampleableDistribution<double[]> priorParameterDistribution;
		private ISampleableDistribution<double[]>[] priorHyperparameterDistribution;
		private string[] priorHyperparameterType;
		private Dictionary<string, ISampleableDistribution<double[]>> priorHyperparameterList;
		private MultivariateContinuousDistribution likelihoodFunction;
		private Func<double[], double[]> model;

		double[] measurementValues;
		double[] measurementError;
		int numParameters;
		int[] numHyperparameters;

		public HierarchicalBayesianUpdate(Func<double[], double[]> model, ISampleableDistribution<double[]> priorParameterDistribution, Dictionary<string, ISampleableDistribution<double[]>> priorHyperparameterList, double[] measurementValues, double[] measurementError)
		{
			this.priorParameterDistribution = priorParameterDistribution;
			numParameters = priorParameterDistribution.Generate().Length;
			var count = 0;
			numHyperparameters = new int[priorHyperparameterList.Count()];
			priorHyperparameterDistribution = new ISampleableDistribution<double[]>[priorHyperparameterList.Count()];
			foreach (var hyperpriorDistribution in priorHyperparameterList.Values)
			{
				priorHyperparameterDistribution[count] = hyperpriorDistribution;
				numHyperparameters[count] = priorHyperparameterDistribution[count].Generate().Length;
				count++;
			}
			count = 0;
			priorHyperparameterType = new string[priorHyperparameterList.Count()];
			foreach (var hyperpriorType in priorHyperparameterList.Keys)
			{
				priorHyperparameterType[count] = hyperpriorType;
				count++;
			}
			this.priorHyperparameterList = priorHyperparameterList;
			this.measurementValues = measurementValues;
			this.measurementError = measurementError;
			this.model = model;
			likelihoodFunction = CreateLikelihoodFunction();
		}

		public double PosteriorModel(double[] sample)
		{
			var posteriorModelValue = likelihoodFunction.LogProbabilityDensityFunction(model(sample.Take(numParameters).ToArray())) + priorParameterDistribution.LogProbabilityFunction(sample.Take(numParameters).ToArray());
			var prevNumHyperParameters = 0;
			for(int i = 0; i < priorHyperparameterDistribution.Length; i++) 
			{
				posteriorModelValue += priorHyperparameterDistribution[i].LogProbabilityFunction(sample.Skip(numParameters + prevNumHyperParameters).Take(numHyperparameters[i]).ToArray());
				prevNumHyperParameters += numHyperparameters[i];
			}
			return posteriorModelValue;
		}

		public double LikelihoodModel(double[] sample) => likelihoodFunction.LogProbabilityDensityFunction(model(sample.Take(numParameters).ToArray()));

		public ISampleableDistribution<double[]>[] PriorHyperparameterDistribution { get => priorHyperparameterDistribution; }

		public ISampleableDistribution<double[]> PriorDistribution { get => priorParameterDistribution; }

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

		public double[] PriorDistributionGenerator()
		{
			var length = priorParameterDistribution.Generate().Length;
			foreach (var hyperpriorDistribution in priorHyperparameterList.Values)
			{
				length += hyperpriorDistribution.Generate().Length;
			}
			var totalPriorSample = new double[length];
			if (priorParameterDistribution is MultivariateNormalDistribution)
			{
				var prior_temp = (MultivariateNormalDistribution)priorParameterDistribution;
				if (priorHyperparameterList.ContainsKey("mean") == true && priorHyperparameterList.ContainsKey("var") == true)
				{
					var mean_hyperprior_sample = priorHyperparameterList["mean"].Generate();
					var var_hyperprior_sample = priorHyperparameterList["var"].Generate();
					var priorSampleDistribution = new MultivariateNormalDistribution(mean_hyperprior_sample, CreateDiagonalArrayFromVector(var_hyperprior_sample));
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(mean_hyperprior_sample).Concat(var_hyperprior_sample).ToArray();
				}
				if (priorHyperparameterList.ContainsKey("mean") == true && priorHyperparameterList.ContainsKey("var") == false)
				{
					var mean_hyperprior_sample = priorHyperparameterList["mean"].Generate();
					var priorSampleDistribution = new MultivariateNormalDistribution(mean_hyperprior_sample, prior_temp.Covariance);
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(mean_hyperprior_sample).ToArray();
				}
				if (priorHyperparameterList.ContainsKey("mean") == false && priorHyperparameterList.ContainsKey("var") == true)
				{
					var var_hyperprior_sample = priorHyperparameterList["var"].Generate();
					var priorSampleDistribution = new MultivariateNormalDistribution(prior_temp.Mean, CreateDiagonalArrayFromVector(var_hyperprior_sample));
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(var_hyperprior_sample).ToArray();
				}
			}
			else if (priorParameterDistribution is MultivariateUniformDistribution)
			{
				var prior_temp = (MultivariateUniformDistribution)priorParameterDistribution;
				if (priorHyperparameterList.ContainsKey("lower") == true && priorHyperparameterList.ContainsKey("upper") == true)
				{
					var lower_hyperprior_sample = priorHyperparameterList["lower"].Generate();
					var upper_hyperprior_sample = priorHyperparameterList["upper"].Generate();
					var priorSampleDistribution = new MultivariateUniformDistribution(lower_hyperprior_sample, upper_hyperprior_sample);
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(lower_hyperprior_sample).Concat(upper_hyperprior_sample).ToArray();
				}
				if (priorHyperparameterList.ContainsKey("lower") == true && priorHyperparameterList.ContainsKey("upper") == false)
				{
					var lower_hyperprior_sample = priorHyperparameterList["lower"].Generate();
					var priorSampleDistribution = new MultivariateUniformDistribution(lower_hyperprior_sample, prior_temp.UpperBound);
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(lower_hyperprior_sample).ToArray();
				}
				if (priorHyperparameterList.ContainsKey("lower") == false && priorHyperparameterList.ContainsKey("upper") == true)
				{
					var upper_hyperprior_sample = priorHyperparameterList["upper"].Generate();
					var priorSampleDistribution = new MultivariateUniformDistribution(prior_temp.LowerBound, upper_hyperprior_sample);
					var priorSample = priorSampleDistribution.Generate();
					totalPriorSample = priorSample.Concat(upper_hyperprior_sample).ToArray();
				}
			}
			return totalPriorSample;
		}

		private double[,] CreateDiagonalArrayFromVector(double[] vector)
		{
			var diagonalArray = new double[vector.Length, vector.Length];
			for (int i = 0; i < vector.Length; i++)
			{
				diagonalArray[i, i] = vector[i];
			}
			return diagonalArray;
		}

	}
}
