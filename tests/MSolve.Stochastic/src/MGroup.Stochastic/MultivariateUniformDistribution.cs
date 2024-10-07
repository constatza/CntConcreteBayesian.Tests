using System;
using Accord.Math.Random;
using Accord.Statistics.Distributions.Fitting;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions;
using Accord;
using Accord.Math;
using Accord.Statistics.Distributions.Sampling;
using System.Reflection.Emit;
using DotNumerics.Optimization;

namespace MGroup.Stochastic
{
	[Serializable]
	public class MultivariateUniformDistribution : MultivariateContinuousDistribution, IFittableDistribution<double[], IFittingOptions>, IFittable<double[], IFittingOptions>, IFittable<double[]>, IFittableDistribution<double[]>, IDistribution<double[]>, IDistribution, ICloneable, ISampleableDistribution<double[]>, IRandomNumberGenerator<double[]>
	{
		private double[] lowerBound;
		private double[] upperBound;
		private double[] mean;
		private double[] variance;
		private double[,] covariance;

		private double area;

		public double[] LowerBound { get => lowerBound; }

		public double[] UpperBound { get => upperBound; }


		public MultivariateUniformDistribution(double[] lowerBound, double[] upperBound) : base(lowerBound.Length)
		{
			this.lowerBound = lowerBound;
			this.upperBound = upperBound;
			this.mean = new double[lowerBound.Length];
			this.variance = new double[lowerBound.Length];
			this.covariance = new double[lowerBound.Length, lowerBound.Length];
			this.area = 1d;
			for (int i = 0; i < lowerBound.Length; i++)
			{
				this.mean[i] = (this.upperBound[i] + this.lowerBound[i]) / 2;
				this.variance[i] = Math.Pow(this.upperBound[i] - this.lowerBound[i], 2) / 12;
				this.covariance[i, i] = this.variance[i];
				this.area *= (this.upperBound[i] - this.lowerBound[i]);
			}
		}

		public override double[] Mean => mean;

		public override double[] Variance => variance;

		public override double[,] Covariance => covariance;

		public override object Clone() => new MultivariateUniformDistribution(lowerBound, upperBound);

		public override void Fit(double[][] observations, double[] weights, IFittingOptions options)
		{
			if (options != null)
			{
				throw new ArgumentException("This method does not accept fitting options.");
			}

			if (weights != null)
			{
				throw new ArgumentException("This distribution does not support weighted samples.");
			}
			double min = double.MaxValue;
			double max = double.MinValue;
			foreach (var observation in observations)
			{
				for (int i = 0; i < observation.Length; i++)
				{
					if (observation[i] < min)
					{
						lowerBound[i] = observation[i];
					}

					if (observation[i] > max)
					{
						upperBound[i] = observation[i];
					}
				}
			}
		}

		public override double ProbabilityDensityFunction(double[] x)
		{
			var probability = 1.0 / area;
			for (int i = 0; i < x.Length; i++)
			{
				if (x[i] > upperBound[i] || x[i] < lowerBound[i])
				{
					probability = 0;
					break;
				}
			}
			return probability;
		}

		public override double LogProbabilityDensityFunction(double[] x)
		{
			var logProbability = 0.0 - System.Math.Log(area);
			for (int i = 0; i < x.Length; i++)
			{
				if (x[i] > upperBound[i] || x[i] < lowerBound[i])
				{
					logProbability = 0;
					break;
				}
			}
			return logProbability;
		}

		public new double[][] Generate(int samples)
		{
			Random random = new Random();
			var totalSamples = new double[samples][];
			for(int i = 0; i < samples; i++)
			{
				totalSamples[i] = new double[Dimension];
				for (int j = 0; j < Dimension; j++)
				{
					totalSamples[i][j] = (upperBound[j] - lowerBound[j]) * random.NextDouble() + lowerBound[j];
				}
			}
			return totalSamples;
		}

		public new double[] Generate()
		{
			return Generate(1)[0];
		}

		public override string ToString(string format, IFormatProvider formatProvider)
		{
			return string.Format(formatProvider, "U(x; lowerBound, upperBound)");
		}
	}
}
