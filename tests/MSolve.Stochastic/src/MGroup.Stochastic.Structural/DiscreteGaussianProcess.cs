using System;
using Accord.Statistics.Distributions;
using Accord.Statistics.Distributions.Multivariate;

namespace MGroup.Stochastic.Structural
{

	public class DiscreteGaussianProcess
	{
		private double[] steps;
		private double std;
		private double corLength;
		private bool zeroInitialValue;
		public Func<double, double> drift;

		public DiscreteGaussianProcess(double[] steps, double std, double corLength, bool zeroInitialValue = false, Func<double, double> drift = null)
		{
			this.steps = steps;
			this.std = std;
			this.corLength = corLength;
			this.zeroInitialValue = zeroInitialValue;
			this.drift = drift;
		}

		public double[] Generate()
		{
			var mean = new double[steps.Length];
			var covariance = new double[steps.Length, steps.Length];
			var response = new double[steps.Length];
			for (int i = 0; i < steps.Length; i++)
			{
				if (drift != null)
				{
					mean[i] = drift(steps[i]);
				}
				else
				{
					mean[i] = 0;
				}

				for (int j = 0; j < steps.Length; j++)
				{
					covariance[i, j] = Math.Pow(std, 2) * Math.Exp(-Math.Pow(steps[i] - steps[j], 2) / (2 * Math.Pow(corLength, 2)));
				}
			}
			//var mvnDistr = new MultivariateNormalDistribution(mean, covariance);
			//response = mvnDistr.Generate();
			response = MultivariateNormalDistribution.Generate(1, mean, covariance)[0];
			var init_response = response[0];
			if (zeroInitialValue == true)
			{
				for (int j = 0; j < steps.Length; j++)
				{
					response[j] = response[j] - init_response;
				}
			}
			return response;
		}
	}
}
