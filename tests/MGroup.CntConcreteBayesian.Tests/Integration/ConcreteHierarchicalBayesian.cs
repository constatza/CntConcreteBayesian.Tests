using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Statistics.Distributions.Multivariate;
using Accord.Statistics.Distributions;
using Accord.Statistics;
using Xunit;
using MGroup.Stochastic;
using MGroup.CntConcreteBayesian.Tests.Commons;
using MGroup.CntConcreteBayesian.Tests.Models;
using MGroup.Stochastic.Bayesian;
using Accord.Statistics.Distributions.Univariate;
using System.IO;

namespace MGroup.CntConcreteBayesian.Tests.Integration
{
	public class ConcreteHierarchicalBayesian
	{
		[Fact]
		public static void RunSimulation()
		{

			var cementSampler = InitializeSimulation();

			var numParameterSamples = 5;

			var cementSamples = GenerateCementSamples(numParameterSamples, cementSampler);
			//var mortarSamples = GenerateMortarSamples(numParameterSamples, mortarSampler);
			//var concreteSamples = GenerateConcreteSamples(numParameterSamples, concreteSampler);

			//var cementSamples = LoadCementSamples();
			//var mortarSamples = LoadMortarSamples();
			//var concreteSamples = LoadConcreteSamples();

			var numHyperparameterSamples = 50000;
			var numNewParameterSamples = 50000;

			var hyperparameterPrior = new MultivariateUniformDistribution(new double[] { 0, 0, 0, 0, 0, 0 }, new double[] { 20, 10, 2, 1, 0.2, 0.1 });
			var hyperparameterSampler = new TransitionalMarkovChainMonteCarlo(6, HyperparameterLikelihoodModel, hyperparameterPrior.Generate, scalingFactor: 0.2);

			(var hyperparameterSamples, var newParameterSamples) = GenerateHyperparameterAndNewParameterSamples(numHyperparameterSamples, numNewParameterSamples, hyperparameterSampler);


			double HyperparameterLikelihoodModel(double[] sample)
			{
				var cementLikelihoodValue = 0d;
				var tempDistribution = new MultivariateUniformDistribution(new double[] { sample[0], sample[2], sample[4] }, new double[] { sample[0] + sample[1], sample[2] + sample[3], sample[4] + sample[5] });
				for (int i = 0; i < cementSamples.GetLength(0); i++)
				{
					var sampleTemp = new double[cementSamples.GetLength(1)];
					for (int j = 0; j < cementSamples.GetLength(1); j++)
					{
						sampleTemp[j] = cementSamples[i, j];
					}
					cementLikelihoodValue += tempDistribution.ProbabilityDensityFunction(sampleTemp);
				}
				cementLikelihoodValue = cementLikelihoodValue / cementSamples.GetLength(0);
				//var mortarLikelihoodValue = 0d;
				//for (int i = 0; i < mortarSamples.GetLength(0); i++)
				//{
				//	var sampleTemp = new double[mortarSamples.GetLength(1)];
				//	for (int j = 0; j < mortarSamples.GetLength(1); j++)
				//	{
				//		sampleTemp[j] = mortarSamples[i, j];
				//	}
				//	mortarLikelihoodValue += tempDistribution.ProbabilityDensityFunction(sampleTemp);
				//}
				//mortarLikelihoodValue = mortarLikelihoodValue / mortarSamples.GetLength(0);
				//var concreteLikelihoodValue = 0d;
				//for (int i = 0; i < concreteSamples.GetLength(0); i++)
				//{
				//	var sampleTemp = new double[concreteSamples.GetLength(1)];
				//	for (int j = 0; j < concreteSamples.GetLength(1); j++)
				//	{
				//		sampleTemp[j] = concreteSamples[i, j];
				//	}
				//	concreteLikelihoodValue += tempDistribution.ProbabilityDensityFunction(sampleTemp);
				//}
				//concreteLikelihoodValue = concreteLikelihoodValue / concreteSamples.GetLength(0);

				//return cementLikelihoodValue*mortarLikelihoodValue*concreteLikelihoodValue;
				return Math.Log(cementLikelihoodValue);
			}
		}

		private static double[,] LoadCementSamples()
		{
			var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
			var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs\last";
			var netPathName = "cementSamples.txt";

			var filePath = Path.Combine(BasePath, SpecPath);
			var filePathName = Path.Combine(filePath, netPathName);

			using (StreamReader reader = new StreamReader(filePathName))
			{
				// Read the entire file into a string
				string fileContent = reader.ReadToEnd();

				// Split the content into lines
				string[] lines = fileContent.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);

				// Initialize a two-dimensional double array to store the data
				double[,] cementSamples = new double[lines.Length, lines[0].Split('\t').Length]; // Assuming tab-separated values

				for (int i = 0; i < lines.Length; i++)
				{
					string[] values = lines[i].Split('\t'); // Assuming tab-separated values

					for (int j = 0; j < values.Length; j++)
					{
						// Parse each value and store it in the array
						if (double.TryParse(values[j], out double parsedValue))
						{
							cementSamples[i, j] = parsedValue;
						}
						else
						{
							Console.WriteLine($"Error parsing value at row {i + 1}, column {j + 1}");
							// Handle parsing error as needed
						}
					}
				}
				return cementSamples;
			}
		}
		//private static double[,] LoadMortarSamples()
		//{
		//	var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
		//	var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs\last";
		//	var netPathName = "mortarSamples.txt";

		//	var filePath = Path.Combine(BasePath, SpecPath);
		//	var filePathName = Path.Combine(filePath, netPathName);

		//	using (StreamReader reader = new StreamReader(filePathName))
		//	{
		//		// Read the entire file into a string
		//		string fileContent = reader.ReadToEnd();

		//		// Split the content into lines
		//		string[] lines = fileContent.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);

		//		// Initialize a two-dimensional double array to store the data
		//		double[,] mortarSamples = new double[lines.Length, lines[0].Split('\t').Length]; // Assuming tab-separated values

		//		for (int i = 0; i < lines.Length; i++)
		//		{
		//			string[] values = lines[i].Split('\t'); // Assuming tab-separated values

		//			for (int j = 0; j < values.Length; j++)
		//			{
		//				// Parse each value and store it in the array
		//				if (double.TryParse(values[j], out double parsedValue))
		//				{
		//					mortarSamples[i, j] = parsedValue;
		//				}
		//				else
		//				{
		//					Console.WriteLine($"Error parsing value at row {i + 1}, column {j + 1}");
		//					// Handle parsing error as needed
		//				}
		//			}
		//		}
		//		return mortarSamples;
		//	}
		//}
		//private static double[,] LoadConcreteSamples()
		//{
		//	var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
		//	var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs\last";
		//	var netPathName = "concreteSamples.txt";

		//	var filePath = Path.Combine(BasePath, SpecPath);
		//	var filePathName = Path.Combine(filePath, netPathName);

		//	using (StreamReader reader = new StreamReader(filePathName))
		//	{
		//		// Read the entire file into a string
		//		string fileContent = reader.ReadToEnd();

		//		// Split the content into lines
		//		string[] lines = fileContent.Split(new char[] { '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);

		//		// Initialize a two-dimensional double array to store the data
		//		double[,] concreteSamples = new double[lines.Length, lines[0].Split('\t').Length]; // Assuming tab-separated values

		//		for (int i = 0; i < lines.Length; i++)
		//		{
		//			string[] values = lines[i].Split('\t'); // Assuming tab-separated values

		//			for (int j = 0; j < values.Length; j++)
		//			{
		//				// Parse each value and store it in the array
		//				if (double.TryParse(values[j], out double parsedValue))
		//				{
		//					concreteSamples[i, j] = parsedValue;
		//				}
		//				else
		//				{
		//					Console.WriteLine($"Error parsing value at row {i + 1}, column {j + 1}");
		//					// Handle parsing error as needed
		//				}
		//			}
		//		}
		//		return concreteSamples;
		//	}
		//}

		private static double[,] GenerateCementSamples(int numParameterSamples, IProbabilityDistributionSampler cementSampler)
		{
			var cementSamples = cementSampler.GenerateSamples(numParameterSamples);
			var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
			var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs";
			var netPathName = "cementSamples.txt";

			var filePath = Path.Combine(BasePath, SpecPath);
			var filePathName = Path.Combine(filePath, netPathName);

			// Create a StreamWriter to write to the file
			using (StreamWriter writer = new StreamWriter(filePathName))
			{
				int rows = cementSamples.GetLength(0);
				int cols = cementSamples.GetLength(1);

				// Loop through the array and write each element to the file
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						writer.Write(cementSamples[i, j]);

						// Add a tab or comma to separate values (you can choose any delimiter)
						if (j < cols - 1)
							writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
					}

					// Add a newline character after each row
					writer.WriteLine();
				}
			}
			return cementSamples;
		}

		//private static double[,] GenerateMortarSamples(int numParameterSamples, IProbabilityDistributionSampler mortarSampler)
		//{
		//	var mortarSamples = mortarSampler.GenerateSamples(numParameterSamples);
		//	var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
		//	var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs";
		//	var netPathName = "mortarSamples.txt";

		//	var filePath = Path.Combine(BasePath, SpecPath);
		//	var filePathName = Path.Combine(filePath, netPathName);

		//	// Create a StreamWriter to write to the file
		//	using (StreamWriter writer = new StreamWriter(filePathName))
		//	{
		//		int rows = mortarSamples.GetLength(0);
		//		int cols = mortarSamples.GetLength(1);

		//		// Loop through the array and write each element to the file
		//		for (int i = 0; i < rows; i++)
		//		{
		//			for (int j = 0; j < cols; j++)
		//			{
		//				writer.Write(mortarSamples[i, j]);

		//				// Add a tab or comma to separate values (you can choose any delimiter)
		//				if (j < cols - 1)
		//					writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
		//			}

		//			// Add a newline character after each row
		//			writer.WriteLine();
		//		}
		//	}
		//	return mortarSamples;
		//}

		//private static double[,] GenerateConcreteSamples(int numParameterSamples, IProbabilityDistributionSampler concreteSampler)
		//{
		//	var concreteSamples = concreteSampler.GenerateSamples(numParameterSamples);
		//	var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
		//	var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs";
		//	var netPathName = "concreteSamples.txt";

		//	var filePath = Path.Combine(BasePath, SpecPath);
		//	var filePathName = Path.Combine(filePath, netPathName);

		//	// Create a StreamWriter to write to the file
		//	using (StreamWriter writer = new StreamWriter(filePathName))
		//	{
		//		int rows = concreteSamples.GetLength(0);
		//		int cols = concreteSamples.GetLength(1);

		//		// Loop through the array and write each element to the file
		//		for (int i = 0; i < rows; i++)
		//		{
		//			for (int j = 0; j < cols; j++)
		//			{
		//				writer.Write(concreteSamples[i, j]);

		//				// Add a tab or comma to separate values (you can choose any delimiter)
		//				if (j < cols - 1)
		//					writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
		//			}

		//			// Add a newline character after each row
		//			writer.WriteLine();
		//		}
		//	}
		//	return concreteSamples;
		//}

		private static (double[,], double[,]) GenerateHyperparameterAndNewParameterSamples(int numHyperparameterSamples, int numNewParameterSamples, IProbabilityDistributionSampler hyperparameterSampler)
		{
			var hyperparameterSamples = hyperparameterSampler.GenerateSamples(numHyperparameterSamples);

			var mixtureComponents = new MultivariateUniformDistribution[numHyperparameterSamples];
			for (int j = 0; j < numHyperparameterSamples; j++)
			{
				mixtureComponents[j] = new MultivariateUniformDistribution(new double[] { hyperparameterSamples[j, 0], hyperparameterSamples[j, 2], hyperparameterSamples[j, 4] }, new double[] { hyperparameterSamples[j, 0] + hyperparameterSamples[j, 1], hyperparameterSamples[j, 2] + hyperparameterSamples[j, 3], hyperparameterSamples[j, 4] + hyperparameterSamples[j, 5] });
			}
			var mixture = new MultivariateMixture<MultivariateUniformDistribution>(mixtureComponents);
			var newParameterSamples = To2D(mixture.Generate(numNewParameterSamples));

			var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
			var SpecPath = @"MsolveOutputs\ConcreteHierarchicalBayesianOutputs\last";
			var netPathName = "hyperparameterSamples.txt";

			var filePath = Path.Combine(BasePath, SpecPath);
			var filePathName = Path.Combine(filePath, netPathName);

			// Create a StreamWriter to write to the file
			using (StreamWriter writer = new StreamWriter(filePathName))
			{
				int rows = hyperparameterSamples.GetLength(0);
				int cols = hyperparameterSamples.GetLength(1);

				// Loop through the array and write each element to the file
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						writer.Write(hyperparameterSamples[i, j]);

						// Add a tab or comma to separate values (you can choose any delimiter)
						if (j < cols - 1)
							writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
					}

					// Add a newline character after each row
					writer.WriteLine();
				}
			}

			netPathName = "newParameterSamples.txt";

			filePath = Path.Combine(BasePath, SpecPath);
			filePathName = Path.Combine(filePath, netPathName);

			// Create a StreamWriter to write to the file
			using (StreamWriter writer = new StreamWriter(filePathName))
			{
				int rows = newParameterSamples.GetLength(0);
				int cols = newParameterSamples.GetLength(1);

				// Loop through the array and write each element to the file
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						writer.Write(newParameterSamples[i, j]);

						// Add a tab or comma to separate values (you can choose any delimiter)
						if (j < cols - 1)
							writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
					}

					// Add a newline character after each row
					writer.WriteLine();
				}
			}
			return (hyperparameterSamples, newParameterSamples);
		}

		private static double[,] To2D(double[][] source)
		{
			int FirstDim = source.Length;
			int SecondDim = source[0].Length;

			var result = new double[FirstDim, SecondDim];
			for (int i = 0; i < FirstDim; i++)
			{
				var tempSecondDim = source[i].Length;
				for (int j = 0; j < tempSecondDim; j++)
					result[i, j] = source[i][j];
				if (tempSecondDim < SecondDim)
					for (int j = tempSecondDim; j < SecondDim; j++)
					{
						result[i, j] = 0;
					}
			}

			return result;
		}

		private static IProbabilityDistributionSampler InitializeSimulation()
		{
			var model = new ConcreteHierarchicalModels();
			model.InitializeModels();
			//var cementResponse = model.FormulateCementProblem(new double[] { 10, 1, 0.1 });
			//var mortarResponse = model.FormulateMortarProblem(new double[] { 10, 1, 0.1 });
			//var concreteResponse = model.FormulateConcreteProblem(new double[] { 10, 1, 0.1 });

			//var mean_hyperprior = new MultivariateUniformDistribution(new double[] { 10 }, new double[] { 30 });
			//var var_hyperprior = new MultivariateUniformDistribution(new double[] { 1 }, new double[] { 10 });
			//var lower_bound_hyperprior = new MultivariateUniformDistribution(new double[] { 0 }, new double[] { 1 });
			//var upper_bound_hyperprior = new MultivariateUniformDistribution(new double[] { 1 }, new double[] { 2 });
			//var hyperprior = new Dictionary<string, ISampleableDistribution<double[]>>();
			//hyperprior.Add("mean", mean_hyperprior);
			//hyperprior.Add("var", var_hyperprior);
			//var measurementValues = new double[3] { 9.012012689089184, 4.440024343431458, 3.1314475332399594 };
			//var measurementError = new double[3] { 0.1, 0.1, 0.1 };
			//var bayesianInstance = new HierarchicalBayesianUpdate(model.FormulateProblem, prior, hyperprior, measurementValues, measurementError);
			//var sampler = new TransitionalMarkovChainMonteCarlo(3, bayesianInstance.LikelihoodModel, bayesianInstance.PriorDistributionGenerator, scalingFactor: 0.2);
			//var samples = sampler.GenerateSamples(10000);
			//var mean = samples.Mean(0);
			//var std = samples.StandardDeviation();
			//Assert.True(Math.Abs(mean - 1.6) < 0.05);
			//Assert.True(Math.Abs(std[0] - 0.45) < 0.1);
			//var model = new ConcreteHierarchicalModels();

			//(new double[] { 29.9, 2.99, 0.299 }, new double[] { 30, 3, 0.3 });
			//(new double[] { 0, 0, 0 }, new double[] { 0.1, 0.01, 0.001 });
			//(new double[] { 0, 0, 0 }, new double[] { 30, 3, 0.3 });
			//(new double[] { 18, 0.5, 0.08 }, new double[] { 18.01, 0.501, 0.0801 });

			var cementPrior = new MultivariateUniformDistribution(new double[] { 0, 0, 0 }, new double[] { 30, 3, 0.3 });
			var cementMeasurementValues = new double[1] { 0.000210 };
			var cementMeasurementError = new double[1] { 0.00000105 * 0.00000105 };
			var cementBayesianInstance = new BayesianUpdate(model.FormulateCementProblem, cementPrior, cementMeasurementValues, cementMeasurementError);
			var cementSampler = new TransitionalMarkovChainMonteCarlo(3, cementBayesianInstance.LikelihoodModel, cementBayesianInstance.PriorDistributionGenerator, scalingFactor: 0.2);

			//var mortarPrior = new MultivariateUniformDistribution(new double[] { 4.99, 0.499, 0.299 }, new double[] { 5,  0.5, 0.3 });
			//var mortarMeasurementValues = new double[1] { 0.0019722238 };
			//var mortarMeasurementError = new double[1] { 0.00000986111 * 0.00000986111 };
			//var mortarBayesianInstance = new BayesianUpdate(model.FormulateMortarProblem, mortarPrior, mortarMeasurementValues, mortarMeasurementError);
			//var mortarSampler = new TransitionalMarkovChainMonteCarlo(3, mortarBayesianInstance.LikelihoodModel, mortarBayesianInstance.PriorDistributionGenerator, scalingFactor: 0.2);

			//var concretePrior = new MultivariateUniformDistribution(new double[] { 0, 0, 0 }, new double[] { 30, 3, 0.3 });
			//var concreteMeasurementValues = new double[1] { 0.0142 };
			//var concreteMeasurementError = new double[1] { 0.000071 * 0.000071 };
			//var concreteBayesianInstance = new BayesianUpdate(model.FormulateConcreteProblem, concretePrior, concreteMeasurementValues, concreteMeasurementError);
			//var concreteSampler = new TransitionalMarkovChainMonteCarlo(3, concreteBayesianInstance.LikelihoodModel, concreteBayesianInstance.PriorDistributionGenerator, scalingFactor: 0.2);

			return cementSampler;
		}

	}
}
