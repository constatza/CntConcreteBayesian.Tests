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
    public class ConcreteMonteCarloTest
    {
        [Fact]
        public static void RunSimulation()
        {
            var model = new ConcreteHierarchicalModels();
            model.InitializeModels();
            Func<double[], double[]> cementModel = model.FormulateCementProblem;
            MultivariateUniformDistribution multivariateUniformDistribution = new MultivariateUniformDistribution(new double[] { 0, 0, 0 }, new double[] { 30, 3, 0.3 });

            var cementSampler = new MonteCarlo(multivariateUniformDistribution.Dimension, cementModel, multivariateUniformDistribution);

            var numParameterSamples = 10;

            var cementSamples = cementSampler.GenerateSamples(numParameterSamples);
            //var mortarSamples = GenerateMortarSamples(numParameterSamples, mortarSampler);
            //var concreteSamples = GenerateConcreteSamples(numParameterSamples, concreteSampler);
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

    }
}
