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
using System.Reflection;

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
            Func<double[], double[]> mortarModel = model.FormulateMortarProblem;
            Func<double[], double[]> concreteModel = model.FormulateConcreteProblem;
            MultivariateUniformDistribution multivariateUniformDistribution = new MultivariateUniformDistribution(new double[] { 0, 0, 0 }, new double[] { 30, 3, 0.3 });

            var cementSampler = new MonteCarlo(multivariateUniformDistribution.Dimension, cementModel, multivariateUniformDistribution);
            var mortarSampler = new MonteCarlo(multivariateUniformDistribution.Dimension, mortarModel, multivariateUniformDistribution);
            var concreteSampler = new MonteCarlo(multivariateUniformDistribution.Dimension, concreteModel, multivariateUniformDistribution);

            var numParameterSamples = 10;

            var cementSamples = GenerateCementSamples(numParameterSamples, cementSampler);
            var mortarSamples = GenerateMortarSamples(numParameterSamples, mortarSampler);
            var concreteSamples = GenerateConcreteSamples(numParameterSamples, concreteSampler);
        }

        private static double[,] GenerateCementSamples(int numParameterSamples, IProbabilityDistributionSampler cementSampler)
        {
            var cementSamples = cementSampler.GenerateSamples(numParameterSamples);
            var BasePath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            var folderName = "DataFiles";
            var inputFileName = "cementSamples.txt";
            var filePathName = Path.Combine(BasePath, folderName, inputFileName);

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

        private static double[,] GenerateMortarSamples(int numParameterSamples, IProbabilityDistributionSampler mortarSampler)
        {
            var mortarSamples = mortarSampler.GenerateSamples(numParameterSamples);
            var BasePath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            var folderName = "DataFiles";
            var inputFileName = "mortarSamples.txt";
            var filePathName = Path.Combine(BasePath, folderName, inputFileName);

            // Create a StreamWriter to write to the file
            using (StreamWriter writer = new StreamWriter(filePathName))
            {
                int rows = mortarSamples.GetLength(0);
                int cols = mortarSamples.GetLength(1);

                // Loop through the array and write each element to the file
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(mortarSamples[i, j]);

                        // Add a tab or comma to separate values (you can choose any delimiter)
                        if (j < cols - 1)
                            writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
                    }

                    // Add a newline character after each row
                    writer.WriteLine();
                }
            }
            return mortarSamples;
        }

        private static double[,] GenerateConcreteSamples(int numParameterSamples, IProbabilityDistributionSampler concreteSampler)
        {
            var concreteSamples = concreteSampler.GenerateSamples(numParameterSamples);
            var BasePath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            var folderName = "DataFiles";
            var inputFileName = "concreteSamples.txt";
            var filePathName = Path.Combine(BasePath, folderName, inputFileName);

            // Create a StreamWriter to write to the file
            using (StreamWriter writer = new StreamWriter(filePathName))
            {
                int rows = concreteSamples.GetLength(0);
                int cols = concreteSamples.GetLength(1);

                // Loop through the array and write each element to the file
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        writer.Write(concreteSamples[i, j]);

                        // Add a tab or comma to separate values (you can choose any delimiter)
                        if (j < cols - 1)
                            writer.Write("\t"); // Use "\t" for tab separation, or "," for comma separation
                    }

                    // Add a newline character after each row
                    writer.WriteLine();
                }
            }
            return concreteSamples;
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
