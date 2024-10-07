using System;
using MGroup.MachineLearning.TensorFlow.NeuralNetworks;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;
using Xunit;
using System.Reflection;
using MGroup.Multiscale;
using MGroup.Constitutive.Structural.Cohesive;
using MGroup.MSolve.MultiscaleAnalysis;
using MiMsolve.SolutionStrategies;

namespace MGroup.Constitutive.Structural.MachineLearning.Tests
{
	public static class NeuralNetworkMaterial3DTest
	{
		[Fact]
		public static void RunTest()
		{
			// these files are used to generate an already trained FeedForwardNeuralNetwork which was created using strain-stress pairs from an ElasticMaterial3D(youngModulus:20, poissonRatio:0.2)
			string initialPath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
			var folderName = "SavedFiles";
			var netPathName = "mortar_network_architecture";
			netPathName = Path.Combine(initialPath, folderName, netPathName);
			var weightsPathName = "mortar_trained_weights";
			weightsPathName = Path.Combine(initialPath, folderName, weightsPathName);
			var normalizationPathName = "mortar_normalization";
			normalizationPathName = Path.Combine(initialPath, folderName, normalizationPathName);

			var neuralNetwork = new FeedForwardNeuralNetwork();
			neuralNetwork.LoadNetwork(netPathName, weightsPathName, normalizationPathName);

			var neuralNetworkMaterial = new NeuralNetworkMaterial3D(neuralNetwork, new double[0]);
			var elasticMaterial = new ElasticMaterial3D(1353000, 0.3);

			CheckNeuralNetworkMaterialAccuracy(neuralNetwork);
		}

		private static void CheckNeuralNetworkMaterialAccuracy(FeedForwardNeuralNetwork neuralNetwork)
		{
			var maxParameterValues = new double[3] { 30, 3, 0.3 };
			var minParameterValues = new double[3] { 0.1, 0.01, 0.001 };

			//analysis properties
			var numOfSolutions = 1;
			var incrementsPerSolution = 20;
			var maxLimitStrain = new double[6] { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
			var minLimitStrain = new double[6] { -0.01, -0.01, -0.01, -0.01, -0.01, -0.01 };

			var strains = new double[numOfSolutions * incrementsPerSolution, 6];
			var stressesNeuralNetwork = new double[numOfSolutions * incrementsPerSolution, 6];
			var constitutiveNeuralNetwork = new Matrix[numOfSolutions * incrementsPerSolution];
			var stressesPlastic = new double[numOfSolutions * incrementsPerSolution, 6];
			var constitutivePlastic = new Matrix[numOfSolutions * incrementsPerSolution];
			var TotalMacroStrain = new double[maxLimitStrain.Length];
			var rnd = new Random(5);
			for (int k = 0; k < numOfSolutions; k++)
			{
				var checkNaN = false;
				var parameterValues = new double[maxParameterValues.Length];
				for (int i = 0; i < maxParameterValues.Length; i++)
				{
					parameterValues[i] = rnd.NextDouble() * (maxParameterValues[i] - minParameterValues[i]) + minParameterValues[i];
				}
				//matrixMaterial = new ElasticMaterial3D(parameterValues[0], parameterValues[1]);
				var matrixMaterial = new DruckerPrager3DFunctional(20, 0.2, 30, 30, 0.025, x => 0.025 + 2 * x);
				var cohesiveMaterial = new BondSlipMaterial(parameterValues[0], parameterValues[1], 100.0, parameterValues[2], new double[2], new double[2], 1e-3);
				var numOfInclusions = 260;
				var rveBuilder = new CntReinforcedElasticNanocomposite(numOfInclusions, matrixMaterial, cohesiveMaterial: cohesiveMaterial);
				rveBuilder.readFromText = true;
				//var matrixMaterial = new NeuralNetworkMaterial3D(neuralNetwork, parameterValues);
				//var inclusionMaterial = new ElasticMaterial3D(60, 0.22);
				//var rveBuilder = new GmshCompositeRveBuilder(matrixMaterial, inclusionMaterial, 100, 100, 100, "..\\..\\RveTemplates\\Input\\Continuum\\concrete20wf.msh");
				var microstructure = new Microstructure3D<SymmetricCscMatrix>(rveBuilder, false, 1, new SuiteSparseSolverPrefernce());
				microstructure.MaxIterations = 20;

				var MacroStrain = new double[maxLimitStrain.Length];
				//MacroStrain[0] = 0.0109521616915730; MacroStrain[1] = 0.00794466169241941; MacroStrain[2] =-0.0165482409336200; MacroStrain[3] = 0.0166872384368614; MacroStrain[4] = 0.0192216267081954; MacroStrain[5] = 0.0169534127761485;
				for (int i = 0; i < maxLimitStrain.Length; i++)
				{
					MacroStrain[i] = rnd.NextDouble() * (maxLimitStrain[i] - minLimitStrain[i]) + minLimitStrain[i];
				}
				var IncrMacroStrain = new double[6];
				for (int i = 0; i < incrementsPerSolution; i++)
				{
					for (int ii = 0; ii < 6; ii++) { IncrMacroStrain[ii] = (i + 1) * MacroStrain[ii] / incrementsPerSolution; }
					microstructure.UpdateConstitutiveMatrixAndEvaluateResponse(new double[6] { IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] });

					for (int j = 0; j < 6; j++)
					{
						stressesPlastic[k * incrementsPerSolution + i, j] = microstructure.Stresses[j];
						strains[k * incrementsPerSolution + i, j] = IncrMacroStrain[j];
					}
					constitutivePlastic[k * incrementsPerSolution + i] = (Matrix)microstructure.ConstitutiveMatrix.Copy();
					microstructure.CreateState();
				}

				//var IncrStrain = new double[6];
				var neuralNetworkMaterial = new NeuralNetworkMaterial3D(neuralNetwork, new double[3] { parameterValues[0], parameterValues[1], parameterValues[2] });
				for (int i = 0; i < incrementsPerSolution; i++)
				{
					for (int ii = 0; ii < 6; ii++) { IncrMacroStrain[ii] = MacroStrain[ii] / incrementsPerSolution; }
					neuralNetworkMaterial.UpdateConstitutiveMatrixAndEvaluateResponse(new double[6] { IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] });

					double[] MacroStress = new double[6] { neuralNetworkMaterial.Stresses[0], neuralNetworkMaterial.Stresses[1], neuralNetworkMaterial.Stresses[2], neuralNetworkMaterial.Stresses[3], neuralNetworkMaterial.Stresses[4], neuralNetworkMaterial.Stresses[5] };

					for (int j = 0; j < 6; j++)
					{
						stressesNeuralNetwork[k * incrementsPerSolution + i, j] = neuralNetworkMaterial.Stresses[j];
					}
					constitutiveNeuralNetwork[k * incrementsPerSolution + i] = (Matrix)neuralNetworkMaterial.ConstitutiveMatrix.Copy();
					neuralNetworkMaterial.CreateState();
				}
			}

			var stressDeviation = 0d;
			var constitutiveDeviation = 0d;
			for (int k = 0; k < numOfSolutions; k++)
			{
				for (int i = 0; i < incrementsPerSolution; i++)
				{
					for (int j1 = 0; j1 < 6; j1++)
					{
						stressDeviation += Math.Pow((stressesNeuralNetwork[k * incrementsPerSolution + i, j1] - stressesPlastic[k * incrementsPerSolution + i, j1]), 2) / Math.Pow(stressesPlastic[k * incrementsPerSolution + i, j1], 2);
						for (int j2 = 0; j2 < 6; j2++)
						{
							constitutiveDeviation += Math.Pow((constitutiveNeuralNetwork[k * incrementsPerSolution + i][j1, j2] - constitutivePlastic[k * incrementsPerSolution + i][j1, j2]), 2) / Math.Pow(constitutivePlastic[k * incrementsPerSolution + i][j1, j2], 2);
						}
					}
				}
			}

			stressDeviation = stressDeviation / (incrementsPerSolution * numOfSolutions * 6);
			stressDeviation = Math.Sqrt(stressDeviation);
			constitutiveDeviation = constitutiveDeviation / (incrementsPerSolution * numOfSolutions * 6 * 6);
			constitutiveDeviation = Math.Sqrt(constitutiveDeviation);

			Assert.True(stressDeviation < 1e-6 && constitutiveDeviation < 2e-1);
		}
	}
}
