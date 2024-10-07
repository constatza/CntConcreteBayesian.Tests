using System;
using System.IO;

using MGroup.Constitutive.Structural.MachineLearning;
using MGroup.MachineLearning.Preprocessing;
using MGroup.MachineLearning.TensorFlow.NeuralNetworks;
using static Tensorflow.KerasApi;

using Xunit;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.MachineLearning.Tests
{
    public static class NeuralNetworkMaterial3DBuilderTest
	{
        [Fact]
		public static void RunTest()
        {
			var neuralNetworkmaterial = new NeuralNetworkMaterialBuilder();
		    neuralNetworkmaterial.GenerateStrainStressData();
			var trainedNetwork = neuralNetworkmaterial.TrainNeuralNetwork();

			var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
			var SpecPath = @"MsolveOutputs\NeuralNetworks";
			var netPathName = "network_architecture";
			var weightsPathName = "trained_weights";
			var normalizationPathName = "normalization";

			var pathName = Path.Combine(BasePath, SpecPath);

			var InputExtension = Path.GetExtension(netPathName);
			var InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(netPathName));
			var netPathFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);

			InputExtension = Path.GetExtension(weightsPathName);
			InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(weightsPathName));
			var weightsPathFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);

			InputExtension = Path.GetExtension(normalizationPathName);
			InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(normalizationPathName));
			var normalizationPathFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);


			var neuralNetwork = new FeedForwardNeuralNetwork();

			neuralNetwork.LoadNetwork(netPathFile, weightsPathFile, normalizationPathFile);

			//var neuralNetworkMaterial1 = new NeuralNetworkMaterial3D(trainedNetwork, new double[0]);

			////neuralNetworkMaterial1.MaterialParameters = new double[] { 20 };
			//for (int i = 0; i < 10; i++)
			//{
			//	var strains = new double[6] { 0.002, 0, 0, 0, 0, 0 };
			//	neuralNetworkMaterial1.UpdateConstitutiveMatrixAndEvaluateResponse(strains);
			//	var stresses = neuralNetworkMaterial1.Stresses;
			//	var constitutiveMatrix = neuralNetworkMaterial1.ConstitutiveMatrix;
			//	neuralNetworkMaterial1.CreateState();
			//}

			var neuralNetworkMaterial = new NeuralNetworkMaterial3D(neuralNetwork, new double[0]);
			var elasticMaterial = new ElasticMaterial3D(20, 0.2);
			//neuralNetworkMaterial.MaterialParameters = new double[] { 20 };
			var stressesNeuralNetwork = new double[10,6];
			var constitutiveNeuralNetwork = new Matrix[10];
			var stressesElastic = new double[10, 6];
			var constitutiveElastic = new Matrix[10];
			for (int i = 0; i < 10; i++)
			{
				var strains = new double[6] { 0.001, 0.001, 0.001, 0.001, 0.001, 0.001 };
				elasticMaterial.UpdateConstitutiveMatrixAndEvaluateResponse(strains);
				neuralNetworkMaterial.UpdateConstitutiveMatrixAndEvaluateResponse(strains);
				for (int j = 0; j < 6; j++)
				{
					stressesNeuralNetwork[i,j] = neuralNetworkMaterial.Stresses[j];
					stressesElastic[i, j] = elasticMaterial.Stresses[j];
				}
				constitutiveNeuralNetwork[i] = (Matrix)neuralNetworkMaterial.ConstitutiveMatrix;
				constitutiveElastic[i] = (Matrix)elasticMaterial.ConstitutiveMatrix;
				neuralNetworkMaterial.CreateState();
				elasticMaterial.CreateState();
			}
		}
	}
}
