using System;
using System.IO;

using MGroup.Constitutive.Structural.MachineLearning;
using MGroup.MachineLearning.Preprocessing;
using MGroup.MachineLearning.TensorFlow.NeuralNetworks;
using static Tensorflow.KerasApi;
using MGroup.CntConcreteBayesian.Tests.Models;
using Xunit;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.CntConcreteBayesian.Tests.Integration
{
    public static class NeuralNetworkMaterial3DBuilderTest
	{
        [Fact]
		public static void RunTest()
        {
			var neuralNetworkmaterial = new NeuralNetworkMaterialBuilder();
		    neuralNetworkmaterial.GenerateStrainStressData();
			var trainedNetwork = neuralNetworkmaterial.TrainNeuralNetwork();
		}
	}
}
