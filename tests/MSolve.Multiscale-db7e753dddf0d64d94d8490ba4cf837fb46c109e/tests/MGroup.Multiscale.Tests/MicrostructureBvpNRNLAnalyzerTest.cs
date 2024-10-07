using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
using MGroup.Multiscale.SupportiveClasses;
using MGroup.Multiscale.Tests.ExampleModels;
using MGroup.Multiscale.Tests.FEMpartB;

using Xunit;

namespace MGroup.Multiscale.Tests
{
	public static class MicrostructureBvpNRNLAnalyzerTest
	{
		[Fact]
		public static void CheckMicrostructureBvpNRNLAnalyzer()
		{
			string results_file1 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MSmicroBvpAnalyzerTest\\U_sunol_micro_3.txt";
			string results_file2 = "..\\..\\..\\MGroup.Multiscale.Tests\\InputFiles\\MSmicroBvpAnalyzerTest\\U_sunol_micro_6.txt";
			double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
			double[] displacements2ndncrement = PrintUtilities.ReadVector(results_file2);

			(double[] uInitialFreeDOFs_state1, double[] uInitialFreeDOFs_state2) = NRNLAnalyzerDevelopExample.SolveDisplLoadsExample();


			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1));
			Assert.True(ComparisonMetods.AreDisplacementsSame(displacements2ndncrement, uInitialFreeDOFs_state2));
		}
	}
}
