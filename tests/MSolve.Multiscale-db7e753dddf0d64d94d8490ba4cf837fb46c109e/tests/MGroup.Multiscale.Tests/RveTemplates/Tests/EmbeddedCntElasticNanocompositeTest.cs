namespace MGroup.Multiscale.Tests.RveTemplates.Tests.CntReinforcedElasticNanocompositeTest
{
	using Xunit;
	using System.Linq;
	using MGroup.Constitutive.Structural.Continuum;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.MultiscaleAnalysis;
	using MiMsolve.SolutionStrategies;
	using System.IO;
	using MGroup.LinearAlgebra.Commons;

	public class EmbeddedCntElasticNanocompositeTest
	{

		[Fact]
		public static void RunTest()
		{
			int numberOfCnts = 1;
			int solutions = 1;
			int increments_per_solution = 3;
			double[][] Input = new double[solutions * increments_per_solution][];
			double[][] Output = new double[solutions * increments_per_solution][];

			var matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
			var inclusionMaterial = new ElasticMaterial3D(youngModulus: 3.5 * 1000, poissonRatio: 0.4);
			var homogeneousRveBuilder1 =
				new EmbeddedCntElasticNanocomposite(numberOfCnts, matrixMaterial, inclusionMaterial); //{ K_el = 20, K_pl = 2, T_max = 0.2, };
			homogeneousRveBuilder1.readFromText = false;
			IContinuumMaterial3D microstructure3 = new Microstructure3D<SymmetricCscMatrix>(homogeneousRveBuilder1, false, 1, new SuiteSparseSolverPrefernce());

			for (int num_solution = 0; num_solution < solutions; num_solution++)
			{
				var maxstrain = 0.01;
				//var MacroStrain = new double[6] { maxstrain, -0.2 * maxstrain, -0.2 * maxstrain, 0.0, 0.0, 0.0 };
				var MacroStrain = new double[6] { maxstrain, 0.0, 0.0, 0.0, 0.0, 0.0 };

				var IncrMacroStrainPrevious = new double[6];
				microstructure3 = new Microstructure3D<SymmetricCscMatrix>(homogeneousRveBuilder1, false, 1, new SuiteSparseSolverPrefernce());
				for (int i = 0; i < increments_per_solution; i++)
				{
					var IncrMacroStrain = new double[6];
					for (int ii = 0; ii < 6; ii++)
					{
						IncrMacroStrain[ii] = IncrMacroStrainPrevious[ii] + (MacroStrain[ii] / increments_per_solution);
						IncrMacroStrainPrevious[ii] = IncrMacroStrain[ii];
					}

					double[] stresses = microstructure3.UpdateConstitutiveMatrixAndEvaluateResponse(new double[6] { IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] });

					double[] IncrMacroStress = new double[6] { stresses[0], stresses[1], stresses[2], stresses[3], stresses[4], stresses[5] };

					microstructure3.CreateState();
					Input[num_solution * increments_per_solution + i] = new double[6] { IncrMacroStrain[0], IncrMacroStrain[1], IncrMacroStrain[2], IncrMacroStrain[3], IncrMacroStrain[4], IncrMacroStrain[5] };
					Output[num_solution * increments_per_solution + i] = new double[6] { IncrMacroStress[0], IncrMacroStress[1], IncrMacroStress[2], IncrMacroStress[3], IncrMacroStress[4], IncrMacroStress[5] }; //homogeneousRveBuilder1.K_el, homogeneousRveBuilder1.K_pl, homogeneousRveBuilder1.T_max,
				}
			}

			double[][] expectedStresses =
			{
				new double[] { 0.114084979796056, 0.0505117120828458, 0.0504880927825645, -4.93933294792444E-05, -0.000208109197775651, -1.82646042480646E-05 },
				new double[] { 0.227567405034224, 0.100787865230969, 0.101071577520994, -2.74541675661799E-05, -0.000394223946003412, 0.000424039268373175 },
				new double[] { 0.343348543430872, 0.150912280313983, 0.150323829694954, -0.000231432098355594, -0.000297071304919592, -0.000574232847650175 },
				new double[] { 0.4554310736263, 0.201244618073806, 0.201115119909942, -0.000227756203607088, -0.00049816369599425, -0.00010106368807711 },
				new double[] { 0.568973804213001, 0.251542314867985, 0.251413101479795, -0.000290378764270461, -0.000595226143332466, -0.000119353161098844 }
			};

			Assert.True(AreStressesSame(expectedStresses, Output));
		}

		public static bool AreStressesSame(double[][] expectedValues, double[][] computedValues, double tol = 1E-9)
		{
			var comparer = new ValueComparer(tol);
			for (int i1 = 0; i1 < expectedValues.Length; i1++)
			{
				for (int i2 = 0; i2 < expectedValues[i1].Length; i2++)
				{
					if (!comparer.AreEqual(expectedValues[i1][i2], computedValues[i1][i2]))
					{
						return false;
					}
				}
			}
			return true;
		}
	}
}
