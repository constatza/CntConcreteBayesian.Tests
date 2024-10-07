using Xunit;
using System.Linq;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.MultiscaleAnalysis;
using MiMsolve.SolutionStrategies;
using System.IO;
using MGroup.LinearAlgebra.Commons;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization.Entities;
using System.Collections.Generic;
using MGroup.Solvers.Direct;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Discretization.Dofs;

namespace MGroup.Multiscale.Tests.RveTemplates.Tests.CntReinforcedElasticNanocompositeTest
{
	public class CntWithoutMicrostructureTest
	{
		[Fact]
		public static void CntWithoutMicrostructureLinearExample()
		{
			string resultsPath = "..\\..\\RveTemplates\\Input\\ExpectedResults\\CntsWithoutMicrostructureLinearResults.txt";
			double[] expectedDisplacements = PrintUtilities.ReadVector(resultsPath);

			int numberOfCnts = 50;
			var matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
			var homogeneousRveBuilder1 =
				new CntReinforcedElasticNanocomposite(numberOfCnts, matrixMaterial);
			homogeneousRveBuilder1.readFromText = true;

			//Boundary Conditions
			var tuple1 = homogeneousRveBuilder1.GetModelAndBoundaryNodes();
			var model = tuple1.Item1;
			var boundaryNodes = tuple1.Item2;

			var loads = new List<NodalLoad[]>();
			var scaleFactor = 1;
			foreach (var nodeID in model.NodesDictionary.Keys)
			{
				bool hasBoundaryCondition = boundaryNodes.TryGetValue(nodeID, out INode boundaryNode);
				if (!hasBoundaryCondition)
				{
					loads.Add(new[]
					{
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationY, amount: -10d * scaleFactor),
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationX, amount: 5d * scaleFactor),
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationZ, amount: 10d * scaleFactor)
					});
					model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(new List<INodalDisplacementBoundaryCondition>(), loads[loads.Count-1]));
				}
				if (loads.Count == 10)
				{
					break;
				}
			}

			foreach (var boundaryNode in boundaryNodes.Values)
			{
				var constraints = new List<INodalDisplacementBoundaryCondition>();
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationX, amount: 0d));
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationY, amount: 0d));
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationZ, amount: 0d));
				model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, new NodalLoad[] { }));
			}

			var solverFactory = new SkylineSolver.Factory();
			var algebraicModel = solverFactory.BuildAlgebraicModel(model);
			var solver = solverFactory.BuildSolver(algebraicModel);

			var problem = new ProblemStructural(model, algebraicModel);
			var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, linearAnalyzer);

			staticAnalyzer.Initialize();
			staticAnalyzer.Solve();
			var solution = solver.LinearSystem.Solution.SingleVector;

			Assert.True(AreDisplacementsSame(expectedDisplacements, solution.CopyToArray()));
		}

		[Fact]
		public static void CntWithoutMicrostructureNonLinearExample()
		{
			string resultsPath = "..\\..\\RveTemplates\\Input\\ExpectedResults\\CntsWithoutMicrostructureNonLinearResults.txt";
			double[] expectedDisplacements = PrintUtilities.ReadVector(resultsPath);
			var expectedResultsList = new List<double[]>();
			int count = 0;
			expectedResultsList.Add(new double[30]);
			foreach (var disp in expectedDisplacements)
			{
				if (count < 30)
				{
					expectedResultsList[expectedResultsList.Count - 1][count] = disp;
				}
				else
				{
					count = 0;
					expectedResultsList.Add(new double[30]);
					expectedResultsList[expectedResultsList.Count - 1][count] = disp;
				}
				count++;
			}

			int numberOfCnts = 50;
			var matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
			var homogeneousRveBuilder1 =
				new CntReinforcedElasticNanocomposite(numberOfCnts, matrixMaterial);
			homogeneousRveBuilder1.readFromText = true;

			//Boundary Conditions
			var tuple1 = homogeneousRveBuilder1.GetModelAndBoundaryNodes();
			var model = tuple1.Item1;
			var boundaryNodes = tuple1.Item2;

			var loads = new List<NodalLoad[]>();
			var scaleFactor = 1;
			foreach (var nodeID in model.NodesDictionary.Keys)
			{
				bool hasBoundaryCondition = boundaryNodes.TryGetValue(nodeID, out INode boundaryNode);
				if (!hasBoundaryCondition)
				{
					loads.Add(new[]
					{
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationY, amount: -10d * scaleFactor),
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationX, amount: 5d * scaleFactor),
						new NodalLoad(model.NodesDictionary[nodeID], StructuralDof.TranslationZ, amount: 10d * scaleFactor),
					});
					model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(new List<INodalDisplacementBoundaryCondition>(), loads[loads.Count - 1]));
				}
				if (loads.Count == 10)
				{
					break;
				}
			}

			foreach (var boundaryNode in boundaryNodes.Values)
			{
				var constraints = new List<INodalDisplacementBoundaryCondition>();
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationX, amount: 0d));
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationY, amount: 0d));
				constraints.Add(new NodalDisplacement(boundaryNode, StructuralDof.TranslationZ, amount: 0d));
				model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, new NodalLoad[] { }));
			}

			var solverFactory = new SkylineSolver.Factory();
			var algebraicModel = solverFactory.BuildAlgebraicModel(model);
			var solver = solverFactory.BuildSolver(algebraicModel);

			var problem = new ProblemStructural(model, algebraicModel);
			var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 2)
			{
				ResidualTolerance = 1E-3,
				MaxIterationsPerIncrement = 1000,
				NumIterationsForMatrixRebuild = 1,
			};
			var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();

			var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);

			var watchDofs = new List<(INode node, IDofType dof)>();
			foreach (var nodeID in model.NodesDictionary.Keys)
			{
				bool hasBoundaryCondition = boundaryNodes.TryGetValue(nodeID, out INode boundaryNode);
				if (!hasBoundaryCondition)
				{
					watchDofs.Add((model.NodesDictionary[nodeID], StructuralDof.TranslationX));
					watchDofs.Add((model.NodesDictionary[nodeID], StructuralDof.TranslationY));
					watchDofs.Add((model.NodesDictionary[nodeID], StructuralDof.TranslationZ));
				}
				if (watchDofs.Count == 30)
				{
					break;
				}
			}
			loadControlAnalyzer.TotalDisplacementsPerIterationLog = new TotalDisplacementsPerIterationLog(watchDofs, algebraicModel);
			
			staticAnalyzer.Initialize();
			staticAnalyzer.Solve();

			Assert.True(AreDisplacementsSame(expectedResultsList, loadControlAnalyzer.TotalDisplacementsPerIterationLog,1E-10));
		}

		public static bool AreDisplacementsSame(double[] expectedValues,
			double[] computedValues, double tol = 1E-10)
		{
			var comparer = new ValueComparer(tol);
			for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
			{
				if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
				{
					return false;
				}
			}
			return true;
		}

		public static bool AreDisplacementsSame(IReadOnlyList<double[]> expectedDisplacements,
			TotalDisplacementsPerIterationLog computedDisplacements, double tolerance)
		{
			var comparer = new ValueComparer(tolerance);
			for (var iter = 0; iter < expectedDisplacements.Count; ++iter)
			{
				for (var i = 0; i < expectedDisplacements[iter].Length; ++i)
				{
					var expected = expectedDisplacements[iter][i];
					(var node, var dof) = computedDisplacements.WatchDofs[i];
					var computed = computedDisplacements.GetTotalDisplacement(iter, node, dof);

					if (!comparer.AreEqual(expected, computed)) return false;
				}
			}
			return true;
		}
	}
}
