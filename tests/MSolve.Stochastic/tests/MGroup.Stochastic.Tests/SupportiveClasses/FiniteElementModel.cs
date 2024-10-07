using System;
using System.Collections.Generic;

using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.Solvers.Direct;

namespace MGroup.Stochastic.Tests.SupportiveClasses
{

	public class FiniteElementModel
	{
		public double[] CreateModel(double[] parameters)
		{
			var nodeData = new double[,] {
			{-0.250000,-0.250000,-1.000000},
			{0.250000,-0.250000,-1.000000},
			{-0.250000,0.250000,-1.000000},
			{0.250000,0.250000,-1.000000},
			{-0.250000,-0.250000,-0.500000},
			{0.250000,-0.250000,-0.500000},
			{-0.250000,0.250000,-0.500000},
			{0.250000,0.250000,-0.500000},
			{-0.250000,-0.250000,0.000000},
			{0.250000,-0.250000,0.000000},
			{-0.250000,0.250000,0.000000},
			{0.250000,0.250000,0.000000},
			{-0.250000,-0.250000,0.500000},
			{0.250000,-0.250000,0.500000},
			{-0.250000,0.250000,0.500000},
			{0.250000,0.250000,0.500000},
			{-0.250000,-0.250000,1.000000},
			{0.250000,-0.250000,1.000000},
			{-0.250000,0.250000,1.000000},
			{0.250000,0.250000,1.000000}
		};
			double correction = 10;// +20;

			var elementData = new int[,] {
			{1,8,7,5,6,4,3,1,2},
			{2,12,11,9,10,8,7,5,6},
			{3,16,15,13,14,12,11,9,10},
			{4,20,19,17,18,16,15,13,14}
		};

			var model = new Model();

			model.SubdomainsDictionary.Add(key: 0, new Subdomain(id: 0));

			for (var i = 0; i < nodeData.GetLength(0); i++)
			{
				var nodeId = i + 1;
				model.NodesDictionary.Add(nodeId, new Node(
					id: nodeId,
					x: nodeData[i, 0] + correction,
					y: nodeData[i, 1] + correction,
					z: nodeData[i, 2] + correction));
			}

			for (var i = 0; i < elementData.GetLength(0); i++)
			{
				var nodeSet = new Node[8];
				for (var j = 0; j < 8; j++)
				{
					var nodeID = elementData[i, j + 1];
					nodeSet[j] = (Node)model.NodesDictionary[nodeID];
				}

				var elementFactory = new ContinuumElement3DFactory(new ElasticMaterial3D(parameters[0], 0.3), new TransientAnalysisProperties(1, 0, 0));
				var element = elementFactory.CreateElement(CellType.Hexa8, nodeSet);
				element.ID = i + 1;

				model.ElementsDictionary.Add(element.ID, element);
				model.SubdomainsDictionary[0].Elements.Add(element);
			}

			var constraints = new List<INodalDisplacementBoundaryCondition>();
			for (var i = 1; i < 5; i++)
			{
				constraints.Add(new NodalDisplacement(model.NodesDictionary[i], StructuralDof.TranslationX, amount: 0d));
				constraints.Add(new NodalDisplacement(model.NodesDictionary[i], StructuralDof.TranslationY, amount: 0d));
				constraints.Add(new NodalDisplacement(model.NodesDictionary[i], StructuralDof.TranslationZ, amount: 0d));
			}

			var loads = new List<INodalLoadBoundaryCondition>();
			for (var i = 17; i < 21; i++)
			{
				loads.Add(new NodalLoad(model.NodesDictionary[i], StructuralDof.TranslationX, amount: 1 * 850d));
			}

			model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

			var computedDisplacements = SolveModel(model);

			return new double[] { computedDisplacements.GetTotalDisplacement(1, computedDisplacements.WatchDofs[4].Item1, StructuralDof.TranslationZ) };
		}

		private IncrementalDisplacementsLog SolveModel(Model model)
		{
			var solverFactory = new SkylineSolver.Factory();
			var algebraicModel = solverFactory.BuildAlgebraicModel(model);
			var solver = solverFactory.BuildSolver(algebraicModel);
			var problem = new ProblemStructural(model, algebraicModel);

			var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 2)
			{
				ResidualTolerance = 1E-8,
				MaxIterationsPerIncrement = 100,
				NumIterationsForMatrixRebuild = 1
			};
			var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();
			var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);

			loadControlAnalyzer.IncrementalDisplacementsLog = new IncrementalDisplacementsLog(
				new List<(INode node, IDofType dof)>()
				{
					(model.NodesDictionary[5], StructuralDof.TranslationX),
					(model.NodesDictionary[8], StructuralDof.TranslationZ),
					(model.NodesDictionary[12], StructuralDof.TranslationZ),
					(model.NodesDictionary[16], StructuralDof.TranslationZ),
					(model.NodesDictionary[20], StructuralDof.TranslationZ)
				}, algebraicModel
			);

			staticAnalyzer.Initialize();
			staticAnalyzer.Solve();

			return loadControlAnalyzer.IncrementalDisplacementsLog;
		}
	}
}
