using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml.Linq;

using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.MachineLearning;
using MGroup.Constitutive.Structural.Transient;
using MGroup.FEM.Structural.Continuum;
using MGroup.FEM.Structural.Line;
using MGroup.MachineLearning.TensorFlow.NeuralNetworks;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Multiscale.SupportiveClasses;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.Solvers.Direct;
using MGroup.CntConcreteBayesian.Tests.Commons;
using MiMsolve.multiScaleSupportiveClasses;

namespace MGroup.CntConcreteBayesian.Tests.Models
{

    public class ConcreteHierarchicalModels
    {

        private int[][,] cementElements;
        private double[][,] cementNodes;
        private int[][] cementSets;
        private int[][] cementRenumbering;
        private FeedForwardNeuralNetwork cementNeuralNetwork;

        private int[][,] mortarElements;
        private double[][,] mortarNodes;
        private int[][] mortarSets;
        private int[][] mortarRenumbering;
        private FeedForwardNeuralNetwork mortarNeuralNetwork;

        private int[][,] concreteElements;
        private double[][,] concreteNodes;
        private int[][] concreteSets;
        private int[][] concreteRenumbering;
        private FeedForwardNeuralNetwork concreteNeuralNetwork;

        public void InitializeModels()
        {
            //import model information from abaqus
            var SpecPath = @"MsolveOutputs\Inputs";
            var InputFileName = "CementCnt.inp";
            var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
            var pathName = Path.Combine(BasePath, SpecPath);

            string InputExtension = Path.GetExtension(InputFileName);
            string InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(InputFileName));
            string inputFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);
            string[] setNames = new string[] { "Constraints", "Loads", "Monitor" };
            (cementNodes, cementElements, cementSets) = AbaqusReader.ReadFile(inputFile, setNames);

            string initialPath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            var folderName = "DataFiles";
            var netPathName = "cement_network_architecture";
            netPathName = Path.Combine(initialPath, folderName, netPathName);
            var weightsPathName = "cement_trained_weights";
            weightsPathName = Path.Combine(initialPath, folderName, weightsPathName);
            var normalizationPathName = "cement_normalization";
            normalizationPathName = Path.Combine(initialPath, folderName, normalizationPathName);

            cementNeuralNetwork = new FeedForwardNeuralNetwork();
            cementNeuralNetwork.LoadNetwork(netPathName, weightsPathName, normalizationPathName);

            SpecPath = @"MsolveOutputs\Inputs";
            InputFileName = "MortarTensile.inp";
            BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
            pathName = Path.Combine(BasePath, SpecPath);

            InputExtension = Path.GetExtension(InputFileName);
            InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(InputFileName));
            inputFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);
            setNames = new string[] { "Constraints", "Loads", "Monitor" };
            (mortarNodes, mortarElements, mortarSets) = AbaqusReader.ReadFile(inputFile, setNames);

            initialPath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            folderName = "DataFiles";
            netPathName = "mortar_network_architecture";
            netPathName = Path.Combine(initialPath, folderName, netPathName);
            weightsPathName = "mortar_trained_weights";
            weightsPathName = Path.Combine(initialPath, folderName, weightsPathName);
            normalizationPathName = "mortar_normalization";
            normalizationPathName = Path.Combine(initialPath, folderName, normalizationPathName);

            mortarNeuralNetwork = new FeedForwardNeuralNetwork();
            mortarNeuralNetwork.LoadNetwork(netPathName, weightsPathName, normalizationPathName);

            SpecPath = @"MsolveOutputs\Inputs";
            InputFileName = "ReinfConcrBeam.inp";
            BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
            pathName = Path.Combine(BasePath, SpecPath);

            InputExtension = Path.GetExtension(InputFileName);
            InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(InputFileName));
            inputFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);
            setNames = new string[] { "Constraints", "Loads", "Monitor" };
            (concreteNodes, concreteElements, concreteSets) = AbaqusReader.ReadFile(inputFile, setNames);

            initialPath = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location).Split(new string[] { "\\bin" }, StringSplitOptions.None)[0];
            folderName = "DataFiles";
            netPathName = "concrete_network_architecture";
            netPathName = Path.Combine(initialPath, folderName, netPathName);
            weightsPathName = "concrete_trained_weights";
            weightsPathName = Path.Combine(initialPath, folderName, weightsPathName);
            normalizationPathName = "concrete_normalization";
            normalizationPathName = Path.Combine(initialPath, folderName, normalizationPathName);

            concreteNeuralNetwork = new FeedForwardNeuralNetwork();
            concreteNeuralNetwork.LoadNetwork(netPathName, weightsPathName, normalizationPathName);
        }

        private Model CreateCementScaleModel(double[] parameters)
        {
            var model = new Model();

            model.SubdomainsDictionary.Add(key: 0, new Subdomain(id: 0));

            for (var i = 0; i < cementNodes[0].GetLength(0); i++)
            {
                var nodeId = i + 1;
                model.NodesDictionary.Add(nodeId, new Node(
                    id: nodeId,
                    x: 0.001 * cementNodes[0][i, 1],
                    y: 0.001 * cementNodes[0][i, 2],
                    z: 0.001 * cementNodes[0][i, 3]));
            }

            for (var i = 0; i < cementElements[0].GetLength(0); i++)
            {
                cementRenumbering = new int[1][];
                cementRenumbering[0] = new int[] { 4, 8, 7, 3, 1, 5, 6, 2 };
                var nodeSet = new Node[cementRenumbering[0].Length];
                for (var j = 0; j < cementRenumbering[0].Length; j++)
                {
                    var nodeID = cementElements[0][i, cementRenumbering[0][j]];
                    nodeSet[j] = (Node)model.NodesDictionary[nodeID];
                }

                //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new ElasticMaterial3D(20, 0.2), new TransientAnalysisProperties(1, 0, 0));
                //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new DruckerPrager3DFunctional(20, 0.2, 30, 30, 0.025, x => 0.025 + 2 * x), new TransientAnalysisProperties(1, 0, 0));
                var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new NeuralNetworkMaterial3D(cementNeuralNetwork, new double[3] { parameters[0], parameters[1], parameters[2] }), new TransientAnalysisProperties(1, 0, 0));
                var element = elementFactory.CreateElement(CellType.Hexa8, nodeSet);
                element.ID = i + 1;

                model.ElementsDictionary.Add(element.ID, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            var constraints = new List<INodalDisplacementBoundaryCondition>();
            for (var i = 0; i < cementSets[0].GetLength(0); i++)
            {
                constraints.Add(new NodalDisplacement(model.NodesDictionary[cementSets[0][i]], StructuralDof.TranslationX, amount: 0d));
                constraints.Add(new NodalDisplacement(model.NodesDictionary[cementSets[0][i]], StructuralDof.TranslationY, amount: 0d));
                constraints.Add(new NodalDisplacement(model.NodesDictionary[cementSets[0][i]], StructuralDof.TranslationZ, amount: 0d));
            }

            var loads = new List<INodalLoadBoundaryCondition>();
            for (var i = 0; i < cementSets[1].GetLength(0); i++)
            {
                loads.Add(new NodalLoad(model.NodesDictionary[cementSets[1][i]], StructuralDof.TranslationY, amount: 3 * 1e-6d));
            }

            model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

            return model;
        }

        //private Model CreateMortarScaleModel(double[] parameters)
        //{
        //    var model = new Model();

        //    model.SubdomainsDictionary.Add(key: 0, new Subdomain(id: 0));

        //    for (var i = 0; i < mortarNodes[0].GetLength(0); i++)
        //    {
        //        var nodeId = i + 1;
        //        model.NodesDictionary.Add(nodeId, new Node(
        //            id: nodeId,
        //            x: 0.001 * mortarNodes[0][i, 1],
        //            y: 0.001 * mortarNodes[0][i, 2],
        //            z: 0.001 * mortarNodes[0][i, 3]));
        //    }

        //    for (var i = 0; i < mortarElements[0].GetLength(0); i++)
        //    {
        //        mortarRenumbering = new int[1][];
        //        mortarRenumbering[0] = new int[] { 1, 2, 3, 4 };
        //        var nodeSet = new Node[mortarRenumbering[0].Length];
        //        for (var j = 0; j < mortarRenumbering[0].Length; j++)
        //        {
        //            var nodeID = mortarElements[0][i, mortarRenumbering[0][j]];
        //            nodeSet[j] = (Node)model.NodesDictionary[nodeID];
        //        }

        //        //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new ElasticMaterial3D(parameters[0], 0.2), new TransientAnalysisProperties(1, 0, 0));
        //        //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new DruckerPrager3DFunctional(20, 0.2, 30, 30, 0.025, x => 0.025 + 2 * x), new TransientAnalysisProperties(1, 0, 0));
        //        var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new NeuralNetworkMaterial3D(mortarNeuralNetwork, new double[3] { parameters[0], parameters[1], parameters[2] }), new TransientAnalysisProperties(1, 0, 0));
        //        var element = elementFactory.CreateElement(CellType.Tet4, nodeSet);
        //        element.ID = i + 1;

        //        model.ElementsDictionary.Add(element.ID, element);
        //        model.SubdomainsDictionary[0].Elements.Add(element);
        //    }

        //    var constraints = new List<INodalDisplacementBoundaryCondition>();
        //    for (var i = 0; i < mortarSets[0].GetLength(0); i++)
        //    {
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[mortarSets[0][i]], StructuralDof.TranslationX, amount: 0d));
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[mortarSets[0][i]], StructuralDof.TranslationY, amount: 0d));
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[mortarSets[0][i]], StructuralDof.TranslationZ, amount: 0d));
        //    }

        //    var loads = new List<INodalLoadBoundaryCondition>();
        //    for (var i = 0; i < mortarSets[1].GetLength(0); i++)
        //    {
        //        loads.Add(new NodalLoad(model.NodesDictionary[mortarSets[1][i]], StructuralDof.TranslationZ, amount: 4 * 13.65e-6d));
        //    }

        //    model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

        //    return model;
        //}

        //private Model CreateConcreteScaleModel(double[] parameters)
        //{
        //    var model = new Model();

        //    model.SubdomainsDictionary.Add(key: 0, new Subdomain(id: 0));

        //    for (var i = 0; i < concreteNodes[0].GetLength(0); i++)
        //    {
        //        var nodeId = i + 1;
        //        model.NodesDictionary.Add(nodeId, new Node(
        //            id: nodeId,
        //            x: 0.001 * concreteNodes[0][i, 1],
        //            y: 0.001 * concreteNodes[0][i, 2],
        //            z: 0.001 * concreteNodes[0][i, 3]));
        //    }

        //    for (var i = 0; i < concreteNodes[1].GetLength(0); i++)
        //    {
        //        var nodeId = concreteNodes[0].GetLength(0) + i + 1;
        //        model.NodesDictionary.Add(nodeId, new Node(
        //            id: nodeId,
        //            x: 0.001 * concreteNodes[1][i, 1],
        //            y: 0.001 * concreteNodes[1][i, 2],
        //            z: 0.001 * concreteNodes[1][i, 3]));
        //    }

        //    for (var i = 0; i < concreteElements[0].GetLength(0); i++)
        //    {
        //        concreteRenumbering = new int[2][];
        //        concreteRenumbering[0] = new int[] { 4, 8, 7, 3, 1, 5, 6, 2 };
        //        var nodeSet = new Node[concreteRenumbering[0].Length];
        //        for (var j = 0; j < concreteRenumbering[0].Length; j++)
        //        {
        //            var nodeID = concreteElements[0][i, concreteRenumbering[0][j]];
        //            nodeSet[j] = (Node)model.NodesDictionary[nodeID];
        //        }

        //        //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new ElasticMaterial3D(20, 0.2), new TransientAnalysisProperties(1, 0, 0));
        //        //var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new DruckerPrager3DFunctional(20, 0.2, 30, 30, 0.025, x => 0.025 + 2 * x), new TransientAnalysisProperties(1, 0, 0));
        //        var elementFactory = new FEM.Structural.Continuum.ContinuumElement3DFactory(new NeuralNetworkMaterial3D(concreteNeuralNetwork, new double[3] { parameters[0], parameters[1], parameters[2] }), new TransientAnalysisProperties(1, 0, 0));
        //        var element = elementFactory.CreateElement(CellType.Hexa8, nodeSet);
        //        element.ID = i + 1;

        //        model.ElementsDictionary.Add(element.ID, element);
        //        model.SubdomainsDictionary[0].Elements.Add(element);
        //    }

        //    for (var i = 0; i < concreteElements[1].GetLength(0); i++)
        //    {
        //        concreteRenumbering[1] = new int[] { 1, 2 };
        //        var rebarNodeSet = new Node[concreteRenumbering[1].Length];
        //        for (var j = 0; j < concreteRenumbering[1].Length; j++)
        //        {
        //            var nodeID = concreteNodes[0].GetLength(0) + concreteElements[1][i, concreteRenumbering[1][j]];
        //            rebarNodeSet[j] = (Node)model.NodesDictionary[nodeID];
        //        }

        //        var rebarElement = new Rod3D(new List<INode>() { rebarNodeSet[0], rebarNodeSet[1] }, youngModulus: 200) { Density = 7.750, SectionArea = 0.00031415926 };
        //        rebarElement.ID = concreteElements[0].GetLength(0) + i + 1;

        //        model.ElementsDictionary.Add(rebarElement.ID, rebarElement);
        //        model.SubdomainsDictionary[0].Elements.Add(rebarElement);
        //    }

        //    var hostElements = concreteElements[0].GetLength(0);
        //    var embeddedElements = concreteElements[1].GetLength(0);

        //    var embeddedGrouping = EmbeddedBeam3DGrouping.CreateFullyBonded(model, model.ElementsDictionary
        //    .Where(x => x.Key <= hostElements).Select(kv => kv.Value).ToArray(), model.ElementsDictionary.Where(x => x.Key > hostElements)
        //    .Select(kv => kv.Value).ToArray(), false);

        //    var constraints = new List<INodalDisplacementBoundaryCondition>();
        //    for (var i = 0; i < concreteSets[0].GetLength(0); i++)
        //    {
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[concreteSets[0][i]], StructuralDof.TranslationX, amount: 0d));
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[concreteSets[0][i]], StructuralDof.TranslationY, amount: 0d));
        //        constraints.Add(new NodalDisplacement(model.NodesDictionary[concreteSets[0][i]], StructuralDof.TranslationZ, amount: 0d));
        //    }

        //    var loads = new List<INodalLoadBoundaryCondition>();
        //    for (var i = 0; i < concreteSets[1].GetLength(0); i++)
        //    {
        //        loads.Add(new NodalLoad(model.NodesDictionary[concreteSets[1][i]], StructuralDof.TranslationY, amount: 4 * 150 * 1e-6d));
        //    }

        //    model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

        //    return model;
        //}

        public double[] FormulateProblem(double[] parameters)
        {
            //create models in MSolve
            var cementModel = CreateCementScaleModel(parameters);

            //var mortarModel = CreateMortarScaleModel(parameters);

            //var concreteModel = CreateConcreteScaleModel(parameters);

            var cementModelOutput = SolveCementModel(cementModel);

            //var mortarModelOutput = SolveMortarModel(mortarModel);

            //var concreteModelOutput = SolveConcreteModel(concreteModel);

            return new double[]
            {
                cementModelOutput.GetTotalDisplacement(1, cementModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY),
                //mortarModelOutput.GetTotalDisplacement(1, mortarModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationZ),
                //concreteModelOutput.GetTotalDisplacement(1, concreteModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY),
            };
        }

        public double[] FormulateCementProblem(double[] parameters)
        {
            //create models in MSolve
            var cementModel = CreateCementScaleModel(parameters);

            var cementModelOutput = SolveCementModel(cementModel);
            var temp = cementModelOutput.GetTotalDisplacement(4, cementModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY);
            return new double[]
            {
                cementModelOutput.GetTotalDisplacement(4, cementModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY),
            };
        }

        //public double[] FormulateMortarProblem(double[] parameters)
        //{
        //    //create models in MSolve
        //    var mortarModel = CreateMortarScaleModel(parameters);

        //    var mortarModelOutput = SolveMortarModel(mortarModel);
        //    var temp = mortarModelOutput.GetTotalDisplacement(4, mortarModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationZ);
        //    return new double[]
        //    {
        //        mortarModelOutput.GetTotalDisplacement(4, mortarModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationZ),
        //    };
        //}

        //public double[] FormulateConcreteProblem(double[] parameters)
        //{
        //    //create models in MSolve
        //    var concreteModel = CreateConcreteScaleModel(parameters);

        //    var concreteModelOutput = SolveConcreteModel(concreteModel);
        //    var temp = concreteModelOutput.GetTotalDisplacement(4, concreteModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY);
        //    return new double[]
        //    {
        //        concreteModelOutput.GetTotalDisplacement(4, concreteModelOutput.WatchDofs[0].Item1, StructuralDof.TranslationY),
        //    };
        //}

        private TotalDisplacementsPerIncrementLog SolveCementModel(Model model)
        {
            var solverFactory = new SuiteSparseSolver.Factory();
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var problem = new ProblemStructural(model, algebraicModel);

            var loadControlAnalyzerBuilder = new MiMsolve.multiScaleSupportiveClasses.LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 5)
            {
                ResidualTolerance = 1E-2,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();
            var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);

            var monitorList = new List<(INode node, IDofType dof)>();
            for (int i = 0; i < cementSets[2].GetLength(0); i++)
            {
                monitorList.Add((model.NodesDictionary[cementSets[2][i]], StructuralDof.TranslationY));
            }
            loadControlAnalyzer.TotalDisplacementsPerIncrementLog = new TotalDisplacementsPerIncrementLog(
                monitorList, algebraicModel
            );

            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            return loadControlAnalyzer.TotalDisplacementsPerIncrementLog;
        }

        //private TotalDisplacementsPerIncrementLog SolveMortarModel(Model model)
        //{
        //    var solverFactory = new SuiteSparseSolver.Factory();
        //    var algebraicModel = solverFactory.BuildAlgebraicModel(model);
        //    var solver = solverFactory.BuildSolver(algebraicModel);
        //    var problem = new ProblemStructural(model, algebraicModel);

        //    var loadControlAnalyzerBuilder = new MiMsolve.multiScaleSupportiveClasses.LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 5)
        //    {
        //        ResidualTolerance = 1E-2,
        //        MaxIterationsPerIncrement = 100,
        //        NumIterationsForMatrixRebuild = 1
        //    };
        //    var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();
        //    var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);

        //    var monitorList = new List<(INode node, IDofType dof)>();
        //    for (int i = 0; i < mortarSets[2].GetLength(0); i++)
        //    {
        //        monitorList.Add((model.NodesDictionary[mortarSets[2][i]], StructuralDof.TranslationZ));
        //    }
        //    loadControlAnalyzer.TotalDisplacementsPerIncrementLog = new TotalDisplacementsPerIncrementLog(
        //        monitorList, algebraicModel
        //    );

        //    staticAnalyzer.Initialize();
        //    staticAnalyzer.Solve();

        //    return loadControlAnalyzer.TotalDisplacementsPerIncrementLog;
        //}

        //private TotalDisplacementsPerIncrementLog SolveConcreteModel(Model model)
        //{
        //    var solverFactory = new SuiteSparseSolver.Factory();
        //    var algebraicModel = solverFactory.BuildAlgebraicModel(model);
        //    var solver = solverFactory.BuildSolver(algebraicModel);
        //    var problem = new ProblemStructural(model, algebraicModel);

        //    var loadControlAnalyzerBuilder = new MiMsolve.multiScaleSupportiveClasses.LoadControlAnalyzer.Builder(algebraicModel, solver, problem, numIncrements: 5)
        //    {
        //        ResidualTolerance = 1E-2,
        //        MaxIterationsPerIncrement = 100,
        //        NumIterationsForMatrixRebuild = 1
        //    };
        //    var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();
        //    var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, loadControlAnalyzer);

        //    var monitorList = new List<(INode node, IDofType dof)>();
        //    for (int i = 0; i < concreteSets[2].GetLength(0); i++)
        //    {
        //        monitorList.Add((model.NodesDictionary[concreteSets[2][i]], StructuralDof.TranslationY));
        //    }
        //    loadControlAnalyzer.TotalDisplacementsPerIncrementLog = new TotalDisplacementsPerIncrementLog(
        //        monitorList, algebraicModel
        //    );

        //    staticAnalyzer.Initialize();
        //    staticAnalyzer.Solve();

        //    return loadControlAnalyzer.TotalDisplacementsPerIncrementLog;
        //}
    }
}
