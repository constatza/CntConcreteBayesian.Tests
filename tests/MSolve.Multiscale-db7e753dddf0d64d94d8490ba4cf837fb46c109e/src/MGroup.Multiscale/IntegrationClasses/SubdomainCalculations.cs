using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using System.Linq;
using MGroup.Solvers.AlgebraicModel;
using MGroup.Constitutive.Structural;
using MGroup.Solvers.Results;

namespace MGroup.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class SubdomainCalculations
    {
        
        public static Dictionary<int,int> GetNodesOrderInDictionary(Dictionary<int, INode> boundaryNodes)
        {
            Dictionary<int, int> boundaryNodesOrder = new Dictionary<int, int>();
            int order = 1;
            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                boundaryNodesOrder.Add(boundaryNode.ID, order);
                order += 1;
            }
            return boundaryNodesOrder;

        }

        
        #region v2 methods
        public static double[][] SubtractConsecutiveVectors(double[][] KppDqVectors, double[][] f3_vectors)
        {
            double[][] f4_vectors = new double[KppDqVectors.GetLength(0)][];
            for (int j1 = 0; j1 < f4_vectors.GetLength(0); j1++)
            {
                f4_vectors[j1] = new double[KppDqVectors[0].GetLength(0)];
                for (int j2 = 0; j2 < f4_vectors[0].GetLength(0); j2++)
                {
                    f4_vectors[j1][j2] = KppDqVectors[j1][j2] - f3_vectors[j1][j2];
                }
            }

            return f4_vectors;

        }

        public static double[,] CalculateDqCondDq(double[][] f4_vectors, IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes)
        {
            double[,] DqCondDq = new double[scaleTransitions.MacroscaleVariableDimension(), scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node boundaryNode in boundaryNodes.Values)
            {
                for (int i1 = 0; i1 < f4_vectors.GetLength(0); i1++)
                {
                    double[] f4DataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                    for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                    {
                        f4DataTriplette[i2] = f4_vectors[i1][scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                    }
                    double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, f4DataTriplette);
                    for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                    {
                        DqCondDq[i3, i1] += contribution[i3];
                    }
                }

            }
            return DqCondDq;
        }

        public static double[] CalculateDqFpp(double[] FppReactionVector, IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes)
        {
            double[] DqFpp = new double[scaleTransitions.MacroscaleVariableDimension()];

            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);

            foreach (Node boundaryNode in boundaryNodes.Values)
            {

                double[] FppDataTriplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                for (int i2 = 0; i2 < scaleTransitions.PrescribedDofsPerNode(); i2++)
                {
                    FppDataTriplette[i2] = FppReactionVector[scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[boundaryNode.ID] - 1) + i2];
                }
                double[] contribution = scaleTransitions.MicroToMacroTransition(boundaryNode, FppDataTriplette);
                for (int i3 = 0; i3 < scaleTransitions.MacroscaleVariableDimension(); i3++)
                {
                    DqFpp[i3] += contribution[i3];
                }


            }

            return DqFpp;

        }

        public static double[] CalculateFppReactionsVector(Subdomain subdomain, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes, List<int> boundaryNodesIds, IGlobalVector solution,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements, IAlgebraicModel algebraicModel )
        {
            
            double[] FppReactionVector;
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);
            FppReactionVector = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            foreach (IElementType element in subdomain.EnumerateElements())
            {
                var isEmbeddedElement = element is Discretization.Embedding.IEmbeddedElement;
                var elStart = DateTime.Now;

                var localSolution = algebraicModel.ExtractElementVector(solution, element); // subdomain.GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
                ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodesIds, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                //double[] localdSolution = subdomain.GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
                double[] f = element.CalculateResponseIntegral(); //element.CalculateForces(element, localSolution, localdSolution);

                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodesIds.Contains(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            FppReactionVector[dofrow_p] += f[iElementMatrixRow + i1];

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return FppReactionVector;

        }

        public static void ImposePrescribedDisplacementsWithInitialConditionSEffect(IElementType element, double[] localSolution, List<int> boundaryNodesIds,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];
                int nodalDofsNumber = elementDOFTypes[j].Count;
                if (boundaryNodesIds.Contains(nodeColumn.ID))
                {
                    Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    int positionOfDofInNode = 0;
                    foreach (IDofType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDofInNode += 1;
                    }
                }
                iElementMatrixColumn += nodalDofsNumber;
            }

        }

        public static double[][] CalculateKpfKffinverseKfpDq(IGlobalVector[] f2_vectors, Subdomain subdomain, IElementMatrixProvider elementProvider,
            IScaleTransitions scaleTransitions, Dictionary<int, INode> boundaryNodes, GlobalAlgebraicModel<Matrix> globalAlgebraicModel)
        {
            var ActiveDofs = new ActiveDofs();
            ActiveDofs.AddDof(StructuralDof.TranslationX);
            ActiveDofs.AddDof(StructuralDof.TranslationY);
            ActiveDofs.AddDof(StructuralDof.TranslationZ);
            ActiveDofs.AddDof(StructuralDof.RotationX);
            ActiveDofs.AddDof(StructuralDof.RotationY);
            ActiveDofs.AddDof(StructuralDof.RotationZ);

            var freeDofs = globalAlgebraicModel.SubdomainFreeDofOrdering.FreeDofs;
            //var dofOrdering = subdomain.FreeDofOrdering; //.1
            //var FreeDofs = subdomain.FreeDofOrdering.FreeDofs;//.1Dictionary<int, Dictionary<DOFType, int>> nodalDOFsDictionary = subdomain.NodalDOFsDictionary;

            double[][] f3_vectors = new double[f2_vectors.GetLength(0)][];
            for (int i1 = 0; i1 < f2_vectors.GetLength(0); i1++)
            {
                f3_vectors[i1] = new double[scaleTransitions.PrescribedDofsPerNode() * boundaryNodes.Count];
            }
            Dictionary<int, int> boundaryNodesOrder = GetNodesOrderInDictionary(boundaryNodes);


            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);

            NodalResults[] f2NodalResults = new NodalResults[f2_vectors.Length];
            for (int j1 = 0; j1 < f2_vectors.Length; j1++) { f2NodalResults[j1]= globalAlgebraicModel.ExtractAllResults(f2_vectors[j1]); }
            

            foreach (IElementType element in subdomain.Elements)
            {
                var isEmbeddedElement = element is Discretization.Embedding.IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix ElementK = elementProvider.Matrix(element);
                times["element"] += DateTime.Now - elStart;

                elStart = DateTime.Now;
                var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    if (boundaryNodes.ContainsKey(nodeRow.ID))
                    {
                        for (int i1 = 0; i1 < scaleTransitions.PrescribedDofsPerNode(); i1++)
                        {
                            int dofrow_p = scaleTransitions.PrescribedDofsPerNode() * (boundaryNodesOrder[nodeRow.ID] - 1) + i1;
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                int dofTypeColumnToNumber = -1;
                                foreach (IDofType dofTypeColumn in elementDOFTypes[j])
                                {
                                    dofTypeColumnToNumber++;
                                    int dofIdINode = ActiveDofs.GetIdOfDof(elementDOFTypes[i][dofTypeColumnToNumber]);
                                    bool isFree = freeDofs.TryGetValue(matrixAssemblyNodes[j].ID, dofIdINode, out int dofColumn);
                                    // v2.4 int dofColumn = nodalDOFsDictionary.ContainsKey(nodeColumn.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeColumn.ID][dofTypeColumn];
                                    if (isFree)// TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                                    {                   // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok

                                        for (int i2 = 0; i2 < f2_vectors.GetLength(0); i2++)
                                        {
                                            f3_vectors[i2][dofrow_p] += ElementK[iElementMatrixRow + i1, iElementMatrixColumn] * f2NodalResults[i2].Data[matrixAssemblyNodes[j].ID, dofIdINode];
                                        }

                                    }
                                    iElementMatrixColumn++;
                                }
                            }

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                times["addition"] += DateTime.Now - elStart;
            }
            var totalTime = DateTime.Now - totalStart;

            return f3_vectors;

        }
        #endregion
    }
}
