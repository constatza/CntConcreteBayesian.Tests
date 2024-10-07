using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Solvers.AlgebraicModel;
using MGroup.Solvers.Assemblers;
using MGroup.Solvers.DofOrdering;
using MGroup.Solvers.LinearSystem;
using System.Collections.Generic;
using System.Linq;

namespace MGroup.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Supportive class  that implements nesessary integration methods associated with FE2 multiscale analysis 
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class GenericSubdomainCalculationsAndAssembly<TMatrix> 
        where TMatrix : class, IMatrix
    {
        //SubdomainCalculationsSimultaneousObje

        private GlobalVector[] KfpDqVectors;
        private double[][] KppDqVectors;
        IScaleTransitions scaleTransitions;
        private SubdomainVectorAssembler vectorAssembler;
        private ISubdomainFreeDofOrdering subdomainFreeDofOrdering;
        Dictionary<int, IElementType> boundaryElements;
        Dictionary<int, INode> boundaryNodes;
        Dictionary<int, int> boundaryNodesOrder;
        private ActiveDofs ActiveDofs;
        private GlobalAlgebraicModel<TMatrix> globalAlgebraicModel;
        private ISolver solver;

       


        public (IGlobalVector[], double[][]) UpdateSubdomainKffAndCalculateKfpDqAndKppDqpMultipleObje(Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions,
            Dictionary<int, INode> boundaryNodes, Dictionary<int, IElementType> boundaryElements, ISolver solver, GlobalAlgebraicModel<TMatrix> globalAlgebraicModel)
        {
            ActiveDofs = new ActiveDofs();
            ActiveDofs.AddDof(StructuralDof.TranslationX);
            ActiveDofs.AddDof(StructuralDof.TranslationY);
            ActiveDofs.AddDof(StructuralDof.TranslationZ);
            ActiveDofs.AddDof(StructuralDof.RotationX);
            ActiveDofs.AddDof(StructuralDof.RotationY);
            ActiveDofs.AddDof(StructuralDof.RotationZ);
            vectorAssembler = new SubdomainVectorAssembler(ActiveDofs);


            subdomainFreeDofOrdering = globalAlgebraicModel.SubdomainFreeDofOrdering;

            
			//two following commented out lines refer to multisubdomain model
            //Dictionary<int, double[][]> KfpDqSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            //Dictionary<int, double[][]> KppDqVectorsSubdomains = new Dictionary<int, double[][]>(model.SubdomainsDictionary.Count);
            this.boundaryElements = boundaryElements;
            this.boundaryNodes = boundaryNodes;
            this.scaleTransitions = scaleTransitions;
            this.globalAlgebraicModel= globalAlgebraicModel;
            this.solver = solver;


            #region Create KfpDq and KppDq vectors 
            KfpDqVectors = new GlobalVector[scaleTransitions.MacroscaleVariableDimension()];
            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
            {
                KfpDqVectors[j1] = globalAlgebraicModel.CreateZeroVector();  // new double[subdomainFreeDofOrdering.NumFreeDofs]; //v2.2 subdomain.TotalDOFs]; 
            }

            KppDqVectors = new double[scaleTransitions.MacroscaleVariableDimension()][];
            boundaryNodesOrder = SubdomainCalculations.GetNodesOrderInDictionary(boundaryNodes);
            for (int j1 = 0; j1 < scaleTransitions.MacroscaleVariableDimension(); j1++)
            {
                KppDqVectors[j1] = new double[boundaryNodesOrder.Count * scaleTransitions.PrescribedDofsPerNode()]; // h allliws subdomain.Forces.GetLength(0)
            }
            #endregion


            var StiffnessProvider = new GenericStiffnessProviderSimu<TMatrix>(this);
            //Dictionary<int, IMatrix> subdomainKs = solver.BuildGlobalMatrices(StiffnessProvider);
            var stiffness = globalAlgebraicModel.BuildGlobalMatrix(StiffnessProvider);
            solver.LinearSystem.Matrix = stiffness;




            return (KfpDqVectors, KppDqVectors);
        }

        

        public void UpdateVectors(IElementType element, IMatrix ElementK)
        {
            //ISubdomain subdomain = element.Subdomain;
            if (boundaryElements.ContainsKey(element.ID))
            {
                //ADDED these lines from another part of UpdateSubdomainKffAndCalculateKfpDqAndKppDqp
                var isEmbeddedElement = element is IEmbeddedElement;
                var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
                var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
                var elementKfpDqvecs = new double[scaleTransitions.MacroscaleVariableDimension()][];
                for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
                {
                    elementKfpDqvecs[i1] = new double[ElementK.NumRows];
                }

                #region KfpDq Multiplication
                int iElementMatrixRow = 0;
                for (int i = 0; i < elementDOFTypes.Count; i++)
                {
                    INode nodeRow = matrixAssemblyNodes[i];
                    int dofTypeRowToNumber = -1; //v2.6
                    foreach (IDofType dofTypeRow in elementDOFTypes[i])
                    {
                        dofTypeRowToNumber++;
                        int dofIdINode = ActiveDofs.GetIdOfDof(elementDOFTypes[i][dofTypeRowToNumber]);
                        bool isFree = subdomainFreeDofOrdering.FreeDofs.TryGetValue(matrixAssemblyNodes[i].ID, dofIdINode, out int dofRow); //v2.6
                        //int dofRow = nodalDOFsDictionary.ContainsKey(nodeRow.ID) == false && isEmbeddedElement ? -1 : nodalDOFsDictionary[nodeRow.ID][dofTypeRow];
                        if (isFree) // TODOGerasimos edw pithanws thelei kai elegxo alliws an den ta exoume afhsei constrained ta p kai einai elefthera px me to an anhkoun sto baoundary nodes
                        {                    // alla etsi einai oti akrivws thewritai kai sto assembly tou Kff opote ok
                            int iElementMatrixColumn = 0;
                            for (int j = 0; j < elementDOFTypes.Count; j++)
                            {
                                INode nodeColumn = matrixAssemblyNodes[j];
                                 int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
                                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                                {
                                    double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 
                                    for (int j1 = 0; j1 < nodalDofsNumber; j1++)                  //scaleTransitions.MicroToMacroTransition ti tha xrhsimopoihsei apo afta analoga pws einai implemented
                                    {
                                        element_Kfp_triplette[j1] = ElementK[iElementMatrixRow, iElementMatrixColumn + j1];
                                    }

                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kfp_triplette);
                                    for (int j2 = 0; j2 < contribution.GetLength(0); j2++)
                                    {
                                        
                                        elementKfpDqvecs[j2][iElementMatrixRow] += contribution[j2]; 
                                        
                                    }

                                }
                                iElementMatrixColumn += nodalDofsNumber;

                            }
                        }
                        iElementMatrixRow++;
                    }
                }


                var elementsEnum = new List<IElementType>() { element };
                var elementContribution = new ElementIntegrationContributionProvider(elementKfpDqvecs[0]);
                for (int j2 = 0; j2 < scaleTransitions.MacroscaleVariableDimension(); j2++)
                {
                    elementContribution.elementContribution=elementKfpDqvecs[j2];
                    vectorAssembler.AddToSubdomainVector(elementsEnum, KfpDqVectors[j2].SingleVector, elementContribution, subdomainFreeDofOrdering);

                }
                


                #endregion

                #region KppDq Multiplications
                iElementMatrixRow = 0;
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
                                
                                int nodalDofsNumber = elementDOFTypes[j].Count;
                                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                                {
                                    double[] element_Kpp_triplette = new double[scaleTransitions.PrescribedDofsPerNode()];
                                    for (int j2 = 0; j2 < scaleTransitions.PrescribedDofsPerNode(); j2++)
                                    {
                                        element_Kpp_triplette[j2] = ElementK[iElementMatrixRow + i1, iElementMatrixColumn + j2]; // mallon iElementMatrixRow + i1
                                    }
                                    double[] contribution = scaleTransitions.MicroToMacroTransition(nodeColumn, element_Kpp_triplette);
                                    for (int j1 = 0; j1 < contribution.GetLength(0); j1++)
                                    {
                                        KppDqVectors[j1][dofrow_p] += contribution[j1];
                                    }
                                }
                                iElementMatrixColumn += nodalDofsNumber;



                            }

                        }
                    }
                    iElementMatrixRow += elementDOFTypes[i].Count;
                }
                #endregion
            }
        }

    }

    public class ElementIntegrationContributionProvider : IElementVectorProvider
    {
        public double[] elementContribution { get; set; }
        public ElementIntegrationContributionProvider(double[] elementContribution)
        {
            this.elementContribution = elementContribution;
        }

        public double[] CalcVector(IElementType element) => elementContribution;
    }

}

