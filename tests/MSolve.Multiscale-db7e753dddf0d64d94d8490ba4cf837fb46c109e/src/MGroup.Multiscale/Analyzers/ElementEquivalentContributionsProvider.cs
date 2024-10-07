using System;
using System.Collections.Generic;
using System.Collections.Generic;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using System.Linq;
using MGroup.Solvers.AlgebraicModel;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DofOrdering;
using MGroup.Solvers;

namespace MiMsolve.intermediateCodeDevelopmentClasses
{
    /// <summary>
    /// Assembles the equivalent contribution of imposed displacements in a bvp with non zero initial conditions.
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class ElementEquivalentContributionsProvider : IElementVectorProvider
    {
        private readonly IAlgebraicModel algebraicModel;
        private IElementMatrixProvider elementProvider;
        private readonly ISubdomain subdomain;
        private readonly IModel model;
        private ActiveDofs ActiveDofs;
        private ISubdomainFreeDofOrdering subdomainFreeDofORdering;
        private IntDofTable freeDofs;
        private Dictionary<int, INode> boundaryNodes;
        private Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements;
        private Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements;
        private int nIncrement;
        private int totalIncrements;
        /// <summary>
        /// ELEMENT provider tha perastei profanws o ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
        /// kai subdomain prosoxh sta ID idia me ta linearsystems
        /// </summary>
        /// <param name="subdomain"></param>
        /// <param name="elementProvider"></param>


        public ElementEquivalentContributionsProvider(IAlgebraicModel algebraicModel, IElementMatrixProvider elementProvider, Dictionary<int, INode> boundaryNodes,
           Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
           int nIncrement, int totalIncrements)
        {
            this.algebraicModel = algebraicModel;
            this.elementProvider = elementProvider;
            ActiveDofs = new ActiveDofs();
            ActiveDofs.AddDof(StructuralDof.TranslationX);
            ActiveDofs.AddDof(StructuralDof.TranslationY);
            ActiveDofs.AddDof(StructuralDof.TranslationZ);
            ActiveDofs.AddDof(StructuralDof.RotationX);
            ActiveDofs.AddDof(StructuralDof.RotationY);
            ActiveDofs.AddDof(StructuralDof.RotationZ);
            this.boundaryNodes = boundaryNodes;
            this.initialConvergedBoundaryDisplacements = initialConvergedBoundaryDisplacements;
            this.totalBoundaryDisplacements = totalBoundaryDisplacements;
            this.nIncrement = nIncrement;
            this.totalIncrements = totalIncrements;

        }

        public double[] CalcVector(IElementType element)
        {
            var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
            var totalElementDofs = elementDOFTypes.Select(x => x.Count()).ToList().Sum();         //TODO if element exei bondary nodes

            double[] uStep_values_orZero_for_free = new double[totalElementDofs];
            bool boundaryNodeFound = false;
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];


                int nodalDofsNumber = elementDOFTypes[j].Count; //TODOGerasimos elegxos oti edw oi ginetai prosvash apo 0:1:megethos 
                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                {
                    boundaryNodeFound = true;


                    double[] element_Kfp_triplette = new double[nodalDofsNumber]; //nodalDofsNumber: giati oxi scaleTransitions.PrescribedDofsPerNode()? Dioti tou ta pairname ola(triplette) kai dialegei to 


                    Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    //double[] uStep_values_orZero_for_free = new double[nodalDofsNumber];

                    int positionOfDof = 0;
                    foreach (IDofType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            uStep_values_orZero_for_free[iElementMatrixColumn + positionOfDof] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * (nIncrement / (double)totalIncrements);
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDof += 1;
                    }





                }
                iElementMatrixColumn += nodalDofsNumber;

            }

            if (boundaryNodeFound)
            {
                IMatrix ElementK = elementProvider.Matrix(element);
                Vector uStep = Vector.CreateFromArray(uStep_values_orZero_for_free);
                var elementContribution = ElementK.Multiply(uStep, false);
                var elementContributionCheck = new double[totalElementDofs];
                for (int j1 = 0; j1 < totalElementDofs; j1++)
                {
                    for (int j2 = 0; j2 < totalElementDofs; j2++)
                    {
                        elementContributionCheck[j1] += ElementK[j1, j2] * uStep[j2];
                    }
                }
                return elementContribution.CopyToArray();
            }
            else
            {
                return new double[totalElementDofs];
            }
        }


    }
}
