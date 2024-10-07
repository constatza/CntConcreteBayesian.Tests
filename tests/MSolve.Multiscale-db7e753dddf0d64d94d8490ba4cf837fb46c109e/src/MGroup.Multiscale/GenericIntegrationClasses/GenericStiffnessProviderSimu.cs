using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;
using MGroup.MSolve.MultiscaleAnalysis;

namespace MGroup.MSolve.MultiscaleAnalysis
{
    /// <summary>
    /// Element stiffness matrix provider for simultaneous calculation of global stiffness matrix macroscopic variables in multiscale FE2 scheme
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class GenericStiffnessProviderSimu<TMatrix> : IElementMatrixProvider
       where TMatrix : class, IMatrix
    {
        #region IElementMatrixProvider Members
        private GenericSubdomainCalculationsAndAssembly<TMatrix> host;

        public GenericStiffnessProviderSimu(GenericSubdomainCalculationsAndAssembly<TMatrix> host)
        {
            this.host = host;
        }

        public IMatrix Matrix(IElementType element) 
        {
            var elementMatrix = element.PhysicsMatrix();
            host.UpdateVectors(element, elementMatrix);
            return elementMatrix;
        }

        #endregion
    }
}



