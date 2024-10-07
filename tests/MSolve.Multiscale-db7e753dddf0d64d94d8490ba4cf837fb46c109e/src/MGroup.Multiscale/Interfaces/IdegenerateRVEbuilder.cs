using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using System.Collections.Generic;


namespace MGroup.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates additional methods that should be implemented by rveBuilders that will be used in FE2 3D to 2D degenerate analysis
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IdegenerateRVEbuilder : IRVEbuilder
    {
        Dictionary<Node, IList<IStructuralDofType>> GetModelRigidBodyNodeConstraints(Model model);
    }
}
