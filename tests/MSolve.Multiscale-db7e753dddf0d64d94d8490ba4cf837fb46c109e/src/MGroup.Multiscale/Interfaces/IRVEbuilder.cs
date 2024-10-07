using MGroup.MSolve.Discretization.Entities;
using System;
using System.Collections.Generic;


namespace MGroup.MSolve.MultiscaleAnalysis.Interfaces
{
    /// <summary>
    /// Indicates the nesessary methods that should be implemented by builders of models intended to be used as rves in Multiscale Problems
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public interface IRVEbuilder
    {
        Tuple<Model, Dictionary<int, INode>,double> GetModelAndBoundaryNodes();
        IRVEbuilder Clone(int a);
    }
}
