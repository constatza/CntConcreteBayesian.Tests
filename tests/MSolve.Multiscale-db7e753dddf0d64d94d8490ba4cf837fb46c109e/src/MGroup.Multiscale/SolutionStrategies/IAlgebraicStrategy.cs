using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.Solvers.AlgebraicModel;

namespace MiMsolve.SolutionStrategies
{
    public interface IAlgebraicStrategy<TMatrix> where TMatrix : class, IMatrix
    {
        public (ISolver, GlobalAlgebraicModel<TMatrix>) GetSolver(Model model);
        //GlobalAlgebraicModel<TMatrix> GetAlgebraicModel(Model model);


    }
}
