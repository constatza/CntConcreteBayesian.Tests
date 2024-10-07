using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.Solvers.AlgebraicModel;
using MGroup.Solvers.Direct;
using MGroup.Solvers.DofOrdering;
using MGroup.Solvers.DofOrdering.Reordering;

namespace MiMsolve.SolutionStrategies
{
    public class SuiteSparseSolverPrefernce : IAlgebraicStrategy<SymmetricCscMatrix>
    {
        
        public (ISolver, GlobalAlgebraicModel<SymmetricCscMatrix>) GetSolver(Model model)
        {
            var solverFactory = new SuiteSparseSolver.Factory();
            solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd());
            var globalAlgebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(globalAlgebraicModel);
            return (solver, globalAlgebraicModel);
        }

        //GlobalAlgebraicModel<SkylineMatrix> IAlgebraicStrategy<SkylineMatrix>.GetAlgebraicModel(Model model)
        //{
        //    throw new NotImplementedException();
        //}

        //ISolver IAlgebraicStrategy<SkylineMatrix>.GetSolver(GlobalAlgebraicModel<SkylineMatrix> algebraicModel)
        //{
        //    throw new NotImplementedException();
        //}
    }

    //public class SkylineSolverPrefernce2<TMatrix> : IAlgebraicStrategy<TMatrix> where TMatrix : class, IMatrix
    //{
    //    public GlobalAlgebraicModel<TMatrix> GetAlgebraicModel(Model model)
    //    {
    //        throw new NotImplementedException();
    //    }

    //    public ISolver GetSolver(GlobalAlgebraicModel<TMatrix> algebraicModel)
    //    {
    //        throw new NotImplementedException();
    //    }
    //}

    //public class metodChecks<TMatrix> where TMatrix : class, IMatrix
    //{
    //    public void GetAlgebraicModel( IAlgebraicStrategy<TMatrix> solverProvider)
    //    {
    //        throw new NotImplementedException();
    //    }

    //    public static void Check()
    //    {
        
    //    }
    //}

    //public static class Checks
    //{
    //    public static void check01()
    //    {
    //        var model = new Model();
    //        var c1 = new SkylineSolverPrefernce();
    //        var c2 = new metodChecks<SymmetricCscMatrix>();
    //        c2.GetAlgebraicModel(c1); // prosoxh to c1 exei mesa tou implied to <SkylineMatrix> giati ulopoiei to interface : IAlgebraicStrategy<SkylineMatrix>
    //        // epishs to c2 tou leme kathara oti prokeitai gia thn tupou <SkylineMatrix> ulopoihsh tou methodChecks<TMatrix>
    //    }
    //}
}
