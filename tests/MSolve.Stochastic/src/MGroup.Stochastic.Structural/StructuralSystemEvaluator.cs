using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization.Entities;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;

namespace MGroup.Stochastic.Structural
{
    public class StructuralSystemEvaluator : ISystemRealizer, ISystemResponseEvaluator
    {
        public double YoungModulus { get; }
        public IStochasticDomainMapper DomainMapper;
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        public GiannisModelBuilder ModelBuilder { get; }
        private Model currentModel;
        int karLoeveTerms = 4;
        double[] domainBounds = new double[2] { 0, 1 };
        double sigmaSquare = .01;
        double meanValue = 1;
        int partition = 21;
        double correlationLength = 1.0;
        bool isGaussian = true;
        int PCorder = 1;
        bool midpointMethod = true;

        /// <summary>Initializes a new instance of the <see cref="StructuralSystemEvaluator"/> class.</summary>
        /// <param name="youngModulus">The young modulus.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        public StructuralSystemEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            DomainMapper = domainMapper;
            ModelBuilder = new GiannisModelBuilder();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider(partition, youngModulus, midpointMethod,
                isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);
        }

        /// <summary>Realizes the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, DomainMapper, iteration);
        }



        /// <summary>Evaluates the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <returns></returns>
        public double[] Evaluate(int iteration)
        {
            var solverFactory = new SkylineSolver.Factory();
			var algebraicModel = solverFactory.BuildAlgebraicModel(currentModel);
			var solver = solverFactory.BuildSolver(algebraicModel);
			var problem = new ProblemStructural(currentModel, algebraicModel);
			var linearAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var staticAnalyzer = new StaticAnalyzer(algebraicModel, problem, linearAnalyzer);

			staticAnalyzer.Initialize();
			staticAnalyzer.Solve();
			return new[] { solver.LinearSystem.Solution.SingleVector[56], solver.LinearSystem.Solution.SingleVector[58] };
        }

    }
}
