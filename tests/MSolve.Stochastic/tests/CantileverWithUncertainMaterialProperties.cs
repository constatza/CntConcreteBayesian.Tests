using System;

namespace MGroup.Stochastic.Structural.Example
{
    public class CantileverWithUncertainMaterialProperties
    {
        public void Solve()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var domainMapper = new CantileverStochasticDomainMapper(new [] { 0d, 0d, 0d });
            var realizer = new StructuralStochasticRealizer(youngModulus, domainMapper);
            var evaluator = new StructuralStochasticEvaluator(youngModulus, domainMapper);
            var m = new MonteCarlo(iterations, realizer, evaluator);
            m.Evaluate();
            Assert.Equal(1.5999674517697445E-08, MonteCarloMeanValue[0], 10);
            Assert.Equal(-2.399309224401548E-08, MonteCarloMeanValue[1], 10);
            Assert.Equal(4.3180461976960273E-12, MonteCarloStandardDeviation[0], 10);
            Assert.Equal(2.6374139668652252E-12, MonteCarloStandardDeviation[1], 10);

        }
    }
}
