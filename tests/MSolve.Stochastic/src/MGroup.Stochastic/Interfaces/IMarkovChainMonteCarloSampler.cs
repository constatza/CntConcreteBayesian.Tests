namespace MGroup.Stochastic
{
    public interface IProbabilityDistributionSampler
    {
        public double[,] GenerateSamples(int numSamples);
    }
}
