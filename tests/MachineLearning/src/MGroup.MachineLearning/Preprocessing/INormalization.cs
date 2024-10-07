namespace MGroup.MachineLearning.Preprocessing
{
    public enum NormalizationDirection
    {
        PerRow,
        PerColumn
    }

	public interface INormalization
	{
        double[] ScalingRatio { get; }

        void Initialize(double[,] data, NormalizationDirection direction);
        double[,] Normalize(double[,] data);
        double[,] Denormalize(double[,] data);
    }
}
