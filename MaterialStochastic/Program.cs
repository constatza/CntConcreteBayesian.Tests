// See https://aka.ms/new-console-template for more information

using MGroup.CntConcreteBayesian.Tests.Integration;

namespace MaterialStochastic
{
	class Program
	{
		static void Main(string[] args)
		{
			var dataPath = args.Length > 0 ? args[0] : "../";
			int numSamples = args.Length > 1 ? int.Parse(args[1]) : 10;
			string whichModel = args.Length > 2 ? args[2] : "cement";
			
			if (!Directory.Exists(dataPath))
			{
				Console.WriteLine("The specified data path does not exist: " + dataPath);
				return;
			}
			Console.WriteLine($"Running {whichModel} Monte-Carlo simulation with {numSamples} samples.");
			ConcreteMonteCarloTest.RunSimulation(dataPath, numSamples);
            Console.WriteLine("Cement samples generated. Analysis complete.");
		}
	}
}