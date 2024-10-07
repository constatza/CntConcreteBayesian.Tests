using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MGroup.Multiscale.SupportiveClasses;

using Xunit;

namespace MGroup.Multiscale.Tests
{
	public static class AbaqusReaderTest
	{
		[Fact]
		public static void Test()
		{
			var SpecPath = @"MsolveOutputs\Inputs";
			var InputFileName = "ReinfConcrBeam.inp";
			var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
			var pathName = Path.Combine(BasePath, SpecPath);

			string InputExtension = Path.GetExtension(InputFileName);
			string InputfileNameOnly = Path.Combine(pathName, Path.GetFileNameWithoutExtension(InputFileName));
			string inputFile = string.Format("{0}{1}", InputfileNameOnly, InputExtension);
			string[] setNames = new string[] { "Constraints", "Loads" };
			(var nodes, var elements, var sets, var renumbering) = AbaqusReader.ReadFile(inputFile, setNames);
		}
	}
}
