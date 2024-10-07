using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http.Headers;

namespace MGroup.CntConcreteBayesian.Tests.Commons
{

	public static class AbaqusReader
	{
		public static (double[][,] nodes, int[][,] elements, int[][] sets) ReadFile(string textFile, string[] setNames = null)
		{
			var text = File.ReadAllText(textFile);
			text = text.ToLower();
			var subtexts = text.Split('*');
			var indAsterisk = AllIndexesOf(text, "*");
			var indNodes = AllIndexesOf(text, "*node");
			var indElements = AllIndexesOf(text, "*element");
			var instances = indNodes.Count;
			var sets = new int[setNames.Length][];
			if (setNames != null)
			{
				for (int i = 0; i < setNames.Length; i++)
				{
					var currentSetPos = AllIndexesOf(text, $"*nset, nset={setNames[i].ToLower()}");
					var currentSet = subtexts[indAsterisk.IndexOf(currentSetPos[0]) + 1];
					var currentSetArray = To2D(GetMatrixFromString(currentSet, "\n", ",")).OfType<double>().ToArray();
					sets[i] = RemoveZeros(currentSetArray);
				}
			}
			var nodes = new double[instances][,];
			var elements = new int[instances][,];
			var renumbering = new int[instances][];
			for (int i = 0; i < instances; i++)
			{
				var elementsSet = subtexts[indAsterisk.IndexOf(indElements[i]) + 1];
				var nodesSet = subtexts[indAsterisk.IndexOf(indNodes[i]) + 1];
				var elementsTemp = To2D(GetMatrixFromString(elementsSet, "\n", ","));
				elements[i] = new int[elementsTemp.GetLength(0), elementsTemp.GetLength(1)];
				for (int k = 0; k < elementsTemp.GetLength(0); k++)
				{
					for (int l = 0; l < elementsTemp.GetLength(1); l++)
					{
						elements[i][k, l] = (int)elementsTemp[k, l];
					}
				}
				nodes[i] = To2D(GetMatrixFromString(nodesSet, "\n", ","));
				var n1 = AllIndexesOf(elementsSet, "type=");
				var n2 = AllIndexesOf(elementsSet, "\n");
				var eT = elementsSet.Substring(n1[0] + 5, 4);
				//renumbering[i] = NodeRenumbering(eT);
			}
			return (nodes, elements, sets);
		}

		private static List<int> AllIndexesOf(this string str, string value)
		{
			if (String.IsNullOrEmpty(value))
				throw new ArgumentException("the string to find may not be empty", "value");
			List<int> indexes = new List<int>();
			for (int index = 0; ; index += value.Length)
			{
				index = str.IndexOf(value, index);
				if (index == -1)
					return indexes;
				indexes.Add(index);
			}
		}

		private static int[] RemoveZeros(double[] source)
		{
			var array = source.Where(x => x != 0).ToArray();
			return array.Select(x => (int)x).ToArray();
		}

		private static double[][] GetMatrixFromString(this string str, string sepRow, string sepCol)
		{
			var strSepRow = str.Split(sepRow);
			var tempMatrix = new double[strSepRow.Length][];
			var flag = new int[strSepRow.Length][];
			var totalNaN = new int[strSepRow.Length];
			var lineFlag = new int[strSepRow.Length];
			for (int i = 0; i < strSepRow.Length; i++)
			{
				var strSepCol = strSepRow[i].Split(sepCol);
				flag[i] = new int[strSepCol.Length];
				tempMatrix[i] = new double[strSepCol.Length];
				for (int j = 0; j < strSepCol.Length; j++)
				{
					double number;
					if (Double.TryParse(strSepCol[j], out number) == true)
						tempMatrix[i][j] = Double.Parse(strSepCol[j]);
					else
						tempMatrix[i][j] = double.NaN;
					if (Double.IsNaN(tempMatrix[i][j]))
					{
						flag[i][j] = 1;
						totalNaN[i]++;
					}
				}
				if (totalNaN[i] < tempMatrix[i].Length)
					lineFlag[i] = 1;
			}
			//var totalRowsKeep = tempMatrix.GetLength(0) - flag.Sum();
			var matrix = new double[lineFlag.Sum()][];
			var count = 0;
			for (int i = 0; i < strSepRow.Length; i++)
			{
				if (lineFlag[i] == 1)
				{
					matrix[count] = new double[tempMatrix[i].Length - totalNaN[i]];
					for (int j = 0; j < tempMatrix[i].Length; j++)
					{
						if (Double.IsNaN(tempMatrix[i][j]) == false)
							matrix[count][j] = tempMatrix[i][j];
					}
					count++;
				}
			}
			return matrix;
		}

		//private static int[] NodeRenumbering(string elementType)
		//{
		//	//renumbering for other element types should be added
		//	var renumbering = new int[8];
		//	if (elementType == "c3d8")
		//	{
		//		renumbering = new int[] { 4, 8, 5, 1, 3, 7, 6, 2 };
		//	}
		//	else
		//	{
		//		renumbering = new int[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		//	}
		//	return renumbering;
		//}

		private static double[,] To2D(double[][] source)
		{
			int FirstDim = source.Length;
			int SecondDim = source[0].Length;

			var result = new double[FirstDim, SecondDim];
			for (int i = 0; i < FirstDim; i++)
			{
				var tempSecondDim = source[i].Length;
				for (int j = 0; j < tempSecondDim; j++)
					result[i, j] = source[i][j];
				if (tempSecondDim < SecondDim)
					for (int j = tempSecondDim; j < SecondDim; j++)
					{
						result[i, j] = 0;
					}
			}

			return result;
		}
	}
}
