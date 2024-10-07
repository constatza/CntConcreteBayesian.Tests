using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using CSparse.Storage;
using System.Xml.Linq;
//using Troschuetz.Random;

namespace MGroup.Multiscale.SupportiveClasses
{
	public class RandomCntGeometryGenerator
	{
		private readonly int numberOfElementsPerCnt;
		private readonly double cntLength;
		private readonly double elementLength;
		private readonly int numberOfCnts;
		private readonly double standardDeviation;
		private readonly double upperAngleBound;
		private readonly double matrixLength;
		private readonly double matrixWidth;
		private readonly double matrixHeight;
		private readonly int numElemLength;
		private readonly int numElemWidth;
		private readonly int numElemHeight;
		private int iNode;
		private double[] randVector;

		public bool periodicInclusions = false;
		public double[,] coordinates;

		public RandomCntGeometryGenerator(int numberOfElementsPerCnt, double cntLength, int numberOfCnts,
		double matrixLength, double matrixWidth, double matrixHeight, int numElemLength, int numElemWidth, int numElemHeight)
		{
			this.numberOfElementsPerCnt = numberOfElementsPerCnt;
			this.cntLength = cntLength;
			if (numberOfElementsPerCnt <= 0)
				throw new ArgumentOutOfRangeException("The number of CNTS must be greater than zero");

			this.elementLength = cntLength / (double)numberOfElementsPerCnt;
			this.numberOfCnts = numberOfCnts;
			this.matrixLength = matrixLength;
			this.matrixWidth = matrixWidth;
			this.matrixHeight = matrixHeight;
			this.numElemLength = numElemLength;
			this.numElemWidth = numElemWidth;
			this.numElemHeight = numElemHeight;
		}

		public (int[] nodeIds, double[,] nodeCoordinates, int[,] elementConnectivity) GenerateCnts(double[,] coordinates = null, bool periodicInclusions = false, int seed = 0)
		{
			var random = new Random(seed);
			var numberOfNodesPerCnt = numberOfElementsPerCnt + 1;
			var nodeIds = new int[numberOfCnts * numberOfNodesPerCnt];
			var elementConnectivity = new int[numberOfCnts * numberOfElementsPerCnt, 2];
			var nodalCoordinates = new double[numberOfCnts * numberOfNodesPerCnt, 3];

			var countNode = (numElemLength + 1) * (numElemWidth + 1) * (numElemHeight + 1);
			var position = 0;
			for (int indexCnt = 0; indexCnt < numberOfCnts; indexCnt++)
			{
				bool state = true;
				var iNode0 = indexCnt * numberOfNodesPerCnt;
				var iNode1 = 2 * indexCnt;
				//nodalCoordinates[iNode0] = new double[3];
				if (coordinates != null)
				{
					nodalCoordinates[iNode0, 0] = coordinates[iNode1, 0];
					nodalCoordinates[iNode0, 1] = coordinates[iNode1, 1];
					nodalCoordinates[iNode0, 2] = coordinates[iNode1, 2];
					randVector = new double[coordinates.GetLength(1)];
					randVector[0] = coordinates[iNode1 + 1, 0] - coordinates[iNode1, 0];
					randVector[1] = coordinates[iNode1 + 1, 1] - coordinates[iNode1, 1];
					randVector[2] = coordinates[iNode1 + 1, 2] - coordinates[iNode1, 2];
					var normValue = Math.Sqrt(randVector[0] * randVector[0] + randVector[1] * randVector[1] + randVector[2] * randVector[2]);
					randVector[0] = randVector[0] / normValue; randVector[1] = randVector[1] / normValue; randVector[2] = randVector[2] / normValue;
				}
				else
				{
					nodalCoordinates[iNode0, 0] = random.NextDouble() * (matrixLength / 2 - (-matrixLength / 2)) + (-matrixLength / 2);
					nodalCoordinates[iNode0, 1] = random.NextDouble() * (matrixHeight / 2 - (-matrixHeight / 2)) + (-matrixHeight / 2);
					nodalCoordinates[iNode0, 2] = random.NextDouble() * (matrixWidth / 2 - (-matrixWidth / 2)) + (-matrixWidth / 2);
					var randNormalVec = new double[3];
					for (int k = 0; k < 3; k++)
					{
						randNormalVec[k] = Math.Sqrt(-2.0 * Math.Log(random.NextDouble())) * Math.Sin(2.0 * Math.PI * random.NextDouble());
					}
					randVector = new double[3] { randNormalVec[0], randNormalVec[1], randNormalVec[2] };
					var normValue = Math.Sqrt(randVector[0] * randVector[0] + randVector[1] * randVector[1] + randVector[2] * randVector[2]);
					randVector[0] = randVector[0] / normValue; randVector[1] = randVector[1] / normValue; randVector[2] = randVector[2] / normValue;
					var xNode = nodalCoordinates[iNode0, 0] + randVector[0] * cntLength;
					var yNode = nodalCoordinates[iNode0, 1] + randVector[1] * cntLength;
					var zNode = nodalCoordinates[iNode0, 2] + randVector[2] * cntLength;
					if (xNode > matrixLength / 2 || xNode < -matrixLength / 2 ||
						yNode > matrixHeight / 2 || yNode < -matrixHeight / 2 ||
						zNode > matrixWidth / 2 || zNode < -matrixWidth / 2)
					{
						indexCnt = indexCnt - 1;
						continue;
					}
				}

				for (int indexElement = 0; indexElement < numberOfElementsPerCnt; indexElement++)
				{
					iNode = indexCnt * numberOfNodesPerCnt + indexElement;
					var xNode = new double(); var yNode = new double(); var zNode = new double();
					if (coordinates != null)
					{
						if (periodicInclusions == false)
						{
							xNode = nodalCoordinates[iNode, 0] + randVector[0] * cntLength / numberOfElementsPerCnt;
							yNode = nodalCoordinates[iNode, 1] + randVector[1] * cntLength / numberOfElementsPerCnt;
							zNode = nodalCoordinates[iNode, 2] + randVector[2] * cntLength / numberOfElementsPerCnt;
						}
						else
						{
							(xNode, yNode, zNode) = CreatePeriodicInclusion(nodalCoordinates, randVector);
						}
					}
					else
					{
						if (periodicInclusions == false)
						{
							xNode = nodalCoordinates[iNode, 0] + randVector[0] * cntLength / numberOfElementsPerCnt;
							yNode = nodalCoordinates[iNode, 1] + randVector[1] * cntLength / numberOfElementsPerCnt;
							zNode = nodalCoordinates[iNode, 2] + randVector[2] * cntLength / numberOfElementsPerCnt;
						}
						else
						{
							(xNode, yNode, zNode) = CreatePeriodicInclusion(nodalCoordinates, randVector);
						}
					}
					//var distance = Math.Sqrt(dx * dx + dy * dy + dz * dz);
					var elementId = indexCnt * numberOfElementsPerCnt + indexElement;
					elementConnectivity[elementId, 0] = countNode;
					elementConnectivity[elementId, 1] = countNode + 1;
					nodalCoordinates[iNode + 1, 0] = xNode;
					nodalCoordinates[iNode + 1, 1] = yNode;
					nodalCoordinates[iNode + 1, 2] = zNode;
					nodeIds[position] = countNode;
					countNode++;
					position++;					
				}
				nodeIds[position] = countNode;
				countNode++;
				position++;
			}
			return (nodeIds, nodalCoordinates, elementConnectivity);
		}

		private (double xNode, double yNode, double zNode) CreatePeriodicInclusion(double[,] startNodeCoordinates, double[] randVector)
		{
			var distFromBoundary = new double[3, 2];
			var auxNodeCoordinates = new double[3] { startNodeCoordinates[iNode, 0], startNodeCoordinates[iNode, 1], startNodeCoordinates[iNode, 2] };
			var auxElementLength = cntLength / numberOfElementsPerCnt;
			var figureElements = new List<double[]>();
			var xNode = auxNodeCoordinates[0] + randVector[0] * cntLength / numberOfElementsPerCnt;
			var yNode = auxNodeCoordinates[1] + randVector[1] * cntLength / numberOfElementsPerCnt;
			var zNode = auxNodeCoordinates[2] + randVector[2] * cntLength / numberOfElementsPerCnt;
			var done = false;
			while (done == false)
			{
				distFromBoundary[0, 0] = (matrixLength / 2 - auxNodeCoordinates[0]) * 1 / ((randVector[0] * auxElementLength) * 1);
				distFromBoundary[0, 1] = (-matrixLength / 2 - auxNodeCoordinates[0]) * 1 / ((randVector[0] * auxElementLength) * 1);
				distFromBoundary[1, 0] = (matrixWidth / 2 - auxNodeCoordinates[1]) * 1 / ((randVector[1] * auxElementLength) * 1);
				distFromBoundary[1, 1] = (-matrixWidth / 2 - auxNodeCoordinates[1]) * 1 / ((randVector[1] * auxElementLength) * 1);
				distFromBoundary[2, 0] = (matrixHeight / 2 - auxNodeCoordinates[2]) * 1 / ((randVector[2] * auxElementLength) * 1);
				distFromBoundary[2, 1] = (-matrixHeight / 2 - auxNodeCoordinates[2]) * 1 / ((randVector[2] * auxElementLength) * 1);
				double minDist = 1;
				for (int i = 0; i < distFromBoundary.GetLength(1); i++)
				{
					for (int j = 0; j < distFromBoundary.GetLength(2); j++)
					{
						if (distFromBoundary[i, j] > 1e-12 && distFromBoundary[i, j] < minDist)
						{
							minDist = distFromBoundary[i, j];
						}
					}
				}
				if (distFromBoundary[0, 0] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] - matrixLength + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else if (distFromBoundary[0, 1] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] + matrixLength + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else if (distFromBoundary[1, 0] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] - matrixWidth + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else if (distFromBoundary[1, 1] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] + matrixWidth + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else if (distFromBoundary[2, 0] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] - matrixHeight + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else if (distFromBoundary[2, 1] == minDist)
				{
					auxNodeCoordinates[0] = auxNodeCoordinates[0] + minDist * randVector[0] * auxElementLength;
					auxNodeCoordinates[1] = auxNodeCoordinates[1] + minDist * randVector[1] * auxElementLength;
					auxNodeCoordinates[2] = auxNodeCoordinates[2] + matrixHeight + minDist * randVector[2] * auxElementLength;
					auxElementLength = auxElementLength * (1 - minDist);
				}
				else
				{
					done = true;
					xNode = auxNodeCoordinates[0] + randVector[0] * auxElementLength;
					yNode = auxNodeCoordinates[1] + randVector[1] * auxElementLength;
					zNode = auxNodeCoordinates[2] + randVector[2] * auxElementLength;
				}
			}
			return (xNode, yNode, zNode);
		}
	}
}
