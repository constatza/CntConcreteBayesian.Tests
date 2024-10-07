//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Globalization;
//using System.IO;
//using System.Linq;
//using System.Threading;

//using MGroup.Constitutive.Structural;
//using MGroup.Constitutive.Structural.BoundaryConditions;
//using MGroup.Constitutive.Structural.Cohesive;
//using MGroup.Constitutive.Structural.Continuum;
//using MGroup.Constitutive.Structural.Line;
//using MGroup.FEM.Structural.Continuum;
//using MGroup.FEM.Structural.Embedding;
//using MGroup.FEM.Structural.Line;
//using MGroup.MSolve.Discretization;
//using MGroup.MSolve.Discretization.Dofs;
//using MGroup.MSolve.Discretization.Entities;
//using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
//using MGroup.MSolve.Numerics.Integration.Quadratures;
//using MGroup.MSolve.Numerics.Interpolation;
//using MGroup.Multiscale.SupportiveClasses;

////using Troschuetz.Random;

//namespace MGroup.Multiscale
//{
//	public class EmbeddedCntElasticNanocomposite : IRVEbuilder
//	{

//		int matrixDiscretization = 15;
//		int cntDiscretization = 10;
//		double matrixLength = 2;
//		int cntNumberY = 2;
//		int cntNumberZ = 2;
//		private double youngModulus = 1.0;//1.051e12; // 5490; // 
//		private double shearModulus = 1.0;//0.45e12; // 871; // 

//		private double[,] connectivity = null;

//		int hexa1;
//		int hexa2;
//		int hexa3;

//		double L01;
//		double L02;
//		double L03;

//		IIsotropicContinuumMaterial3D matrixMaterial;
//		IIsotropicContinuumMaterial3D inclusionMaterial;
//		ICohesiveZoneMaterial cohesiveMaterial;
//		int hostElements;
//		int embeddedElements;
//		int hostNodes;
//		int embeddedNodes;
//		// cnt paramaters
//		double cntLength;
//		double cntRadius;
//		IIsotropicContinuumMaterial3D CntMaterial;
//		int numberOfCnts;
//		// define mechanical properties
//		double poissonRatio;  //2.15; // 0.034;
//		double area; // 1776.65;  // CNT(20,20)-LinearEBE-TBT-L = 10nm
//		double inertiaY; //1058.55;
//		double inertiaZ; //1058.55;1058.55;
//		double torsionalInertia; //496.38;
//		double effectiveAreaY;
//		double effectiveAreaZ;
//		int subdomainID = 0;

//		// Cohesive Zone mechanical properties
//		//double t_max = 0.05;
//		//double K_coh = 10.0;
//		//public double K_el { get; set; }
//		//public double K_pl { get; set; }
//		//public double T_max { get; set; }

//		private double[] constParameters;

//		public Model model { get; private set; }
//		List<INode> elementNodesClone;
//		List<INode> elementNodesBeam;

//		public bool readFromText;

//		public EmbeddedCntElasticNanocomposite(int numberOfCnts, IIsotropicContinuumMaterial3D matrixMaterial, IIsotropicContinuumMaterial3D inclusionMaterial = null, ICohesiveZoneMaterial cohesiveMaterial = null)
//		{
//			var K_el = 10; var K_pl = 1; var T_max = 0.1;
//			//this.matrixMaterial = new ElasticMaterial3D(youngModulus: 3.5, poissonRatio: 0.4);
//			//constParameters = new double[3] { K_el, K_pl, T_max };
//			this.matrixMaterial = matrixMaterial;
//			this.inclusionMaterial = inclusionMaterial;
//			this.cohesiveMaterial = cohesiveMaterial;
//			CalculateInclusionProperties();
//			//this.matrixMaterial = new NeuralNetworkTrainedMaterial();// { ConstParameters = this.constParameters };

//			if (inclusionMaterial == null)
//				this.CntMaterial = new ElasticMaterial3D(youngModulus, (double)((youngModulus / (2 * shearModulus)) - 1));
//			else
//				this.CntMaterial = inclusionMaterial;

//			if (cohesiveMaterial == null)
//				this.cohesiveMaterial = new BondSlipMaterial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);
//				//this.cohesiveMaterial = new BenzeggaghKenaneCohesiveMaterial(10, 1, 10, 10, 1, 10, 1.4);
//			else
//				this.cohesiveMaterial = cohesiveMaterial;

//			double effectiveAreaY = area;
//			double effectiveAreaZ = area;

//			this.numberOfCnts = numberOfCnts;
//		}

//		private void CalculateInclusionProperties() 
//		{
//			L01 = matrixLength; L02 = matrixLength; L03 = matrixLength;
//			hexa1 = matrixDiscretization; hexa2 = matrixDiscretization; hexa3 = matrixDiscretization;
//			cntLength = 0.82 * matrixLength;
//			cntRadius = 0.015 * matrixLength;
//			area = Math.PI * Math.Pow(cntRadius, 2);
//			effectiveAreaY = area;
//			effectiveAreaZ = area;
//			inertiaY = 1 / 4 * Math.PI * Math.Pow(cntRadius, 4);
//			inertiaZ = 1 / 4 * Math.PI * Math.Pow(cntRadius, 4);
//			torsionalInertia = inertiaY + inertiaZ;

//			var nCNTy = cntNumberY - 1;
//			var nCNTz = cntNumberZ - 1;
//			var nInclusions = (nCNTy) * (nCNTz);
//			var Ldiscry = matrixLength / cntNumberY;
//			var Ldiscrz = matrixLength / cntNumberZ;

//			connectivity = new double[2*nInclusions, 6];
//			var t = 0;
//			var zor = -0.5 * matrixLength + Ldiscrz;
//			for (var i1 = 0; i1 < nCNTz; i1++)
//			{
//				var yor = -0.5 * matrixLength + Ldiscry;
//				for (var i2 = 0; i2 < nCNTy; i2++)
//				{
//					connectivity[t, 0] = -0.5 * cntLength;
//					connectivity[t, 1] = yor;
//					connectivity[t, 2] = zor;
//					connectivity[t + 1, 0] = 0.5 * cntLength;
//					connectivity[t + 1, 1] = yor;
//					connectivity[t + 1, 2] = zor;
//					t = t + 2;

//					yor = yor + Ldiscry;
//				}
//				zor = zor + Ldiscrz;
//			}

//			//For i2 In { 1:nCNTy}

//			//// THIS CHANGES ACCORDING TO nInclusions
//			//Cylinder(1 + t) = { -0.5 * Lc,  yor, zor, Lc, 0, 0, r, 2 * Pi};
//			//Physical Volume(t) = { 1 + t};
//			//t = t + 1;

//			//yor = yor + Ldiscry;
//			//EndFor
//			//zor = zor + Ldiscrz;
//			//EndFor

//			//connectivity = new double[,] { { -0.5 * cntLength, 0, 0 }, { 0.5 * cntLength, 0, 0 } };
//		}

//		public Tuple<Model, Dictionary<int, INode>, double> GetModelAndBoundaryNodes()
//		{
//			var boundaryNodes = new Dictionary<int, INode>();
//			if (model == null)
//			{
//				(int[] NodeIds, double[,] node_coords) = GetHexaRveNodesData();
//				int[,] elementConnectivity = GetHexaRveConnectivity();

//				model = new Model();
//				model.SubdomainsDictionary[0] = new Subdomain(0);
//				AddHexaElements(model, NodeIds, node_coords, elementConnectivity);

//				//var hostElements = model.EnumerateElements().Count() - 1;
//				//var hostNodes = model.EnumerateNodes().Count();

//				hostNodes = model.NodesDictionary.Count;
//				embeddedNodes = 2 * numberOfCnts;
//				hostElements = model.ElementsDictionary.Count;
//				embeddedElements = numberOfCnts;

//				(int[] cntNodeIds, double[,] cntNodeCoords, int[,] cntElementConnectivity) = GetCntBeamsNodesData(hostNodes, hostElements, readFromText);

//				AddCntBeamElements(model, cntNodeIds, cntNodeCoords, cntElementConnectivity);
//				var embeddedGrouping = EmbeddedBeam3DGrouping.CreateFullyBonded(model, model.ElementsDictionary
//				.Where(x => x.Key < hostElements).Select(kv => kv.Value).ToArray(), model.ElementsDictionary.Where(x => x.Key >= hostElements)
//				.Select(kv => kv.Value).ToArray(), false);
//				//AddCohesiveBeamElements(model, cntNodeIds, cntNodeCoords, cntElementConnectivity);
//				//var embeddedGrouping = EmbeddedBeam3DGrouping.CreateCohesive(model, model.ElementsDictionary
//				//			.Where(x => x.Key < hostElements).Select(kv => kv.Value).ToArray(), model.ElementsDictionary.Where(x => x.Key >= hostElements + embeddedElements)
//				//			.Select(kv => kv.Value).ToArray(), true);

//				//var paraviewEmbedded =
//				//new ParaviewEmbedded3D(model, null, Path.Combine(Directory.GetCurrentDirectory(), "ParaviewCNT"));
//				//paraviewEmbedded.CreateParaviewFile();

//				var boundaryNodesIds = GetBoundaryNodeIds();
//				boundaryNodesIds.Sort();
//				boundaryNodes = boundaryNodesIds.ToDictionary(t => t, t => model.NodesDictionary[t]);
//			}
//			else
//			{
//				var boundaryNodesIds = GetBoundaryNodeIds();
//				boundaryNodesIds.Sort();
//				boundaryNodes = boundaryNodesIds.ToDictionary(t => t, t => model.NodesDictionary[t]);
//			}
//			return new Tuple<Model, Dictionary<int, INode>, double>(model, boundaryNodes, L01 * L02 * L03);
//		}

//		private List<int> GetBoundaryNodeIds()
//		{
//			var boundaryNodes = new List<int>();
//			for (int k = 0; k < hexa3 + 1; k++)
//			{
//				var indexZeta = (hexa1 + 1) * (hexa2 + 1) * k;
//				if (k == 0 || k == hexa3)
//				{
//					for (int j = 0; j < hexa2 + 1; j++)
//					{
//						for (int i = 0; i < hexa1 + 1; i++)
//						{
//							boundaryNodes.Add(i + (hexa1 + 1) * j + indexZeta);
//						}
//					}
//				}
//				else
//				{
//					for (int i = 0; i < hexa1 + 1; i++)
//					{
//						boundaryNodes.Add(i + indexZeta);
//						boundaryNodes.Add(i + (hexa1 + 1) * hexa2 + indexZeta);
//					}

//					for (int j = 1; j < hexa2; j++)
//					{
//						boundaryNodes.Add(indexZeta + (hexa1 + 1) * j);
//						boundaryNodes.Add(indexZeta + (hexa1 + 1) * j + hexa1);
//					}
//				}
//			}
//			//return boundaryNodes.Distinct().ToList();
//			return boundaryNodes;
//		}

//		private void AddCntBeamElements(Model model, int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity)
//		{
//			var hostElements = model.ElementsDictionary.Count;
//			for (int i = 0; i < cntNodeIds.Length; i++)
//			{
//				var nodeID = cntNodeIds[i];
//				var coordX = cntNodeCoordinates[i, 0];
//				var coordY = cntNodeCoordinates[i, 1];
//				var coordZ = cntNodeCoordinates[i, 2];

//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: coordX, y: coordY, z: coordZ));
//			}

//			var beamSection =
//				new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

//			for (int i = 0; i < cntElementConnectivity.GetLength(0); i++)
//			{
//				var elementNodes = new List<INode>
//				{
//					model.NodesDictionary[cntElementConnectivity[i, 0]],
//					model.NodesDictionary[cntElementConnectivity[i, 1]]
//				};
//				var cntMaterialClone = (ElasticMaterial3D)this.CntMaterial;
//				var CntMaterial = new ElasticMaterial3D(cntMaterialClone.YoungModulus, cntMaterialClone.PoissonRatio);
//				//IElementType beam_1 = new Beam3DCorotationalQuaternion(elementNodes, CntMaterial.YoungModulus, CntMaterial.PoissonRatio, 7.85, beamSection);
//				IElementType beam_1 = new Rod3D(elementNodes, CntMaterial.YoungModulus) {SectionArea = area, Density = 7.85 };
//				beam_1.ID = i + hostElements;
//				//var beamElement = new Element { ID = i + hostElements, ElementType = beam_1 };

//				//beamElement.AddNode(elementNodes[0]);
//				//beamElement.AddNode(elementNodes[1]);

//				model.ElementsDictionary.Add(beam_1.ID, beam_1);
//				model.SubdomainsDictionary[0].Elements.Add(beam_1);
//			}

//		}

//		private void AddCohesiveBeamElements(Model model, int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity)
//		{
//			// define mechanical properties

//			double mi = 8.0;
//			double ni = 8.0;
//			double thickness_CNT = 0.34;
//			double a = 0.241;
//			double diameter_CNT = (a / Math.PI) * Math.Sqrt(Math.Pow(ni, 2) + ni * mi + Math.Pow(mi, 2));
//			double radius_CNT = diameter_CNT / 2.0;
//			double radius_CNT_outer = radius_CNT + (thickness_CNT / 2);
//			double CntPerimeter = 2.0 * Math.PI * radius_CNT_outer;

//			for (int i = 0; i < cntNodeIds.Length; i++)
//			{
//				var nodeID = i + hostNodes + embeddedNodes;
//				var coordX = cntNodeCoordinates[i, 0];
//				var coordY = cntNodeCoordinates[i, 1];
//				var coordZ = cntNodeCoordinates[i, 2];
//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: coordX, y: coordY, z: coordZ));
//			}

//			// Create Cohesive Material
//			//var cohesiveMaterial = new BondSlipMaterial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);
//			//var cohesiveMaterial = new BondSlipCohMatUniaxial(K_el, K_pl, 100.0, T_max, new double[2], new double[2], 1e-3);

//			// Create Beam3D Section
//			var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

//			// element nodes
//			for (int i = 0; i < cntElementConnectivity.GetLength(0); i++)
//			{

//				int elementID = i + (hostElements + embeddedElements); // matrixElements + CNTelements
//				var node1 = model.NodesDictionary[cntElementConnectivity[i, 0] + embeddedNodes].ID;
//				var node2 = model.NodesDictionary[cntElementConnectivity[i, 1] + embeddedNodes].ID;  // matrixNodes + CNTnodes

//				elementNodesClone = new List<INode>();
//				elementNodesClone.Add(model.NodesDictionary[cntElementConnectivity[i, 0]]);
//				elementNodesClone.Add(model.NodesDictionary[cntElementConnectivity[i, 1]]);
//				// element nodes beam
//				elementNodesBeam = new List<INode>();
//				elementNodesBeam.Add(model.NodesDictionary[cntElementConnectivity[i, 0]]);
//				elementNodesBeam.Add(model.NodesDictionary[cntElementConnectivity[i, 1]]);

//				List<INode> nodeSet = new List<INode>()
//				{
//					model.NodesDictionary[node1 - embeddedNodes],
//					model.NodesDictionary[node2 - embeddedNodes],
//					model.NodesDictionary[node1],
//					model.NodesDictionary[node2]
//				};
//				IElementType cohesiveElement = new CohesiveBeam3DToBeam3D(nodeSet, cohesiveMaterial, GaussLegendre1D.GetQuadratureWithOrder(2), elementNodesBeam, elementNodesClone, matrixMaterial, 1, beamSection, CntPerimeter);
//				cohesiveElement.ID = elementID;
//				//{
//				//	ID = elementID,
//				//	//ElementType = new Beam3DCorotationalQuaternion(elementNodesClone, CntMaterial, 7.85, beamSection),
//				//	ElementType = new CohesiveBeam3DToBeam3D(cohesiveMaterial, GaussLegendre1D.GetQuadratureWithOrder(2), elementNodesBeam,
//				//	elementNodesClone, matrixMaterial, 1, beamSection, CntPerimeter)
//				//};

//				// Add beam element to the element and subdomains dictionary of the model
//				model.ElementsDictionary.Add(cohesiveElement.ID, cohesiveElement);
//				// Add Cohesive Element Nodes (!)
//				//model.ElementsDictionary[cohesiveElement.ID].AddNode(model.NodesDictionary[node1 - embeddedNodes]);
//				//model.ElementsDictionary[cohesiveElement.ID].AddNode(model.NodesDictionary[node2 - embeddedNodes]);
//				//model.ElementsDictionary[cohesiveElement.ID].AddNode(model.NodesDictionary[node1]);
//				//model.ElementsDictionary[cohesiveElement.ID].AddNode(model.NodesDictionary[node2]);
//				// Add Cohesive Element in Subdomain
//				model.SubdomainsDictionary[0].Elements.Add(cohesiveElement);
//			}
//		}

//		private (int[] cntNodeIds, double[,] cntNodeCoordinates, int[,] cntElementConnectivity) GetCntBeamsNodesData(int hostNodes, int hostElements, bool readFromText)
//		{
//			if (readFromText == false)
//			{
//				var cntNodeIds = new int[numberOfCnts * 2];
//				var cntNodeCoordinates = new double[numberOfCnts * 2, 3];
//				var cntElementConnectivity = new int[numberOfCnts, 2];
				
//				var cntGenerator = new RandomCntGeometryGenerator(cntDiscretization, cntLength, numberOfCnts, L01, L02, L03, hexa1, hexa2, hexa3);
//				(cntNodeIds, cntNodeCoordinates, cntElementConnectivity) = cntGenerator.GenerateCnts(connectivity);
//				return (cntNodeIds, cntNodeCoordinates, cntElementConnectivity);
//			}
//			else
//			{
//				var cntNodeIds = new int[numberOfCnts * 2];
//				var cntNodeCoordinates = new double[numberOfCnts * 2, 3];
//				var cntElementConnectivity = new int[numberOfCnts, 2];

//				//var BasePath = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
//				var SpecPath = "E:\\Desktop\\MsolveOutputs\\matlabGeneratedCNTs\\RVE_configurations";
//				//var workingDirectory = Path.Combine(BasePath, SpecPath);
//				//string workingDirectory = @"C:\Users\stefp\OneDrive\Desktop\MsolveOutputs\matlabGeneratedCNTs";

//				string CNTgeometryFileName = "nodes_1_wf.txt";
//				string CNTconnectivityFileName = "connectivity_1_wf.txt";

//				string fileNameOnlyCNTgeometryFileName = Path.Combine(SpecPath, Path.GetFileNameWithoutExtension(CNTgeometryFileName));
//				string fileNameOnlyCNTconnectivityFileName = Path.Combine(SpecPath, Path.GetFileNameWithoutExtension(CNTconnectivityFileName));
//				string extension = Path.GetExtension(CNTgeometryFileName);
//				string extension_2 = Path.GetExtension(CNTconnectivityFileName);

//				string currentCNTconnectivityFileName = string.Format("{0}{1}", fileNameOnlyCNTconnectivityFileName, extension_2);

//				string currentCNTgeometryFileName = string.Format("{0}{1}", fileNameOnlyCNTgeometryFileName, extension);

//				//string currentCNTconnectivityFileName = string.Format("{0}_{1}{2}", fileNameOnlyCNTconnectivityFileName, noStochasticSimulation, extension_2);


//				int CNTNodes = File.ReadLines(currentCNTgeometryFileName).Count();
//				int CNTElems = File.ReadLines(currentCNTconnectivityFileName).Count();

//				// Geometry
//				using (TextReader reader = File.OpenText(currentCNTgeometryFileName))
//				{
//					string CultureName = Thread.CurrentThread.CurrentCulture.Name;
//					CultureInfo ci = new CultureInfo(CultureName);
//					ci.NumberFormat.NumberDecimalSeparator = ".";
//					Thread.CurrentThread.CurrentCulture = ci;
//					for (int i = 0; i < CNTNodes; i++)
//					{
//						string text = reader.ReadLine();
//						string[] bits = text.Split(' ');
//						int nodeID = int.Parse(bits[0]); // matrixNodes
//						double nodeX = double.Parse(bits[1]);
//						double nodeY = double.Parse(bits[2]);
//						double nodeZ = double.Parse(bits[3]);
//						cntNodeIds[i] = nodeID;
//						cntNodeCoordinates[i, 0] = nodeX; cntNodeCoordinates[i, 1] = nodeY; cntNodeCoordinates[i, 2] = nodeZ;
//						//TEMP
//						cntNodeIds[i] = hostNodes + i;
//						//TEMP

//					}
//				}
//				using (TextReader reader = File.OpenText(currentCNTconnectivityFileName))
//				{
//					for (int i = 0; i < CNTElems; i++)
//					{
//						string text = reader.ReadLine();
//						string[] bits = text.Split(' ');
//						int node1 = int.Parse(bits[0]); // matrixNodes + CNTnodes
//						int node2 = int.Parse(bits[1]); // matrixNodes + CNTnodes
//						cntElementConnectivity[i, 0] = node1; cntElementConnectivity[i, 1] = node2;
//						//TEMP
//						cntElementConnectivity[i, 0] = hostNodes + 2 * i; cntElementConnectivity[i, 1] = hostNodes + 2 * i + 1;
//						//TEMP
//					}
//				}
//				return (cntNodeIds, cntNodeCoordinates, cntElementConnectivity);
//			}
//		}

//		private void AddHexaElements(Model model, int[] nodeIds, double[,] node_coords, int[,] elementConnectivity)
//		{
//			for (int i1 = 0; i1 < nodeIds.Length; i1++)
//			{
//				int nodeID = nodeIds[i1];
//				double nodeCoordX = node_coords[i1, 0];
//				double nodeCoordY = node_coords[i1, 1];
//				double nodeCoordZ = node_coords[i1, 2];

//				model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ));
//			}

//			var renumbering = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
//			for (int i1 = 0; i1 < elementConnectivity.GetLength(0); i1++)
//			{
//				List<INode> nodeSet = new List<INode>();
//				for (int j = 0; j < 8; j++)
//				{
//					int nodeID = elementConnectivity[i1, renumbering[j]];
//					nodeSet.Add(model.NodesDictionary[nodeID]);
//				}
//				var elementFactory = new Multiscale.SupportiveClasses.ContinuumElement3DFactory((IContinuumMaterial3D)matrixMaterial.Clone(), commonDynamicProperties: null);
//				IElementType e1 = elementFactory.CreateElement(CellType.Hexa8, nodeSet);
//				e1.ID = i1;

//				model.ElementsDictionary.Add(e1.ID, e1);
//				model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
//			}
//		}



//		private int[,] GetHexaRveConnectivity()
//		{
//			int numberOfElements = hexa1 * hexa2 * hexa3;
//			int[,] elementConnectivity = new int[numberOfElements, 8];

//			int counterElement = 0;
//			for (int i = 0; i < hexa3; i++)
//			{
//				for (int j = 0; j < hexa2; j++)
//				{
//					for (int k = 0; k < hexa1; k++)
//					{
//						elementConnectivity[counterElement, 0] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * i;
//						elementConnectivity[counterElement, 1] = elementConnectivity[counterElement, 0] + 1;
//						elementConnectivity[counterElement, 2] = elementConnectivity[counterElement, 1] + (hexa1 + 1);
//						elementConnectivity[counterElement, 3] = elementConnectivity[counterElement, 2] - 1;

//						elementConnectivity[counterElement, 4] = k + (hexa1 + 1) * j + (hexa1 + 1) * (hexa2 + 1) * (i + 1);
//						elementConnectivity[counterElement, 5] = elementConnectivity[counterElement, 4] + 1;
//						elementConnectivity[counterElement, 6] = elementConnectivity[counterElement, 5] + (hexa1 + 1);
//						elementConnectivity[counterElement, 7] = elementConnectivity[counterElement, 6] - 1;

//						counterElement++;
//					}
//				}
//			}

//			return elementConnectivity;
//		}

//		private (int[] NodeIds, double[,] node_coords) GetHexaRveNodesData()
//		{
//			int numberOfNodes = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
//			int[] NodeIds = new int[numberOfNodes];
//			double[,] node_coords = new double[numberOfNodes, 3];

//			int counterNode = 0;

//			for (int i = 0; i <= hexa3; i++)
//			{
//				for (int j = 0; j <= hexa2; j++)
//				{
//					for (int k = 0; k <= hexa1; k++)
//					{
//						NodeIds[counterNode] = counterNode;
//						node_coords[counterNode, 0] = k * (L01 / hexa1) - L01 * 0.5;
//						node_coords[counterNode, 1] = j * (L02 / hexa2) - L02 * 0.5;
//						node_coords[counterNode, 2] = i * (L03 / hexa3) - L03 * 0.5;
//						counterNode++;
//					}
//				}
//			}

//			return (NodeIds, node_coords);
//		}

//		public NodalDisplacement[] GetModelRigidBodyNodeConstraints(Model model)
//		{
//			NodalDisplacement[] rigidBodyNodeConstraints = new NodalDisplacement[9]
//			{
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationZ, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[1], StructuralDof.TranslationZ, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationX, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationY, amount: 0),
//				new NodalDisplacement(model.NodesDictionary[2], StructuralDof.TranslationZ, amount: 0),
//			};
//			//Dictionary<Node, IList<IDofType>> RigidBodyNodeConstraints = new Dictionary<Node, IList<IDofType>>();
//			//var rigidNodes = RigidNodes;

//			//RigidBodyNodeConstraints.Add(model.NodesDictionary[rigidNodes[0]], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
//			//RigidBodyNodeConstraints.Add(model.NodesDictionary[rigidNodes[1]], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });
//			//RigidBodyNodeConstraints.Add(model.NodesDictionary[rigidNodes[2]], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });

//			return rigidBodyNodeConstraints;
//		}

//		private List<int> RigidNodes => new List<int>
//		{
//			{ 0},
//			{hexa1 },
//			{(hexa1+1)*hexa2 }
//		};

//		public IRVEbuilder Clone(int a)
//		{
//			//    var random= new TRandom();
//			//    var randomCnts = (int)random.Normal(numberOfCnts, 100);
//			//    var cnts=randomCnts<0?1:randomCnts;
//			var cnts = numberOfCnts;
//			return new EmbeddedCntElasticNanocomposite(cnts, matrixMaterial, inclusionMaterial, cohesiveMaterial);
//		}

//		public object Clone()
//		{
//			var rveBuilder = new EmbeddedCntElasticNanocomposite(numberOfCnts, matrixMaterial, inclusionMaterial, cohesiveMaterial)
//			{
//			};
//			rveBuilder.model = (Model)this.model.Clone();
//			return rveBuilder;
//		}
//	}
//}
