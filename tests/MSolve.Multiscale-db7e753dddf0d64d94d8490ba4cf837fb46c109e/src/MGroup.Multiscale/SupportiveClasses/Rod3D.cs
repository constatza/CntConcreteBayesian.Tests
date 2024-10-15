using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.MSolve.Constitutive;

namespace MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses
{
	public class Rod3D : IStructuralElementType, IEmbeddedElement
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double youngModulus;
		private double[] currentDisplacements; 
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();

		public double Density { get; set; }
		public double SectionArea { get; set; }
		public IList<EmbeddedNode> EmbeddedNodes { get { return embeddedNodes; } }

		public Rod3D(IReadOnlyList<INode> nodes, double youngModulus)
		{
			this.youngModulus = youngModulus;
			Nodes = nodes;
		}

		public Rod3D(IReadOnlyList<INode> nodes,  double youngModulus, IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Line2;
		public IReadOnlyList<IConstitutiveLawWithGenericState> MaterialsAtGaussPoints { get ; }
		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		//TODO: this should be either cached, or even better the calculations should be incorporated into Stiffness()
		public IMatrix TransformationMatrix()
		{
			double x2 = Math.Pow(Nodes[1].X - Nodes[0].X, 2);
			double y2 = Math.Pow(Nodes[1].Y - Nodes[0].Y, 2);
			double z2 = Math.Pow(Nodes[1].Z - Nodes[0].Z, 2);
			double L = Math.Sqrt(x2 + y2 + z2);
			double c = (Nodes[1].X - Nodes[0].X) / L;
			double s = (Nodes[1].Y - Nodes[0].Y) / L;
			double q = (Nodes[1].Z - Nodes[0].Z) / L;

			// T = [ cos sin 0 0; 0 0 cos sin]
			var transformation = Matrix.CreateZero(2, 6);
			transformation[0, 0] = c;
			transformation[0, 1] = s;
			transformation[0, 2] = q;
			transformation[1, 3] = c;
			transformation[1, 4] = s;
			transformation[1, 5] = q;
			return transformation;
		}

		/// <summary>
		/// Stress0         Stress1
		/// -> ------------ ->
		/// </summary>
		/// <param name="element"></param>
		/// <param name="localDisplacements"></param>
		/// <param name="local_d_Displacements"></param>
		/// <returns></returns>
		public double CalculateAxialStress( double[] localDisplacements, double[] local_d_Displacements)
		{
			double[] globalStresses = CalculateResponse(localDisplacements).Item2; // item1 = strains
			IMatrix transformation = TransformationMatrix();
			double[] localStresses = transformation.Multiply(globalStresses); // In local natural system there are 2 dofs
																			  // If Stress1 = localStresses[1] > 0 => tension. Else compression
			return localStresses[1];
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		public IMatrix StiffnessMatrix()
		{
			double x2 = Math.Pow(Nodes[1].X - Nodes[0].X, 2);
			double y2 = Math.Pow(Nodes[1].Y - Nodes[0].Y, 2);
			double z2 = Math.Pow(Nodes[1].Z - Nodes[0].Z, 2);
			double L = Math.Sqrt(x2 + y2 + z2);
			double c = (Nodes[1].X - Nodes[0].X) / L;
			double c2 = c * c;
			double s = (Nodes[1].Y - Nodes[0].Y) / L;
			double s2 = s * s;
			double q = (Nodes[1].Z - Nodes[0].Z) / L;
			double q2 = q * q;
			double cs = c * s;
			double cq = c * q;
			double sq = s * q;
			double E = this.youngModulus;
			double A = SectionArea;
			return dofEnumerator.GetTransformedMatrix(
				Matrix.CreateFromArray(new double[,]
				{
					{A*E*c2/L, A*E*cs/L, A*E*cq/L, -A*E*c2/L, -A*E*cs/L, -A*E*cq/L, },
					{A*E*cs/L, A*E*s2/L, A*E*sq/L, -A*E*cs/L, -A*E*s2/L, -A*E*sq/L, },
					{A*E*cq/L, A*E*sq/L, A*E*q2/L, -A*E*cq/L, -A*E*sq/L, -A*E*q2/L, },
					{-A*E*c2/L, -A*E*cs/L, -A*E*cq/L, A*E*c2/L, A*E*cs/L, A*E*cq/L, },
					{-A*E*cs/L, -A*E*s2/L, -A*E*sq/L, A*E*cs/L, A*E*s2/L, A*E*sq/L, },
					{-A*E*cq/L, -A*E*sq/L, -A*E*q2/L, A*E*cq/L, A*E*sq/L, A*E*q2/L, },
				}));
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IMatrix MassMatrix()
		{
			double x2 = Math.Pow(Nodes[1].X - Nodes[0].X, 2);
			double y2 = Math.Pow(Nodes[1].Y - Nodes[0].Y, 2);
			double z2 = Math.Pow(Nodes[1].Z - Nodes[0].Z, 2);
			double L = Math.Sqrt(x2 + y2 + z2);

			double totalMassOver2 = Density * SectionArea * L / 2.0;

			// Lumped mass: M = [m/2 0 0 0 0 0; 0 m/2 0 0 0 0; 0 0 m/2 0 0 0; 0 0 0 m/2 0 0; 0 0 0 0 m/2 0 ; 0 0 0 0 0 m/2]
			int order = 6;
			var lumpedMass = Matrix.CreateZero(order, order);
			for (int i = 0; i < order; ++i) lumpedMass[i, i] = totalMassOver2;
			return lumpedMass;
		}

		public IMatrix DampingMatrix() => throw new NotImplementedException();

		private double[] CalculateResponseIntegral(double[] localdDisplacements)
		{
			if (localdDisplacements != null)
			{
				IMatrix stiffness = StiffnessMatrix();
				return stiffness.Multiply(localdDisplacements);
			}
			else
			{
				return null;
			}
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CalculateResponseIntegral(local_Displacements);
			double[] stresses = Array.ConvertAll(forces, x => x / SectionArea);
			if (currentDisplacements == null || currentDisplacements.Length != local_Displacements.Length)
			{
				currentDisplacements = new double[local_Displacements.Length];
			}

			Array.Copy(local_Displacements, currentDisplacements, local_Displacements.Length);

			return new Tuple<double[], double[]>(strains, stresses);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public double[] CalculateResponseIntegral() => CalculateResponseIntegral( currentDisplacements);

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	var accelerations = new double[4];
		//	IMatrix massMatrix = MassMatrix(element);

		//	int index = 0;
		//	foreach (MassAccelerationLoad load in loads)
		//		foreach (IDofType[] nodalDOFTypes in dofs)
		//			foreach (IDofType dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}

		//	return massMatrix.Multiply(accelerations);
		//}

		public void SaveConstitutiveLawState(IHaveState externalState) { }

		#endregion

		#region IFiniteElement Members


		public bool ConstitutiveLawModified => false;

		public void ResetConstitutiveLawModified() { }

		#endregion

		#region IFiniteElement Members

		public void ClearConstitutiveLawState() { }

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();

		#endregion

		public Dictionary<IDofType, int> GetInternalNodalDOFs(IElementType element, INode node)
		{
			throw new NotImplementedException();
		}

		public double[] GetLocalDOFValues(IElementType hostElement, double[] hostDOFValues)
		{
			throw new NotImplementedException();
		}

		public IEnumerable<IEnumerable<double>> InterpolateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities, IEnumerable<double[]> coordinates) => throw new NotImplementedException();
		public IEnumerable<double[]> IntegrateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities) => throw new NotImplementedException();
	}
}
