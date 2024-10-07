using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Providers;
using MGroup.MachineLearning;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGroup.Constitutive.Structural.MachineLearning
{
	public class NeuralNetworkPartialElasticMaterial3D : IIsotropicContinuumMaterial3D
	{
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";

		private readonly double[] strains = new double[6];
		private readonly double[] stresses = new double[6];
		private Matrix constitutiveMatrix = null;
		public double YoungModulus { get; }
		public double PoissonRatio { get; set; }
		private double[] stressesNew = new double[6];
		private double[] strainsNew = new double[6];
		private double[] incrementalStrains = new double[6];
		private INeuralNetwork neuralNetwork;
		private double[] materialParameters;
		private double elasticStrainNorm;
		//public double[] ConstParameters { get; set; }
		private GenericConstitutiveLawState currentState;

		public double[] MaterialParameters
		{
			get { return materialParameters; }
			set { materialParameters = value; }
		}

		public NeuralNetworkPartialElasticMaterial3D(INeuralNetwork neuralNetwork, double elasticModulus, double poissonRatio, double elasticStrainNorm, double[]? materialParameters = null)
		{
			this.neuralNetwork = neuralNetwork;
			if (materialParameters != null)
				this.materialParameters = materialParameters;
			else
				this.materialParameters = new double[0];
			this.YoungModulus = elasticModulus;
			this.PoissonRatio = poissonRatio;
			this.elasticStrainNorm = elasticStrainNorm;
		}

		private Matrix GetConstitutiveMatrix()
		{
			var totalStrains = strains.Copy();
			for (int i = 0; i < 6; i++)
			{
				totalStrains[i] += incrementalStrains[i];
			}
			if (totalStrains.Norm2() <= elasticStrainNorm)
			{
				return ElasticConstitutiveMatrix();
			}
			else
			{
				var neuralNetworkInput = new double[1, totalStrains.Length + materialParameters.Length];
				for (int i = 0; i < totalStrains.Length; i++)
				{
					neuralNetworkInput[0, i] = totalStrains[i];
				}
				for (int i = totalStrains.Length; i < neuralNetworkInput.Length; i++)
				{
					neuralNetworkInput[0, i] = materialParameters[i - totalStrains.Length];
				}
				var costitutiveMatrix = Matrix.CreateFromArray(neuralNetwork.EvaluateResponseGradients(neuralNetworkInput)[0]).GetSubmatrix(new int[] { 0, 1, 2, 3, 4, 5 }, new int[] { 0, 1, 2, 3, 4, 5 });
				//TEMP
				for (int i = 0; i < costitutiveMatrix.NumRows; i++)
				{
					for (int j = i; j < costitutiveMatrix.NumColumns; j++)
					{
						var temp2 = 0.5 * (costitutiveMatrix[i, j] + costitutiveMatrix[j, i]);
						costitutiveMatrix[i, j] = temp2;
						costitutiveMatrix[j, i] = temp2;
					}
				}
				//TEMP
				return costitutiveMatrix;
			}
		}

		private void CalculateNextStressStrainPoint()
		{
			var totalStrains = strains.Copy();
			for (int i = 0; i < 6; i++)
			{
				totalStrains[i] += incrementalStrains[i];
			}
			if (totalStrains.Norm2() <= elasticStrainNorm)
			{
				var elasticConstitutiveMatrix = ElasticConstitutiveMatrix();
				var stressesElastic = new double[6];
				for (int i = 0; i < 6; i++)
				{
					stressesElastic[i] = this.stresses[i];
					for (int j = 0; j < 6; j++)
					{
						stressesElastic[i] += elasticConstitutiveMatrix[i, j] * incrementalStrains[j];
					}
				}

				this.stressesNew = stressesElastic;
				this.strainsNew = totalStrains.Copy();
			}
			else
			{
				var neuralNetworkInput = new double[1, totalStrains.Length + materialParameters.Length];
				for (int i = 0; i < totalStrains.Length; i++)
				{
					neuralNetworkInput[0, i] = totalStrains[i];
				}
				for (int i = totalStrains.Length; i < neuralNetworkInput.Length; i++)
				{
					neuralNetworkInput[0, i] = materialParameters[i - totalStrains.Length];
				}
				var stressesTotal = neuralNetwork.EvaluateResponses(neuralNetworkInput);
				this.stressesNew = new double[6] { stressesTotal[0, 0], stressesTotal[0, 1], stressesTotal[0, 2], stressesTotal[0, 3], stressesTotal[0, 4], stressesTotal[0, 5] };
				this.strainsNew = totalStrains.Copy();
			}
		}

		private Matrix ElasticConstitutiveMatrix()
		{
			double fE1 = YoungModulus / (double)(1 + PoissonRatio);
			double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
			double fE3 = fE1 + fE2;
			double fE4 = fE1 * 0.5;
			var afE = Matrix.CreateZero(6, 6);
			afE.MatrixSymmetry = MatrixSymmetry.Symmetric;
			afE[0, 0] = fE3;
			afE[0, 1] = fE2;
			afE[0, 2] = fE2;
			afE[1, 0] = fE2;
			afE[1, 1] = fE3;
			afE[1, 2] = fE2;
			afE[2, 0] = fE2;
			afE[2, 1] = fE2;
			afE[2, 2] = fE3;
			afE[3, 3] = fE4;
			afE[4, 4] = fE4;
			afE[5, 5] = fE4;

			return afE;
		}

		#region IFiniteElementMaterial Members

		public int ID { get; set; }

		public bool Modified => false;

		public void ResetModified() { }

		#endregion

		#region IFiniteElementMaterial3D Members

		public double[] Stresses => stressesNew;

		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[6]);
				return constitutiveMatrix;
			}
		}

		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			//throw new NotImplementedException();
			this.incrementalStrains.CopyFrom(strainsIncrement);
			constitutiveMatrix = GetConstitutiveMatrix();
			this.CalculateNextStressStrainPoint();

			return stressesNew;
		}

		public void ClearStresses()
		{
			stresses.Clear();
			stressesNew.Clear();
		}

		#endregion

		#region ICloneable Members
		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		object ICloneable.Clone() => Clone();

		public NeuralNetworkMaterial3D Clone()
		{
			return new NeuralNetworkMaterial3D(neuralNetwork, materialParameters);
		}

		public void ClearState()
		{
			//constitutiveMatrix.Clear();
			incrementalStrains.Clear();
			stresses.Clear();
			strains.Clear();
			stressesNew.Clear();
			strainsNew.Clear();
		}

		public GenericConstitutiveLawState CreateState()
		{
			stresses.CopyFrom(stressesNew);
			strains.CopyFrom(strainsNew);
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, stresses[0]),
				(STRESS_Y, stresses[1]),
				(STRESS_Z, stresses[2]),
				(STRESS_XY, stresses[3]),
				(STRESS_XZ, stresses[4]),
				(STRESS_YZ, stresses[5]),
			});

			return currentState;
		}
		IHaveState ICreateState.CreateState() => CreateState();
		public GenericConstitutiveLawState CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}

		#endregion
	}
}
