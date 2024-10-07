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
	public class NeuralNetworkVonMisesMaterial3D : IIsotropicContinuumMaterial3D
	{
		private const string EQUIVALENT_STRAIN = "Equivalent strain";
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
		public double yieldStress;
		public double hardeningRatio;
		//public double[] ConstParameters { get; set; }
		private GenericConstitutiveLawState currentState;

		public double[] MaterialParameters
		{
			get { return materialParameters; }
			set { materialParameters = value; }
		}

		public NeuralNetworkVonMisesMaterial3D(INeuralNetwork neuralNetwork, double elasticModulus, double poissonRatio, double yieldStress, double hardeningRatio, double[]? materialParameters = null)
		{
			this.neuralNetwork = neuralNetwork;
			if (materialParameters != null)
				this.materialParameters = materialParameters;
			else
				this.materialParameters = new double[0];
			this.YoungModulus = elasticModulus;
			this.PoissonRatio = poissonRatio;
			this.yieldStress = yieldStress;
			this.hardeningRatio = hardeningRatio;
			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			this.bulkModulus = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			strainsEquivalentPrev = 0;
		}

		private Matrix GetConstitutiveMatrix()
		{
			var totalStrains = strains.Copy();
			for (int i = 0; i < 6; i++)
			{
				totalStrains[i] += incrementalStrains[i];
			}
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

		private void CalculateNextStressStrainPoint()
		{
			//calculate trial variables
			var stressTrial = new double[incrementalStrains.Length];
			var strainDeviatoricTrial = new double[incrementalStrains.Length];
			var normStrainDeviatoricTrial = new double();

			var totalStrains = strains.Copy();
			for (int i = 0; i < 6; i++)
			{
				totalStrains[i] += incrementalStrains[i];
			}
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

			var J2 = GetDeviatorSecondStressInvariant(stressesNew);
			double vonMisesStress = Math.Sqrt(3 * J2);
			double vonMisesStressMinusYieldStress = vonMisesStress -
													(this.yieldStress + (this.hardeningRatio * this.strainsEquivalentPrev));

			bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;

			if (materialIsInElasticRegion)
			{
				strainsEquivalent = strainsEquivalentPrev;
			}
			else
			{
				double deltastrainsEquivalent = vonMisesStressMinusYieldStress /
											((3 * this.shearModulus) + this.hardeningRatio);
				this.strainsEquivalent = this.strainsEquivalentPrev + deltastrainsEquivalent;
				//this.BuildTangentialConstitutiveMatrix();
			}
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
			this.CalculateNextStressStrainPoint();
			constitutiveMatrix = GetConstitutiveMatrix();

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

		public NeuralNetworkVonMisesMaterial3D Clone()
		{
			return new NeuralNetworkVonMisesMaterial3D(neuralNetwork, YoungModulus, PoissonRatio, yieldStress, hardeningRatio, materialParameters);
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
			strainsEquivalentPrev = strainsEquivalent;
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(EQUIVALENT_STRAIN, strainsEquivalent),
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
				strainsEquivalent = currentState.StateValues[EQUIVALENT_STRAIN];
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}

		#endregion

		private void BuildConsistentTangentialConstitutiveMatrix(double[] strainDeviatoricTrial, double normStrainDeviatoricTrial, double vonMisesStress, double deltastrainsEquivalent)
		{
			var N = new double[incrementalStrains.Length];
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				N[i] = strainDeviatoricTrial[i] / normStrainDeviatoricTrial;
			}
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					constitutiveMatrix[i, j] = 2 * shearModulus * (1 - (3 * shearModulus * vonMisesStress) / deltastrainsEquivalent) * DeviatoricProjection[i, j] +
						6 * shearModulus * shearModulus * (vonMisesStress / deltastrainsEquivalent - 1 / (3 * shearModulus + hardeningRatio)) * (N[i] * N[j]) + bulkModulus * VolumetricProjection[i, j];
				}
			}
		}

		private static readonly double[,] DeviatoricProjection = new[,]
	{
				{  2.0/3.0, -1.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, 2.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, -1.0/3.0, 2.0/3.0, 0,   0,   0  },
				{  0,  0,  0, 1.0, 0,   0   },
				{  0,  0,  0, 0,   1.0, 0   },
				{  0,  0,  0, 0,   0,   1.0 }
			};

		private static readonly double[,] VolumetricProjection = new[,]
		{
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{ 1.0, 1.0,  1.0, 0,   0,   0   },
			{  0,  0,  0, 0, 0,   0   },
			{  0,  0,  0, 0,   0, 0   },
			{  0,  0,  0, 0,   0,   0}
		};

		private readonly double shearModulus;

		private readonly double bulkModulus;

		/// <summary>
		///   The previously converged equivalent/accumulated plastic strain vector.
		/// </summary>
		private double strainsEquivalentPrev;

		/// <summary>
		///   The current equivalent/accumulated plastic strain vector.
		/// </summary>
		private double strainsEquivalent;


		public double GetThirdStressInvariant(double[] stresses)
	=> (stresses[0] * stresses[1] * stresses[2]) + (2 * stresses[5] * stresses[3] * stresses[4])
	- (Math.Pow(stresses[5], 2) * stresses[2]) - (Math.Pow(stresses[3], 2) * stresses[0])
	- (Math.Pow(stresses[4], 2) * stresses[1]);

		/// <summary>
		///   Returns the first stress invariant of the stress deviator tensor (J1), which is zero.
		/// </summary>
		/// <returns> The first stress invariant of the stress deviator tensor (J1). </returns>
		public double GetDeviatorFirstStressInvariant(double[] stresses) => 0;

		/// <summary>
		///   Calculates and returns the second stress invariant of the stress deviator tensor (J2).
		/// </summary>
		/// <returns> The second stress invariant of the stress deviator tensor (J2). </returns>
		public double GetDeviatorSecondStressInvariant(double[] stresses)
		{
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);

			double j2 = (1 / 3d * Math.Pow(i1, 2)) - i2;
			return j2;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant of the stress deviator tensor (J3).
		/// </summary>
		/// <returns> The third deviator stress invariant (J3). </returns>
		public double GetDeviatorThirdStressInvariant(double[] stresses)
		{
			double i1 = this.GetFirstStressInvariant(stresses);
			double i2 = this.GetSecondStressInvariant(stresses);
			double i3 = this.GetThirdStressInvariant(stresses);

			double j3 = (2 / 27 * Math.Pow(i1, 3)) - (1 / 3 * i1 * i2) + i3;
			return j3;
		}

		public double GetFirstStressInvariant(double[] stresses) => stresses[0] + stresses[1] + stresses[2];

		/// <summary>
		///   Calculates and returns the mean hydrostatic stress.
		/// </summary>
		/// <returns> The mean hydrostatic stress.</returns>
		public double GetMeanStress(double[] stresses) => GetFirstStressInvariant(stresses) / 3.0;

		/// <summary>
		///   Calculates and returns the second stress invariant (I2).
		/// </summary>
		/// <returns> The second stress invariant (I2).</returns>
		public double GetSecondStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2])
			- Math.Pow(stresses[5], 2) - Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);

		/// <summary>
		///   Calculates and returns the stress deviator tensor in vector form.
		/// </summary>
		/// <returns> The stress deviator tensor in vector form.</returns>
		public double[] GetStressDeviator(double[] stresses)
		{
			var hydrostaticStress = this.GetMeanStress(stresses);
			var stressDeviator = new double[]
			{
				stresses[0] - hydrostaticStress,
				stresses[1] - hydrostaticStress,
				stresses[2] - hydrostaticStress,
				stresses[3],
				stresses[4],
				stresses[5]
			};

			return stressDeviator;
		}
	}
}
