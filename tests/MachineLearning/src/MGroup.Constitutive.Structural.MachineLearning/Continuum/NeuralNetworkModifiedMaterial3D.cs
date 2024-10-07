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
	public class NeuralNetworkModifiedMaterial3D : IIsotropicContinuumMaterial3D
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
		private double plasticStrainNorm;
		public double yieldStress;
		public double hardeningRatio;
		//public double[] ConstParameters { get; set; }
		private GenericConstitutiveLawState currentState;

		public double[] MaterialParameters
		{
			get { return materialParameters; }
			set { materialParameters = value; }
		}

		public NeuralNetworkModifiedMaterial3D(INeuralNetwork neuralNetwork, double elasticModulus, double poissonRatio, double elasticStrainNorm, double yieldStress, double hardeningRatio, double plasticStrainNorm, double[]? materialParameters = null)
		{
			this.neuralNetwork = neuralNetwork;
			if (materialParameters != null)
				this.materialParameters = materialParameters;
			else
				this.materialParameters = new double[0];
			this.YoungModulus = elasticModulus;
			this.PoissonRatio = poissonRatio;
			this.elasticStrainNorm = elasticStrainNorm;
			this.plasticStrainNorm = plasticStrainNorm;
			this.yieldStress = yieldStress;
			this.hardeningRatio = hardeningRatio;
			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			this.bulkModulus = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			strainsElasticPrev = new double[6];
			strainsPlasticPrev = new double[6];
			strainsEquivalentPrev = 0;
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
			else if (totalStrains.Norm2() <= plasticStrainNorm)
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
			else
			{
				return constitutiveMatrix;
			}
		}

		private void CalculateNextStressStrainPoint()
		{
			//calculate trial variables
			var strainTrial = new double[incrementalStrains.Length];
			var stressTrial = new double[incrementalStrains.Length];
			var strainDeviatoricTrial = new double[incrementalStrains.Length];
			var stressVolumetricTrial = new double[incrementalStrains.Length];
			var stressDeviatoricTrial = new double[incrementalStrains.Length];
			var normStrainDeviatoricTrial = new double();
			for (int i = 0; i < incrementalStrains.Length; i++)
				strainTrial[i] = strainsElasticPrev[i] + incrementalStrains[i];
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					strainDeviatoricTrial[i] = strainDeviatoricTrial[i] + DeviatoricProjection[i, j] * strainTrial[j];
					if (i >= 0 && i <= 2)
						stressTrial[i] += (2 * shearModulus * DeviatoricProjection[i, j] * strainTrial[j] + bulkModulus * VolumetricProjection[i, j] * strainTrial[j]);
					else
						stressTrial[i] += 0.5 * (2 * shearModulus * DeviatoricProjection[i, j] * strainTrial[j] + bulkModulus * VolumetricProjection[i, j] * strainTrial[j]);
				}
			}
			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				for (int j = 0; j < incrementalStrains.Length; j++)
				{
					stressVolumetricTrial[i] += (double)1 / 3 * VolumetricProjection[i, j] * stressTrial[j];
					stressDeviatoricTrial[i] += DeviatoricProjection[i, j] * stressTrial[j];
				}
				normStrainDeviatoricTrial += strainDeviatoricTrial[i] * strainDeviatoricTrial[i];
			}
			normStrainDeviatoricTrial = Math.Sqrt(normStrainDeviatoricTrial);
			var J2 = GetDeviatorSecondStressInvariant(stressTrial);
			double vonMisesStress = Math.Sqrt(3 * J2);
			double vonMisesStressMinusYieldStress = vonMisesStress -
													(this.yieldStress + (this.hardeningRatio * this.strainsEquivalentPrev));

			bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;
			var stressVolumetric = new double[incrementalStrains.Length];
			var stressDeviatoric = new double[incrementalStrains.Length];

			if (materialIsInElasticRegion)
			{
				stressesNew.CopyFrom(stressTrial);
				stressVolumetric.CopyFrom(stressVolumetricTrial);
				stressDeviatoric.CopyFrom(stressDeviatoricTrial);
				strainsEquivalent = strainsEquivalentPrev;
				constitutiveMatrix = ElasticConstitutiveMatrix();
			}
			else
			{
				double deltastrainsEquivalent = vonMisesStressMinusYieldStress /
											((3 * this.shearModulus) + this.hardeningRatio);
				this.strainsEquivalent = this.strainsEquivalentPrev + deltastrainsEquivalent;

				for (int i = 0; i < incrementalStrains.Length; i++)
				{
					stressVolumetric[i] = stressVolumetricTrial[i];
					stressDeviatoric[i] = (1 - (3 * shearModulus * deltastrainsEquivalent) / vonMisesStress) * stressDeviatoricTrial[i];
					stressesNew[i] = stressVolumetric[i] + stressDeviatoric[i];
				}
				BuildConsistentTangentialConstitutiveMatrix(strainDeviatoricTrial, normStrainDeviatoricTrial, deltastrainsEquivalent, vonMisesStress);
				//this.BuildTangentialConstitutiveMatrix();
			}

			for (int i = 0; i < incrementalStrains.Length; i++)
			{
				if (i >= 0 && i <= 2)
					strainsElastic[i] = 1 / (2 * shearModulus) * stressDeviatoric[i] + 1 / (3 * bulkModulus) * stressVolumetric[i];
				else
					strainsElastic[i] = 2 * (1 / (2 * shearModulus) * stressDeviatoric[i] + 1 / (3 * bulkModulus) * stressVolumetric[i]);
				if (vonMisesStress > 0)
					strainsPlastic[i] = strainsPlasticPrev[i] + strainTrial[i] - strainsElastic[i];
				strainsNew[i] = strainsElastic[i] + strainsPlastic[i];
			}

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
			else if (totalStrains.Norm2() <= plasticStrainNorm)
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

		public NeuralNetworkModifiedMaterial3D Clone()
		{
			return new NeuralNetworkModifiedMaterial3D(neuralNetwork, YoungModulus, PoissonRatio, elasticStrainNorm, yieldStress, hardeningRatio, plasticStrainNorm, materialParameters);
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
			strainsElasticPrev.CopyFrom(strainsElastic);
			strainsPlasticPrev.CopyFrom(strainsPlastic);
			strainsEquivalentPrev = strainsEquivalent;
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

		private double[] strainsElasticPrev;

		/// <summary>
		///   The current elastic strain vector.
		/// </summary>
		private double[] strainsElastic = new double[6];

		private double[] strainsPlasticPrev;

		private double[] lastConvergedStress;

		/// <summary>
		///   The current plastic strain vector.
		/// </summary>
		private double[] strainsPlastic = new double[6];



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
