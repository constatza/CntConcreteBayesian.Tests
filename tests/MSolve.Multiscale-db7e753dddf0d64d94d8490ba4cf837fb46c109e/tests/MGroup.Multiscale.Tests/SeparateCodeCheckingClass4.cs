using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.MultiscaleAnalysis;
using MGroup.MSolve.MultiscaleAnalysis.Interfaces;
using MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Multiscale.RveTemplates;

using MiMsolve.SolutionStrategies;
using Xunit;

namespace MGroup.Multiscale.Tests.FEMpartB
{
    public class SeparateCodeCheckingClass4
    {
        [Fact]
	        public static void CheckElasticMicrostructure()
        {
            (double[] stressesCheck1, double[] stressesCheck2, double[] stressesCheck3, double[] stressesCheck4, double[] uInitialFreeDOFs_state1, double[] uInitialFreeDOFs_state2) = SeparateCodeCheckingClass4.Check05bStressIntegrationObje_Integration();

            string results_file1 = "..\\..\\InputFiles\\MicroStructureAndIntegrationTestsA\\U_sunol_micro_1.txt";
            string results_file2 = "..\\..\\InputFiles\\MicroStructureAndIntegrationTestsA\\U_sunol_micro_2.txt";
            double[] displacements1sIncrement = PrintUtilities.ReadVector(results_file1);
            double[] displacements2ndIncrement = PrintUtilities.ReadVector(results_file2);

            var ch01 = AreDisplacementsSame(displacements1sIncrement, uInitialFreeDOFs_state1);
            var ch02 = AreDisplacementsSame(displacements2ndIncrement, uInitialFreeDOFs_state2);

            Assert.True(ch01);
            Assert.True(ch02);

        }

        public static (double[], double[], double[], double[], double[], double[]) Check05bStressIntegrationObje_Integration()
        {
            //Origin: SeparateCodeCheckingClass.Check05bStressIntegration
            //modifications: tha xrhsimopoithei h nea microstructure me obje kapoia subdomainCalculations

            double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
            var material1 = new ElasticMaterial3D(youngModulus: E_disp, poissonRatio: ni_disp);
            double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateConstitutiveMatrixAndEvaluateResponse(GLVec);
            //double[] stressesCheck1 = material1.Stresses;
            double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
                material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
            DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
            GLVec = Transform_DGtr_to_GLvec(DGtr);
            material1.UpdateConstitutiveMatrixAndEvaluateResponse(GLVec);
            material1.CreateState();
            double[] stressesCheck2 = material1.Stresses;

            //VectorExtensions.AssignTotalAffinityCount();
            IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderNonLinear();
            //IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

            //IContinuumMaterial3DDefGrad
            var microstructure3 = new MicrostructureDefGrad3D<SkylineMatrix>(homogeneousRveBuilder1,
                 false, 1, new SkylineSolverPrefernce());
            //IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
            double[,] consCheck1 = new double[6, 6];
            for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }

            microstructure3.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck3 = microstructure3.Stresses;
            microstructure3.CreateState();
            var uInitialFreeDOFs_state1 = microstructure3.RetrieveDisplacementsOfFreeDofs();
            


            microstructure3.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.20, 1, 1, 0, 0, 0, 0, 0, 0 });
            double[] stressesCheck4 = microstructure3.Stresses;
            var uInitialFreeDOFs_state2 = microstructure3.RetrieveDisplacementsOfFreeDofs();
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state1.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state1.txt");
            //PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state2.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state2.txt");

            return (stressesCheck1, stressesCheck2, stressesCheck3, stressesCheck4, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);

        }

		public static (double[], double[], double[,], double[], double[]) Check_Graphene_rve_Obje_Integration()
		{
			//Origin: SeparateCodeCheckingClass4.Check05bStressIntegrationObje_Integration parontos
			//modifications: tha xrhsimopoithei o GrapheneBuilder...35...v2 gia epilush enos paradeigmatos GrapheneReinforcedRVEBuilderCHECK
			//gia elegxo twn newn domwn
			//PROSOXH gia na elegxei kai h defterh iteration u_sunol_micro_2 prepei na valoume ston graphenebuilder Addgraphenesheet xwris to bondslip.

			double E_disp = 3.5; /*Gpa*/ double ni_disp = 0.4; // stather Poisson
			var material1 = new ElasticMaterial3D(E_disp, ni_disp);
			double[,] DGtr = new double[3, 3] { { 1.10, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
			double[] GLVec = Transform_DGtr_to_GLvec(DGtr);
			material1.UpdateConstitutiveMatrixAndEvaluateResponse(GLVec);
			//double[] stressesCheck1 = material1.Stresses;
			double[] stressesCheck1 = new double[6] {material1.Stresses[0], material1.Stresses[1], material1.Stresses[2],
				material1.Stresses[3],material1.Stresses[4],material1.Stresses[5] };
			DGtr = new double[3, 3] { { 1.20, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
			GLVec = Transform_DGtr_to_GLvec(DGtr);
			material1.UpdateConstitutiveMatrixAndEvaluateResponse(GLVec);
			material1.CreateState();
			double[] stressesCheck2 = material1.Stresses;

			//VectorExtensions.AssignTotalAffinityCount();
			IRVEbuilder homogeneousRveBuilder1 = new RveGrShOne(1);
			//IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderCheckEnaHexa();

			//IContinuumMaterial3DDefGrad 
			var microstructure3 = new MicrostructureDefGrad3D<SkylineMatrix>(homogeneousRveBuilder1, false, 1, new SkylineSolverPrefernce());
			//IContinuumMaterial3DDefGrad microstructure3copyConsCheck = new Microstructure3copyConsCheckEna(homogeneousRveBuilder1);
			double[,] consCheck1 = new double[6, 6];
			for (int i1 = 0; i1 < 6; i1++) { for (int i2 = 0; i2 < 6; i2++) { consCheck1[i1, i2] = microstructure3.ConstitutiveMatrix[i1, i2]; } }
			microstructure3.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.05, 1, 1, 0, 0, 0, 0, 0, 0 });
			double[] stressesCheck3 = microstructure3.Stresses;
			microstructure3.CreateState();
			var uInitialFreeDOFs_state1 = microstructure3.RetrieveDisplacementsOfFreeDofs();

			microstructure3.UpdateConstitutiveMatrixAndEvaluateResponse(new double[9] { 1.10, 1, 1, 0, 0, 0, 0, 0, 0 });
			double[] stressesCheck4 = microstructure3.Stresses;
			var uInitialFreeDOFs_state2 = microstructure3.RetrieveDisplacementsOfFreeDofs();

			//PrintUtilities.WriteToFileVector(stressesCheck3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck3.txt");
			//PrintUtilities.WriteToFileVector(stressesCheck4, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\stressesCheck4.txt");
			//PrintUtilities.WriteToFile(consCheck1, @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\consCheck1.txt");
			//PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state1.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state1.txt");
			//PrintUtilities.WriteToFileVector(uInitialFreeDOFs_state2.CopyToArray(), @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\uInitialFreeDOFs_state2.txt");


			return (stressesCheck3, stressesCheck4, consCheck1, uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);
		}

		#region transformation methods
		public static double[] Transform_DGtr_to_GLvec(double[,] DGtr)
        {
            double[,] GL = new double[3, 3];

            //
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0;
                    for (int p = 0; p < 3; p++)
                    {
                        GL[m, n] += DGtr[m, p] * DGtr[n, p];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                GL[m, m] += -1;
            }
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    GL[m, n] = 0.5 * GL[m, n];
                }
            }

            double[] GLvec = new double[6];
            //
            for (int m = 0; m < 3; m++)
            {
                GLvec[m] = GL[m, m];
            }
            GLvec[3] = 2 * GL[0, 1];
            GLvec[4] = 2 * GL[1, 2];
            GLvec[5] = 2 * GL[2, 0];

            return GLvec;
        }

        public static bool AreDisplacementsSame(double[] expectedValues,
            double[] computedValues, double tol = 1E-14)
        {
            var comparer = new ValueComparer(tol);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
                {
                    return false;
                }
            }
            return true;
        }
        #endregion
    }
}
