using MGroup.LinearAlgebra.Commons;

namespace MGroup.Multiscale.Tests.FEMpartB
{
    public static  class ComparisonMetods //Origin: programElegxoiDdm opou eixan ginei comment out kai den htan updated apo ekdosh feat/prosthiki_allagwn 
    {        
        
        public static bool AreDisplacementsSame(double[] expectedValues,
            double[] computedValues, double tol=1E-14)
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

        public static bool AreDisplacementsSame(double[,] expectedValues,
            double[,] computedValues, double tol = 1E-14)
        {
            var comparer = new ValueComparer(tol);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {
                for (int i2 = 0; i2 < expectedValues.GetLength(0); i2++)
                {
                    if (!comparer.AreEqual(expectedValues[i1,i2], computedValues[i1,i2]))
                    {
                        return false;
                    }
                }
            }
            return true;
        }
    }
}
