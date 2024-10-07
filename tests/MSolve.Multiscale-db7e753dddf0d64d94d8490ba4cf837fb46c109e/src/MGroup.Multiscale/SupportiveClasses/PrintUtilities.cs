using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;

using MGroup.LinearAlgebra.Matrices;

namespace MGroup.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    /// <summary>
    /// Supportive class with print methods for results output creation or data input
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class PrintUtilities
    {
        public static void WriteToFile(double[,] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
            writer.Dispose();
        }
        public static void WriteToFile(int[,] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
            writer.Dispose();
        }

        public static void WriteToFileMsolveInput(double[,] array, string path)
        {
            var writer = new StreamWriter(path);
            writer.Write('{');
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer.Write('{');
                for (int j = 0; j < array.GetLength(1) - 1; ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(',');
                }

                writer.Write(array[i, array.GetLength(1) - 1]);

                writer.Write('}');
                if (!(i == array.GetLength(0) - 1)) writer.Write(',');
                if (!(i == array.GetLength(0) - 1)) writer.WriteLine();
            }
            writer.Write('}');
            writer.Flush();
            writer.Dispose();

            double[,] a1 = new double[2, 2] { { 0, 0 }, { 0, 0 } };
        }

        public static void WriteToFileMsolveInput(int[,] array, string path)
        {
            var writer = new StreamWriter(path);
            writer.Write('{');
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer.Write('{');
                for (int j = 0; j < array.GetLength(1)-1; ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(',');
                }

                writer.Write(array[i, array.GetLength(1) - 1]);

                writer.Write('}');
                if(!(i==array.GetLength(0)-1)) writer.Write(',');
                if (!(i == array.GetLength(0) - 1)) writer.WriteLine();
            }
            writer.Write('}');
            writer.Flush();
            writer.Dispose();

            double[,] a1 = new double[2, 2] { { 0, 0 }, { 0, 0 } };
        }

        public static void WriteToFileDictionaryMsolveInput(Dictionary<int, int[]> subdBoundariesNodes, string generalPath, string fileAddedPath)
        {
            string boundNodesPath = generalPath + fileAddedPath;
            var writer = new StreamWriter(boundNodesPath);
            writer.WriteLine("new Dictionary<int, int[]>{");
            foreach (int subdBoundNodeID in subdBoundariesNodes.Keys)
            {
                int[] connectedsubdsIDs = subdBoundariesNodes[subdBoundNodeID];
                writer.Write("{"+$"{subdBoundNodeID}"+", new int[]{"); //TODOGer12 Console.WriteLine("o{"+$"{b1}");
                for (int j = 0; j < connectedsubdsIDs.GetLength(0) - 1; ++j)
                {
                    writer.Write(connectedsubdsIDs[j]);
                    writer.Write(',');
                }
                writer.Write(connectedsubdsIDs[connectedsubdsIDs.GetLength(0) - 1]);

                if (!(subdBoundNodeID==subdBoundariesNodes.Keys.ElementAt(subdBoundariesNodes.Count-1))) writer.WriteLine("} },");
                else writer.WriteLine("} }};");
                
            }
            
            writer.Flush();
            writer.Dispose();


        }

        public static void WriteToFile(SkylineMatrix Mat, int i1, int j1, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < 40; ++i)
            {
                for (int j = 0; j < 40; ++j)
                {
                    writer.Write(Mat[i + i1, j + j1]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
            writer.Dispose();
        }

        public static double[] ReadVector(string path)
        {
			string CultureName = Thread.CurrentThread.CurrentCulture.Name;
			CultureInfo ci = new CultureInfo(CultureName);
			ci.NumberFormat.NumberDecimalSeparator = ".";
			Thread.CurrentThread.CurrentCulture = ci;
			var reader = new StreamReader(path);
            var lines = File.ReadLines(path).Count();
            double[] data = new double[lines];
            for (int i = 0; i < lines; ++i)
            {
                data[i] = Convert.ToDouble(reader.ReadLine());

            }
            reader.Close();
            reader.Dispose();
            return data;
        }

        public static int[] ReadIntVector(string path)
        {
			string CultureName = Thread.CurrentThread.CurrentCulture.Name;
			CultureInfo ci = new CultureInfo(CultureName);
			ci.NumberFormat.NumberDecimalSeparator = ".";
			Thread.CurrentThread.CurrentCulture = ci;
			var reader = new StreamReader(path);
            var lines = File.ReadLines(path).Count();
            int[] data = new int[lines];
            for (int i = 0; i < lines; ++i)
            {
                data[i] = Convert.ToInt32(reader.ReadLine());

            }
            reader.Close();
            reader.Dispose();
            return data;
        }

        public static void WriteToFileVector(double[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();
            writer2.Dispose();

        }

        public static void WriteToFileVector(int[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();
            writer2.Dispose();

        }

        public static void WriteToFileVectorMsolveInput(int[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            writer2.Write('{');
            for (int i = 0; i < array.GetLength(0)-1; ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(',');
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }

            writer2.Write(array[array.GetLength(0) - 1]);
            writer2.Write('}');

            writer2.Flush();
            writer2.Dispose();

        }

        public static void WriteToFileVectorMsolveInput(double[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            writer2.Write('{');
            for (int i = 0; i < array.GetLength(0) - 1; ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(',');
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }

            writer2.Write(array[array.GetLength(0) - 1]);
            writer2.Write('}');

            writer2.Flush();
            writer2.Dispose();

        }

        public static void ConvertAndWriteToFileVector(double[][] array, string path3)
        {
            int length1 = array.GetLength(0);
            int length2 = array[0].GetLength(0);
            int vec_length = length1 * length2;
            double[] vector = new double[vec_length];
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array[0].GetLength(0); ++j)
                {
                    vector[length2 * i + j] = array[i][j];
                }
            }
            WriteToFileVector(vector, path3);
        }

        public static void SeparateAndWriteToFile(double[,] array, string path_A, string path_B)
        {

            int length1 = array.GetLength(0);
            int length2 = array.GetLength(1);
            double[,] array_A;
            double[,] array_B;
            array_A = new double[40, length2];
            array_B = new double[length1 - 40, length2];
            // opou 24 --> length1-40
            for (int n = 0; n < 40; n++)
            {
                for (int p = 0; p < length2; p++)
                {
                    array_A[n, p] = array[n, p];
                }
            }
            for (int n = 0; n < length1 - 40; n++)
            {
                for (int p = 0; p < length2; p++)
                {
                    array_B[n, p] = array[n + 40, p];
                }
            }

            WriteToFile(array_A, path_A);
            WriteToFile(array_B, path_B);

        }

        public static Dictionary<int, int[]> ConvertArrayToDictionary(int[] array1)
        {
            Dictionary<int, int[]> dictionary1 = new Dictionary<int, int[]>();
            int thesi = 0;
            while (thesi < array1.Length)
            {
                int ID = array1[thesi];
                thesi++;
                int[] data = new int[array1[thesi]];
                thesi++;
                for (int i1 = 0; i1 < data.Length; i1++)
                {
                    data[i1] = array1[thesi];
                    thesi++;
                }
                dictionary1.Add(ID, data);

            }

            return dictionary1;
        }
    }
}
