using System;
using System.Collections.Generic;
using System.Diagnostics;
using System;
using System.Diagnostics;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.AnalysisWorkflow.Logging;
using MiMsolve.intermediateCodeDevelopmentClasses;
using MGroup.Constitutive.Structural.Providers;
using System.Linq;
using MGroup.MSolve.DataStructures;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Iterative;

namespace MiMsolve.multiScaleSupportiveClasses
{
    /// <summary>
    /// Newton Raphson Non-Linear Analyzer for the solution of boundary value problems that arise 
    /// from a FE2 multiscale Analysis and are related to the microstructure model.
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public class DisplacementBvpNRNLAnalyzer : IChildAnalyzer
    {
        protected readonly IAlgebraicModel algebraicModel;
        protected readonly IModel model;
        private readonly int increments;
        private int maxSteps = 1000;
        private int stepsForMatrixRebuild = 0;
        private readonly double tolerance = 1e-2;
        private double rhsNorm;
        private INonLinearParentAnalyzer parentAnalyzer = null;
        private readonly ISolver solver;
        private readonly ProblemStructural provider;
        private readonly IGlobalVector rhs ;//comment MS2:apothikevetai se afto h timh (externalLoads/increments) gia kathe subdomain kai apo ekei pernietai opou xreiasthei (p.x. subdomain.RHS)
        private IGlobalVector u ;
        private IGlobalVector du ;
        private IGlobalVector uPlusdu ;
        Dictionary<int, INode> boundaryNodes;
        Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements;
        Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements;
        //private readonly ElementEquivalentContributionsProvider equivalentContributionsAssemblers;
        private IGlobalVector globalRhs;
        private readonly Dictionary<int, LinearAnalyzerLogFactory> logFactories = new Dictionary<int, LinearAnalyzerLogFactory>();
        private readonly Dictionary<int, IAnalysisWorkflowLog[]> logs = new Dictionary<int, IAnalysisWorkflowLog[]>();
		private const string CURRENTSOLUTION = "Current solution";
		private GenericAnalyzerState currentState;

		#region necessary logs for dispalcements and forces
		public Dictionary<int, TotalLoadsDisplacementsPerIncrementLog> IncrementalLogs { get; }
            = new Dictionary<int, TotalLoadsDisplacementsPerIncrementLog>();

        

        #endregion

        public DisplacementBvpNRNLAnalyzer(IModel model, ISolver solver,
			ProblemStructural provider, int increments, IGlobalVector uInitialFreeDOFDisplacementsPerSubdomain,
            Dictionary<int, INode> boundaryNodes, Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
           /*ElementEquivalentContributionsProvider equivalentContributionsAssemblers,*/ IAlgebraicModel algebraicModel)
        {
            this.model = model;
            this.solver = solver;
            this.algebraicModel = algebraicModel;
            this.provider = provider;
            this.increments = increments;
            //this.globalRhs = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs); 
            u = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS, TODO:possibly check for compatibility elements format: u.Add(subdomain.ID, new Vector(subdomain.RHS.Length));
                                                               //this.uPlusdu = uInitialFreeDOFDisplacementsPerSubdomain; //prosthiki MS: commented out possible pass by reference
            this.boundaryNodes = boundaryNodes;
            this.initialConvergedBoundaryDisplacements = initialConvergedBoundaryDisplacements;
            this.totalBoundaryDisplacements = totalBoundaryDisplacements;
            //this.equivalentContributionsAssemblers = equivalentContributionsAssemblers;

            rhs = algebraicModel.CreateZeroVector();


            //InitializeInternalVectors();
        }

        public int SetMaxIterations
        {
            set
            {
                if (value > 0) { maxSteps = value; }
                else { throw new Exception("Iterations number cannot be negative or zero"); }
            }
        }


        public int SetIterationsForMatrixRebuild
        {
            set
            {
                if (value > 0) { stepsForMatrixRebuild = value; }
                else { throw new Exception("Iterations number for matrix rebuild cannot be negative or zero"); }
            }
        }


        private void InitializeLogs()
        {
            logs.Clear();
            foreach (int id in logFactories.Keys) logs.Add(id, logFactories[id].CreateLogs());
        }

        

        public Dictionary<int, LinearAnalyzerLogFactory> LogFactories { get { return logFactories; } }

        #region IAnalyzer Members

        
        public IAnalysisWorkflowLog[] Logs { get; set; } = new IAnalysisWorkflowLog[0];

        public IParentAnalyzer ParentAnalyzer
        {
            get => parentAnalyzer;
            set => parentAnalyzer = (INonLinearParentAnalyzer)value; //TODO: remove this cast. Now it only serves as a check
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return null; }
            set { throw new InvalidOperationException("Newton-Raphson analyzer cannot contain an embedded analyzer."); }
        }

        public void InitializeInternalVectors()
        {
            globalRhs = algebraicModel.CreateZeroVector(); 
            rhs.Clear();
            rhs.CopyFrom(solver.LinearSystem.RhsVector);
            rhs.ScaleIntoThis(1 / (double)increments);
            //u.Clear(); prosthiki MS
            if (du is null)
            {
                du = algebraicModel.CreateZeroVector();
                uPlusdu = algebraicModel.CreateZeroVector();
                uPlusdu.AddIntoThis(u);
            }
            else
            {
                du.Clear();
                uPlusdu.Clear(); 
                uPlusdu.AddIntoThis(u);
            }
           
            globalRhs.AddIntoThis(solver.LinearSystem.RhsVector);
           


            rhsNorm = provider.CalculateRhsNorm(globalRhs);
        }

        public void InitializeLinearSystemRhsVector()
        {
            IGlobalVector forces = algebraicModel.CreateZeroVector();
            algebraicModel.LinearSystem.RhsVector = forces;
        }

        private void UpdateInternalVectors()
        {
            globalRhs.Clear();
            rhs.Clear();
            rhs.CopyFrom(solver.LinearSystem.RhsVector);
            rhs.ScaleIntoThis(1 / (double)increments);
            globalRhs.AddIntoThis(solver.LinearSystem.RhsVector);

           
            rhsNorm = provider.CalculateRhsNorm(globalRhs);
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            InitializeInternalVectors();
            
        }

        private void UpdateRHS(int step)
        {
            solver.LinearSystem.RhsVector.CopyFrom(rhs);
            
        }

        public void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();//this divides the externally applied load by the number of increments 
            for (int increment = 0; increment < increments; increment++)
            {
                
                double errorNorm = 0;
                ClearIncrementalSolutionVector();//this sets du to 0
                UpdateRHS(increment);//comment MS2: apo to rhs[subdomain.ID] pernaei sto subdomain.RHS h fixed timh (externalLoads/increments) (ginetai copy kai oxi add)  AFTO thewreitai RHS sthn prwth iteration
                UpdateRHSForLinearizationContributions(increment + 1, increments);// opote edw pouu uparxei to neo rhs tou kuklou epanalhpsewn prepei na afairountai oi draseis
                double firstError = 0;
                int step = 0;
                for (step = 0; step < maxSteps; step++)
                {
                    
                    solver.Solve();
                    errorNorm = rhsNorm != 0 ? CalculateInternalRHS(increment, step, increments) / rhsNorm : 0;//comment MS2: to subdomain.RHS lamvanei thn timh nIncrement*(externalLoads/increments)-interanalRHS me xrhsh ths fixed timhs apo to rhs[subdomain.ID]
                    Debug.WriteLine("INFO NR step {0}, iteration {1}, norm: {2}", increment, step, errorNorm);
                    if (step == 0) firstError = errorNorm;
                    if (errorNorm < tolerance)
                    {
                        
                        int subdomainID = 0;
                        int c1 = 0;
                        foreach (var subdomainLogPair in IncrementalLogs)
                        {
                            if (c1 == 0) { subdomainID = subdomainLogPair.Key; }
                            TotalLoadsDisplacementsPerIncrementLog log = subdomainLogPair.Value;
                            log.LogTotalDataForIncrement(increment, step, errorNorm,
                                uPlusdu, globalRhs);
                            c1++;
                        }
                        

                        break;
                    }

                    SplitResidualForcesToSubdomains();// scatter residuals to subdomains
                    if ((step + 1) % stepsForMatrixRebuild == 0)
                    {
                        provider.Reset();
                        BuildMatrices();
                       
                    }
                }
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", step, firstError, errorNorm);
                UpdateSolution();
                //SaveMaterialState();
            }
            CopySolutionToSubdomains();
            DateTime end = DateTime.Now;

           
        }

        private double CalculateInternalRHS(int currentIncrement, int step, int totalIncrements)
        {
            globalRhs.Clear();
            

            if (currentIncrement == 0 && step == 0)
            {
                du.AddIntoThis(solver.LinearSystem.Solution);
                uPlusdu.Clear();
                uPlusdu.AddIntoThis(u);
                uPlusdu.AddIntoThis(du);
            }
            else
            {
                du.AddIntoThis(solver.LinearSystem.Solution);
                uPlusdu.Clear();
                uPlusdu.AddIntoThis(u);
                uPlusdu.AddIntoThis(du);
            }

            var boundaryNodesIds = boundaryNodes.Keys.ToList();
            var internalRhs = provider.GetRHSFromSolutionWithInitialDisplacemntsEffect(uPlusdu, boundaryNodes,
             initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, currentIncrement+1, totalIncrements);

            provider.ProcessInternalRhs(uPlusdu, internalRhs);// this does nothing

            if (parentAnalyzer != null)
            {
                var otherRhsComponents = parentAnalyzer.GetOtherRhsComponents(uPlusdu);
                internalRhs.AddIntoThis(otherRhsComponents);// this does nothing for the static problem
            }


            solver.LinearSystem.RhsVector.Clear();
            for (int j = 0; j <= currentIncrement; j++) solver.LinearSystem.RhsVector.AddIntoThis(rhs);
            solver.LinearSystem.RhsVector.SubtractIntoThis(internalRhs);
            globalRhs.AddIntoThis(solver.LinearSystem.RhsVector);
           
            double providerRHSNorm = provider.CalculateRhsNorm(globalRhs);
            return providerRHSNorm;
        }

        private void ClearIncrementalSolutionVector()
        {
            du.Clear();
            
        }

        private void SplitResidualForcesToSubdomains()
        {
            solver.LinearSystem.RhsVector.Clear();
            solver.LinearSystem.RhsVector.CopyFrom(globalRhs);
                       
        }

        private void UpdateSolution()
        {
            u.AddIntoThis(du);
        }

        public void SaveMaterialState()
        {
			//subdomainUpdaters.UpdateState(); this comment is replaced by the following two lines
			CreateState();
			provider.UpdateState(currentState);
		}

        public IGlobalVector GetConvergedSolutionVectorsOfFreeDofs()
        {
            return u;
           
        }

        public IGlobalVector GetConvergedIncrementalSolutionVectorsOfFreeDofs()
        {
            return du;
        }


        private void CopySolutionToSubdomains()
        {
                
        }

        private void ClearMaterialStresses()
        {
            //subdomainUpdaters.ResetState(); <--- se progoumenh ekdosh ginotan etsi
            
        }

        public void BuildMatrices()
        {
            if (parentAnalyzer == null)
                throw new InvalidOperationException("This Newton-Raphson non-linear analyzer has no parent.");

            parentAnalyzer.BuildMatrices();
            //solver.Initialize();
        }

        private void UpdateRHSForLinearizationContributions(int nIncrement, int increments)
        {
            globalRhs.Clear();

            
            var elementStiffnessProvider = new ElementStructuralStiffnessProvider();
            var elementContributionVecprovider = new ElementEquivalentContributionsProvider(algebraicModel, elementStiffnessProvider, boundaryNodes,
             initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, increments);
            
            IGlobalVector conttribution = algebraicModel.CreateZeroVector();
            algebraicModel.AddToGlobalVector(conttribution, elementContributionVecprovider);

            solver.LinearSystem.RhsVector.SubtractIntoThis(conttribution);
            globalRhs.AddIntoThis(solver.LinearSystem.RhsVector);

            
            rhsNorm = provider.CalculateRhsNorm(globalRhs);
            //TODOMS: possibly it is not nessesary to update globalRHS (and of course not clear it before the loop) as it is not updated in 177-UpdateRHS(increment) either.
            // and it is used when it is recalculated in CalculateInternalRHS......
        }

		GenericAnalyzerState CreateState()
		{
			currentState = new GenericAnalyzerState(this, new[]
			{
				(CURRENTSOLUTION, u),
			});

			return currentState;
		}

		IHaveState ICreateState.CreateState() => CreateState();
		GenericAnalyzerState IAnalyzer.CreateState() => CreateState();

		GenericAnalyzerState IAnalyzer.CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = false;

				u.CopyFrom(currentState.StateVectors[CURRENTSOLUTION]);

				currentState.StateVectors[CURRENTSOLUTION].CheckForCompatibility = true;
			}
		}

		public IGlobalVector CurrentAnalysisLinearSystemRhs { get => solver.LinearSystem.RhsVector; }
		public IGlobalVector CurrentAnalysisResult { get => u; }

		public IterativeStatistics AnalysisStatistics => throw new NotImplementedException();

		#endregion
	}
}
