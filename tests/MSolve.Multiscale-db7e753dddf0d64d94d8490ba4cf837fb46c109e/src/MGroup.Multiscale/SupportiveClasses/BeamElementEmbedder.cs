using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Embedding;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Multiscale.SupportiveClasses;

namespace MGroup.Multiscale.SupportiveClasses
{
	public class BeamElementEmbedder : BeamElementEmbedderBase
	{
		public BeamElementEmbedder(Model model, IElementType embeddedElement, IEmbeddedDOFInHostTransformationVector transformation)
			: base(embeddedElement, transformation)
		{
		}

		protected void CalculateTransformationMatrix()
		{
			base.CalculateTransformationMatrix();
			
		}
	}
}
