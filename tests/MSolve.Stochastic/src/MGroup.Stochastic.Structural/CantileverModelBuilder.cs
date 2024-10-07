using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.FEM.Structural.Line;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Stochastic.Interfaces;
using System;
using System.Collections.Generic;

namespace MGroup.Stochastic.Structural
{
    public class CantileverModelBuilder
    {
        public CantileverModelBuilder()
        {
        }

        public Model GetModel(IStochasticProcessRealizer stochasticRealizer, IStochasticDomainMapper domainMapper, int iteration)
        {


            var model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            model.NodesDictionary.Add(0, new Node(id: 0, x: 0.0, y: 0, z: 0));
            model.NodesDictionary.Add(1, new Node(id: 1, x: 0.1, y: 0, z: 0));
            model.NodesDictionary.Add(2, new Node(id: 2, x: 0.2, y: 0, z: 0));
            model.NodesDictionary.Add(3, new Node(id: 3, x: 0.3, y: 0, z: 0));
            model.NodesDictionary.Add(4, new Node(id: 4, x: 0.4, y: 0, z: 0));
            model.NodesDictionary.Add(5, new Node(id: 5, x: 0.5, y: 0, z: 0));
            model.NodesDictionary.Add(6, new Node(id: 6, x: 0.6, y: 0, z: 0));
            model.NodesDictionary.Add(7, new Node(id: 7, x: 0.7, y: 0, z: 0));
            model.NodesDictionary.Add(8, new Node(id: 8, x: 0.8, y: 0, z: 0));
            model.NodesDictionary.Add(9, new Node(id: 9, x: 0.9, y: 0, z: 0));
            model.NodesDictionary.Add(10, new Node(id: 10, x: 1, y: 0, z: 0));

            for (int i = 0; i < model.NodesDictionary.Count - 1; i++)
            {
				var element = new EulerBeam3D(new List<INode>() { model.NodesDictionary[i], model.NodesDictionary[i + 1] }, 
					stochasticRealizer.Realize(iteration, domainMapper, new[]
						{
							(model.NodesDictionary[i + 1].X + model.NodesDictionary[i].X)/2,
							(model.NodesDictionary[i + 1].Y + model.NodesDictionary[i].Y)/2,
							(model.NodesDictionary[i + 1].Z + model.NodesDictionary[i].Z)/2,
						}), 0.3)
				{
					ID = i,
					Density = 7.85,
					SectionArea = 1,
					MomentOfInertiaY = 1,
					MomentOfInertiaZ = 1,
					MomentOfInertiaPolar = 1,
				};

                model.ElementsDictionary.Add(i, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

			model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(
				new[]
				{
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationX, amount: 0d),
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationY, amount: 0d),
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.TranslationZ, amount: 0d),
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.RotationX, amount: 0d),
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.RotationY, amount: 0d),
					new NodalDisplacement(model.NodesDictionary[0], StructuralDof.RotationZ, amount: 0d)
				},
				new[]
				{
					new NodalLoad(model.NodesDictionary[10], StructuralDof.TranslationZ, amount: 10d)
				}
			));
            //model.ConnectDataStructures();

            return model;


        }
    }
}
