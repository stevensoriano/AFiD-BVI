# AFiD-BVI
A highly parallel application for Blade-Vortex Interactions (BVI). 

AFiD-BVI is a code for highly parallel direct numerical simulations (DNS) of Blade-Vortex Interactions (BVI) built from https://github.com/PhysicsofFluids/AFiD. This code uses energy-conserving centered finite differences to spatially discretize the domain. Time marching is done through a third-order Runge-Kutta for the non-linear terms and a second-order Crank-Nicolson for the viscous terms. The domain is triply periodic.

The code uses the Immersed Boundary Method (IBM) using a moving-least-squares formulation to incorporate the body into the flow field. 

