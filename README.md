# AFiD-BVI
A highly parallel application for Blade-Vortex Interactions (BVI). 

AFiD-BVI is a code for direct numerical simulations (DNS) of Blade-Vortex Interactions (BVI) built from https://github.com/PhysicsofFluids/AFiD. The code will solve the fluid flow of an incompressible Newtonian fluid coupled with fluid-structure interactions.

The code uses energy-conserving centered finite differences to discretize the domain spatially. Time marching is done third-order Runge-Kutta for the non-linear terms and second-order Crank-Nicolson for the viscous terms. Using a moving least squares formulation, the Immersed Boundary Method (IBM) is used to incorporate the body into the flow field. The Vortex is initialised using a Lamb–Oseen vortex distribution. 

The computational domain is a triply periodic box that allows a body of any geometry to be incorporated into the flow field. The body is represented in a GTS (GNU Triangulated Surface) file format. The body can be programmed to move within the flow field. 

