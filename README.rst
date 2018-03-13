=========================
Geometric Fitting Example
=========================

This example demonstrates geometric fitting. The initial undeformed geometry is a unit cube. The mesh is a single tri-cubic Hermite element. A number of data points on the surface of a sphere of unit radius are then generated on the top face of the cube. Geometric fitting is then used to fit the undeformed geometry. Sobolov smoothing is used. 

Command Line arguments
======================

Up to four command line arguments can be specified. They are (in order):
* number of data points
* tau Sobolv smoothing parameter
* kappa Sobolov smoothing parameter
* number of iterations


