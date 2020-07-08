# multi-shooting-methods
Multi-shooting methods are useful to change a boundary value problem into an initial value problem in order to guarantee the existence and uniqueness of the solution.

The idea is that in order to change the nature of the problem it is necessary to introduce a parameter. Once that the parameter has been computed via root-finding or optimization algorithms, the initial value problem can be solved with high order numerical methods.

In my Bachelor's thesis, I have studied in details multi-shooting methods and I have dealt with different examples as the brachistochrone and the pendulum problem.

Folder description: 
# function
Functions that implement the shooting methods and are useful to determine the best parameter choice.
For example, secant method, Newton method, bisection method.
# fisica
Brachistochrone and pendulum problems solved with different shooting methods. 
