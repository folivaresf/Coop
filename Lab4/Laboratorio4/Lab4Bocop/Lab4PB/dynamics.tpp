// Function for the dynamics of the problem
// dy/dt = dynamics(y,u,z,p)

// The following are the input and output available variables 
// for the dynamics of your optimal control problem.

// Input :
// time : current time (t)
// normalized_time: t renormalized in [0,1]
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization parameters
// constants : vector of constants

// Output :
// state_dynamics : vector giving the expression of the dynamic of each state variable.

// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])

// Tdouble variables correspond to values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during optimization use standard types (double, int, ...).

#include "header_dynamics"
{
	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
	double a = constants[0];
	double m = constants[1];
	double M = constants[2];
	double b = constants[3];

	Tdouble x = state[0];
	Tdouble y = state[1];
	Tdouble z = state[2];

	Tdouble u = control[0];
	Tdouble v = control[1];

	Tdouble c = 0.0288 + (0.0613 - 0.0288)*(time - 5)*(time - 5)/25;
	Tdouble d = 0.287 + (0.767 - 0.287)*(time - 5)*(time - 5)/25;


	state_dynamics[0] = control[0] + control[1] - a*x;
	state_dynamics[1] = -b*y*z + sqrt((M - control[0])*(m + control[0]));
	state_dynamics[2] = z*(c*y - d) - control[1];
}


