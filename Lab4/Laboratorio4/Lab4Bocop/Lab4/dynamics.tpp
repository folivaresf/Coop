// Function for the dynamics of the problem
// dy/dt = dynamics(y,u,z,p)

// Input :
// time : current time (t)
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
	// DYNAMICS FOR Lab4 PROBLEM
	// dr/dt = v
	// dv/dt = (Thrust(u) - Drag(r,v)) / m - grav(r)
	// dm/dt = -b*|u|

	double a = constants[0];
	double m = constants[1];
	double M = constants[2];

	Tdouble x = state[0];
	Tdouble J = state[1];

	state_dynamics[0] = control[0] - a*x;
	state_dynamics[1] = x*x - sqrt((M - control[0])*(m + control[0]));
}

