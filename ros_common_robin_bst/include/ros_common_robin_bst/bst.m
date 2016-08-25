%%% BST B-Spline Toolbox MATLAB Interface description:
%% Initialization: result is a spline structure (see below)
% - spline = bst(degree, control_points):
%   initialization of spline with given degree for uniform knots in interval [0,1] will be
%   generated
%   degree: (positive integer) maximum degree of spline
%   control_points: (double vector) vector of control points
% - spline = bst(degree, control_points, knots)
%   initialization of a spline with given degree for parameter interval [min(knots),max(knots)]
%	The equation length(knots) - degree - 1 == length(control_points) must
%	hold.
%   degree: (positive integer) maximum degree of spline
%   control_points: (double vector) vector of control points
%   knots: (double vector) vector of knots
% - spline = bst(degree, control_points, par_start, par_end)
%   initialization of spline with given degree for uniform knots in interval [par_start,par_end] will be
%   generated
%   degree: (positive integer) maximum degree of spline
%   control_points: (double vector) vector of control points
% 	par_start: (double) start parameter of spline
% 	par_end: (double) end parameter of spline
% - spline = bst(degree, approx_struct, n_control_points)
%   initialization of a spline with given degree using an approximation
%   structure (see below) for uniform knots in interval
%   [min(approx_struct.par),max(approx_struct.par)]
%   degree: (positive integer) maximum degree of spline
%   approx_struct: (approximation structure)
%   n_control_points: (positive integer) expected number of control points
%
% spline structure fields:
% - degree: (positive integer) maximum degree of spline
% - ctrl_pts: (double vector) vector of control points
% - par_start: (double) start parameter of spline
% - par_end: (double) end parameter of spline
% - knots: (double vector) vector of knots
% - coeffs: (3D double matrix) matrix of coefficients (dimensions 1,2) for each
%           interval (dimension 3)
%
% approximation structure fields (vectors of same size):
% - par: (double vector) vector of parameter values
% - der: (positive integer vector) vector of derivative orders
% - val: (double vector) vector of values of derivative of order
%        (corresponding field in der)
%
%% Warm start:
% - values = bst(spline, parameters)
%   returns the spline value corresponding to the given parameters
%   spline: (spline structure)
%   parameters: (double vector) vector of parameters
% - values = bst(spline, parameters, derivatives)
%   returns the value of the requested derivatives corresponding to the 
%   given parameters
%   spline: (spline structure)
%   parameters: (double vector) vector of parameters
%   derivatives: (zero/positive integer vector) vector of derivative orders
%