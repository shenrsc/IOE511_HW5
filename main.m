% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Script to run code

% close all figures, clear all variables from workspace and clear command
% window
close all
clear all
clc

problem.name = 'p1';
problem.x0 = [2;2];
problem.n = length(problem.x0);
x_star = [-1;-1];
problem.lambda = 1;


% problem.name = 'p2';
% problem.x0 = [-1.8;1.7;1.9;-0.8;-0.8];
% problem.n = length(problem.x0);
% x_star = [-1.71; 1.59; 1.82; -0.763; -0.763];
% problem.lambda = ones(3,1);


% set method (minimal requirement: name of method)
method.name = 'SQP';
% method.name = 'Newton';
% method.options.step_type = 'Constant';
method.options.step_type = 'Backtracking';
method.options.constant_step_size = 1;
method.options.tao = 0.5;
method.options.c1 = 1e-4;

% set options
options.term_tol = 1e-5;
options.max_iterations = 40;
options.x_star = x_star; %give the minimum


%set problem
problem.mu = 10e-4;
problem.gama = 10;

% run method and return x^* and f^*
[x,f,norm_c,k_array,f_array] = optSolverConst_Yi_Shen(problem, method, options);
semilogy(k_array,f_array+1e-26);