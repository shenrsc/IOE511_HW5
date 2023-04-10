% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that: (1) computes the GD step; (2) updates the iterate; and, 
%                (3) computes the function and gradient at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%
function [x_new,f_new,g_f_new,g_c_new,ceq_new,lambda,H_new] = SQP(x,f,g_f,g_c,ceq,H,problem,method,options)

    len = length(g_c(1,:));
    A = [H ,  g_c;
        g_c.',zeros(len)];
    b = -[g_f; ceq];
    ans = A\b;
    d = ans(1:size(H,1));
    l = ans(size(H,1)+1:end);
    x_new = x+d;
    lambda = l;
    
    [f_new,L_new,Fy_new,ceq_new] = problem.compute_f(problem,x_new);
    [g_f_new, g_c_new] = problem.compute_g(problem,x_new);
    H_new = problem.compute_H(problem,x_new);



end

