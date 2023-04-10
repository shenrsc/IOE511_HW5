% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that: (1) computes the GD step; (2) updates the iterate; and, 
%                (3) computes the function and gradient at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%
function [x,f,g,ceq,mu,lambda] = QP(x,f,g,ceq,problem,method,options)

    k = 0;
    norm_g = norm(g,inf);
    eps = 1.0/problem.mu;
    method.alpha = method.options.constant_step_size; %reset to value
    while  norm_g>eps && k<100
        % take step according to a chosen method

        
        [x_new,f_new,g_new,d,alpha,ceq_new] = GD(x,f,g,ceq,problem,method,options);
        % update old and new function values
        x_old = x; f_old = f; g_old = g; norm_g_old = norm_g;
        x = x_new; f = f_new; g = g_new; norm_g = norm(g,inf);
        method.alpha = alpha; ceq = ceq_new; norm_c = norm(ceq,inf);

        % increment iteration counter
        k = k + 1;
    end
    
    mu = problem.gama*problem.mu;
    [g_f, g_c, g_new, g_L] = problem.compute_g(problem,x);
    lambda = -g_c\g_f;

end

