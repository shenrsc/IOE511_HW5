% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that: (1) computes the GD step; (2) updates the iterate; and, 
%                (3) computes the function and gradient at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%
function [x_new,f_new,g_new,d,alpha,ceq_new] = GD(x,f,g,ceq,problem,method,options)

    % search direction is -g
    d = -g;
    % determine step size

    % eps = 1.0/problem.mu;
    switch method.options.step_type
        case 'Constant'
            alpha = method.options.constant_step_size;
            x_new = x + alpha*d;
            [~,L,f_new,ceq_new] = problem.compute_f(problem,x_new);
            [g_f, g_c, g_new] = problem.compute_g(problem,x_new);
            % f_new = problem.compute_f(problem,x_new);
            % g_new = problem.compute_g(problem,x_new);
            
        case 'Backtracking'
            alpha = method.alpha;
            x_new = x + alpha*d;
            [~,L,f_new,ceq_new] = problem.compute_f(problem,x_new);
            [g_f, g_c, g_new] = problem.compute_g(problem,x_new);
            % f_new = problem.compute_f(problem,x_new);
            % g_new = problem.compute_g(problem,x_new);
            while f_new>f+method.options.c1*alpha*g.'*d;
                alpha = method.options.tao*alpha;
                x_new = x + alpha*d;
                [~,L,f_new,ceq_new] = problem.compute_f(problem,x_new);
                [g_f, g_c, g_new] = problem.compute_g(problem,x_new);
                % f_new = problem.compute_f(problem,x_new);
                % g_new = problem.compute_g(problem,x_new);
            end

    end
    % if norm(g_new,inf)> eps
    %     x_new = x;
    %     f_new = f;
    %     g_new = g;
    %     alpha = method.alpha;
    %     ceq_new = ceq;
    % end
    % problem.mu = problem.gama*problem.mu;
end

