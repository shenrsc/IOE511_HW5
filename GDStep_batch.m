% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that: (1) computes the GD step; (2) updates the iterate; and, 
%                (3) computes the function and gradient at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%
function [x_new,f_new,g_new,d,alpha] = GDStep_batch(x,f,g,h,x_train,y_train,problem,method,options)

    % search direction is -g
    d = -g;

    % determine step size
    switch method.options.step_type
        case 'Constant'
            alpha = method.options.constant_step_size;
            x_new = x + alpha*d;
            f_new = problem.compute_f(problem,x_new,x_train,y_train);
            g_new = problem.compute_g(problem,x_new,x_train,y_train);
            
        case 'Backtracking'
            alpha = method.alpha;
            x_new = x + alpha*d;
            f_new = problem.compute_f(problem,x_new,x_train,y_train);
            g_new = problem.compute_g(problem,x_new,x_train,y_train);
            while f_new>f+method.options.c1*alpha*g.'*d;
                alpha = method.options.tao*alpha;
                x_new = x + alpha*d;
                f_new = problem.compute_f(problem,x_new,x_train,y_train);
                g_new = problem.compute_g(problem,x_new,x_train,y_train);
            end
            
        otherwise
            error('method not defined!!!')

    end
end

