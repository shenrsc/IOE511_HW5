function [x_new,f_new,g_new,h_new,d,alpha,s_all,y_all] = LBFGS(x,f,g,h0k,m,s_all,y_all,problem,method,options)
    %fuction for BFGS method, get a new step
    if ~isfield(method.options,'eps')
        warning('eps for L_BFGS not specified!!! Setting a new one as 10e-6')
        eps = 10e-6;
    else
        eps = method.options.eps;
    end
    len = length(x);
    % H0 = eye(len);

    alpha = method.options.constant_step_size;
    
    %
    [d] =  two_loop_LBFGS(m,s_all,y_all,f,g,h0k);

    %now h is inverse Hessian matrix(not hessian matrix)
    x_new = x+alpha*d;
    f_new = problem.compute_f(problem,x_new);
    g_new = problem.compute_g(problem,x_new);
    h_new = problem.compute_H(problem,x_new);
    

    %satisfiy Armijo-Wolfe condition
    while f_new>f+method.options.c1*alpha*g.'*d
        alpha = method.options.tao*alpha;
        x_new = x + alpha*d;
        f_new = problem.compute_f(problem,x_new);
        g_new = problem.compute_g(problem,x_new);
        h_new = problem.compute_H(problem,x_new);
    end

    sk = x_new-x;
    yk = g_new-g;

    size_s = size(s_all);

    %only update if sk'*yk is sufficiently positive
    if(sk.'*yk>eps*norm(sk)*norm(yk))
        if(size_s(2)<m)
            s_all = [s_all, sk];
            y_all = [y_all, yk];
        else
            s_all(:,1:end-1) = s_all(:,2:end);
            s_all(:,end) = sk;
            y_all(:,1:end-1) = y_all(:,2:end);
            y_all(:,end) = yk;
        end
    end
end