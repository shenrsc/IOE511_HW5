function [x_new,f_new,g_new,h_new,d,alpha] = BFGS(x,f,g,h,problem,method,options)
%fuction for BFGS method, get a new step

    if ~isfield(method.options,'eps')
        warning('eps for BFGS not specified!!! Setting a new one as 10e-6')
        eps = 10e-6;
    else
        eps = method.options.eps;
    end
    len = length(x);
    % H0 = eye(len);

    alpha = method.options.constant_step_size;
    %now h is inverse Hessian matrix(not hessian matrix)
    d = -h*g;
    x_new = x+alpha*d;
    f_new = problem.compute_f(problem,x_new);
    g_new = problem.compute_g(problem,x_new);
    
    while f_new>f+method.options.c1*alpha*g.'*d;
        alpha = method.options.tao*alpha;
        x_new = x + alpha*d;
        f_new = problem.compute_f(problem,x_new);
        g_new = problem.compute_g(problem,x_new);
    end
    sk = x_new-x;
    yk = g_new-g;
    if(sk.'*yk<=eps*norm(sk)*norm(yk))
        disp("skipped")
        h_new = h;
    else
        h_new = (eye(size(h))-sk*(yk.')/(sk.'*yk))*h*(eye(size(h))-yk*(sk.')/(sk.'*yk))...
                +sk*(sk.')/(sk.'*yk);
    end
end