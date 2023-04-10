function [x_new,f_new,g_new,h_new,d,alpha] = newton_step(x,f,g,h,problem,method,options)
%fuction for newton method, get a new step
    alpha = method.options.constant_step_size;
    d = -h\g;
    x_new = x+alpha*d;
    f_new = problem.compute_f(problem,x_new);
    g_new = problem.compute_g(problem,x_new);
    
    while f_new>f+method.options.c1*alpha*g.'*d;
        alpha = method.options.tao*alpha;
        x_new = x + alpha*d;
        f_new = problem.compute_f(problem,x_new);
        g_new = problem.compute_g(problem,x_new);
    end
    h_new = problem.compute_H(problem,x_new);
end