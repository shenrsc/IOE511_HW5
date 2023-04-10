% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that runs a chosen algorithm on a chosen problem
%           Inputs: problem, method, options (structs)
%           Outputs: final iterate (x) and final function value (f)
function [x,f,norm_c,k_array,f_array] = optSolverConst_Yi_Shen(problem,method,options)

    % set problem, method and options
    [problem] = setProblem(problem);
    [method] = setMethod(method);
    [options] = setOptions(options);

    % compute initial function/gradient/Hessian
    x = problem.x0;
    [f,~,Fy,ceq] = problem.compute_f(problem,x);
    [g_f, g_c, g,g_L] = problem.compute_g(problem,x);
    problem.lambda = -g_c\g_f;
    % eps = 1.0/problem.mu;
    % while norm(g,inf)>eps
    %     problem.mu
    % end
    % f = problem.compute_f(problem,x);
    % g = problem.compute_g(problem,x);
    [~,~,H] = problem.compute_H(problem,x);
    norm_c = norm(ceq,inf);
    norm_g = norm(g_L,inf);
    method.alpha = method.options.constant_step_size;

    % set initial iteration counter
    k = 0;
    norm0 = norm_g;
    normc0 = norm_c;

    k_array = k;
    f_array = f;

    while k<options.max_iterations && (norm_g>options.term_tol*max(1,norm0) || norm_c>options.term_tol*max(1,normc0) )
        
        % take step according to a chosen method
        switch method.name
            case 'QP'
                [x_new,Fy_new,g_new,ceq_new,mu,lambda] = QP(x,Fy,g,ceq,problem,method,options);
                problem.mu = mu;
                problem.lambda = lambda;
                %recalculate Fy_new and g_new so that distory next sub-optimization 
                [f,~,Fy_new,ceq_new] = problem.compute_f(problem,x_new);
                [g_f, g_c, g_new,g_L] = problem.compute_g(problem,x_new);

                x_old = x; f_old = Fy; g_old = g; norm_g_old = norm_g;
                x = x_new; Fy = Fy_new; g = g_new; norm_g = norm(g_L,inf);
                ceq = ceq_new; norm_c = norm(ceq,inf);
                
                
                
            % case 'Newton'
            %     [x_new,f_new,g_new,d,alpha] = newton_step(x,f,g,h,problem,method,options);
            case 'SQP'
                [x_new,f_new,g_f_new,g_c_new,ceq_new,lambda,H_new] = SQP(x,f,g_f,g_c,ceq,H,problem,method,options);
                %x_old = x; f_old = Fy; g_old = g; norm_g_old = norm_g;
                x = x_new; f = f_new; g_f = g_f_new; g_c = g_c_new; norm_g = norm(g_L,inf);
                ceq = ceq_new; norm_c = norm(ceq,inf);
                problem.lambda = lambda;
            otherwise
                
                error('Method not implemented yet!')
                
        end
        
        % update old and new function values
        

        % increment iteration counter
        k = k + 1;


        %store the list of results
        k_array = [k_array,k];
        f_array = [f_array,f];
        
        
    end
    [v_star,L_star,~,ceq_star] = problem.compute_f(problem,options.x_star);
    f_array = abs(f_array-v_star);
end