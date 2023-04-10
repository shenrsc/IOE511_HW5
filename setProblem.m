% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that specifies the problem. Specifically, a way to compute: 
%    (1) function values; (2) gradients; and, (3) Hessians (if necessary).
%
%           Input: problem (struct), required (problem.name)
%           Output: problem (struct)
%
% Error(s): 
%       (1) if problem name not specified;
%       (2) function handles (for function, gradient, Hessian) not found
%
function [problem] = setProblem(problem)
    % check is problem name available
    if ~isfield(problem,'name')
        error('Problem name not defined!!!')
    end

    % set function handles according the the selected problem
    switch problem.name

        case 'p1'
            problem.compute_f = @p1_func;
            problem.compute_g = @p1_grad;
            problem.compute_H = @p1_H;
        
        case 'p2'
            problem.compute_f = @p2_func;
            problem.compute_g = @p2_grad;
            problem.compute_H = @p2_H;
            
        otherwise
            error('Problem not defined!!!')
    end

end


function [f, L, Fy,ceq] = p1_func(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    %L is the lagrangian function
    %Fy is the Quadratic penalty function
    f = x(1)+x(2);
    lambda = problem.lambda;
    mu = problem.mu;
    ceq = [x(1)^2+x(2)^2-2];
    L = f+lambda.'*ceq;
    % L=0;
    Fy = f+mu/2*norm(ceq)^2;
end

function [g_f, g_c, g_Fy,g_L] = p1_grad(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    %L is the lagrangian function
    %Fy is the Quadratic penalty function
    g_f = [1;1];
    lambda = problem.lambda;
    mu = problem.mu;
    g_c = [2*x(1);
            2*x(2)];
    g_L = g_f+g_c*lambda;
    g_Fy = g_f+mu*(x(1)^2+x(2)^2-2)*g_c;
end

function [H_f,H_c,H_L] = p1_H(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    %L is the lagrangian function
    %Fy is the Quadratic penalty function
    H_f = zeros(2,2);
    H_c = 2*eye(2);
    lambda = problem.lambda;
    % g_L = g_f+g_c*lambda;
    H_L = H_f+H_c*lambda;
end


function [f,L,Fy,ceq] = p2_func(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    f = exp(x(1)*x(2)*x(3)*x(4)*x(5))-1/2*(x(1)^3+x(2)^3+1)^2;
    lambda = problem.lambda;
    mu = problem.mu;
    ceq = [norm(x)^2-10;
            x(2)*x(3)-5*x(4)*x(5);
            x(1)^3+x(2)^3+1];
    L = f+lambda.'*ceq;
    % L=0;
    Fy = f+mu/2*norm(ceq)^2;
end

function [g_f, g_c, g_Fy,g_L] = p2_grad(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    %L is the lagrangian function
    %Fy is the Quadratic penalty function
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    lambda = problem.lambda;
    mu = problem.mu;
    c =  [norm(x)^2-10;
    x(2)*x(3)-5*x(4)*x(5);
    x(1)^3+x(2)^3+1];
    g_f = [x2*x3*x4*x5*exp(x1*x2*x3*x4*x5)-(x1^3+x2^3+1)*3*x1^2;
           x1*x3*x4*x5*exp(x1*x2*x3*x4*x5)-(x1^3+x2^3+1)*3*x2^2;
           x1*x2*x4*x5*exp(x1*x2*x3*x4*x5);
           x1*x2*x3*x5*exp(x1*x2*x3*x4*x5);
           x1*x2*x3*x4*exp(x1*x2*x3*x4*x5)];
    g_c = [2*x1,    0,      3*x1^2;
           2*x2,    x3,     3*x2^2;
           2*x3,    x2,     0;
           2*x4,    -5*x5,   0;
           2*x5,    -5*x4,   0];
    g_L = g_f+g_c*lambda;
    g_Fy = g_f+mu*g_c*c;
end

function [H_f,H_c,H_L] = p2_H(problem,x)
    %x is 2x1 dimension
    %f is a scale, value of ojective function
    %ceq is value of equality constraints(should be 0 if constraints satisfied)
    %c is value of inequality constraints(should be less than 0 is satisfied)
    %L is the lagrangian function
    %Fy is the Quadratic penalty function
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    lambda = problem.lambda;
    
    % H_c = zeros(5,3,5);
    H_c_1 = [2, 0,  0,  0,  0;
            0,  0,  0,  0,  0;
            6*x1,0, 0,  0,  0];

    H_c_2 = [0, 2,  0,  0,  0;
            0,  0,  1,  0,  0;
            0,  6*x2,0, 0,  0];
    
    H_c_3 = [0, 0,  2,  0,  0;
            0,  1,  0,  0,  0;
            0,  0,  0,  0,  0];
    
    H_c_4 = [0, 0,  0,  2,  0;
            0,  0,  0,  0,  -5;
            0,  0,  0,  0,  0];

    H_c_5 = [0, 0,  0,  0,  2;
            0,  0,  0,  -5, 0;
            0,  0,  0,  0,  0];
    H_c = [lambda'*H_c_1;
            lambda'*H_c_2;
            lambda'*H_c_3;
            lambda'*H_c_4;
            lambda'*H_c_5];

    H_f = zeros(5,5);
    exp_item = exp(x1*x2*x3*x4*x5);
    %row 1
    H_f(1,1) = (x2*x3*x4*x5)^2*exp_item-3*(5*x1^4+2*x1+2*x1*x2^3);
    H_f(1,2) = x3*x4*x5*exp_item+x1*x2*(x3*x4*x5)^2*exp_item-3*x1^2*(3*x2^2);
    H_f(1,3) = x2*x4*x5*exp_item+x1*x3*(x2*x4*x5)^2*exp_item;
    H_f(1,4) = x2*x3*x5*exp_item+x1*x4*(x2*x3*x5)^2*exp_item;
    H_f(1,5) = x2*x3*x4*exp_item+x1*x5*(x2*x3*x4)^2*exp_item;

    %row 2
    H_f(2,1) =  x3*x4*x5*exp_item+x1*x2*(x3*x4*x5)^2*exp_item-3*x2^2*(3*x1^2);
    H_f(2,2) =  (x1*x3*x4*x5)^2*exp_item-3*(5*x2^4+2*x2+2*x2*x1^3);
    H_f(2,3) = x1*x4*x5*exp_item+x2*x3*(x1*x4*x5)^2*exp_item;
    H_f(2,4) = x1*x3*x5*exp_item+x2*x4*(x1*x3*x5)^2*exp_item;
    H_f(2,5) =  x1*x3*x4*exp_item+x2*x5*(x1*x3*x4)^2*exp_item;

    %row 3
    H_f(3,1) = H_f(1,3);
    H_f(3,2) = H_f(2,3);
    H_f(3,3)=(x2*x1*x4*x5)^2*exp_item;
    H_f(3,4)= x2*x1*x5*exp_item+x3*x4*(x2*x1*x5)^2*exp_item;
    H_f(3,5)= x2*x1*x4*exp_item+x3*x5*(x2*x1*x4)^2*exp_item;

    %row 4
    H_f(4,1) = H_f(1,4);
    H_f(4,2) = H_f(2,4);
    H_f(4,3) = H_f(3,4);
    H_f(4,4) = (x1*x2*x3*x5)^2*exp_item;
    H_f(4,5) =  x2*x3*x1*exp_item+x4*x5*(x2*x3*x1)^2*exp_item;
    
    %row 5
    H_f(5,1) = H_f(1,5);
    H_f(5,2) = H_f(2,5);
    H_f(5,3) = H_f(3,5);
    H_f(5,4) = H_f(4,5);
    H_f(5,5) = (x2*x3*x4*x1)^2*exp_item;

    H_L = H_f+H_c;
end