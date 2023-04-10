% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas

% Function that specifies the method and method specific options. 
% 
%           Input: method (struct); required (method.name)
%           Output: method (struct)
%
% Error(s): 
%   (1) if method not specified
%
function [method] = setMethod(method)

    % check is method name specified
    if ~isfield(method,'name')
        error('Method not specified!!!')
    end
    
    %check if option specified
    if ~isfield(method,'options')
        warning('option not specified!!! Setting a new one')
        method.options;
    end

    if ~isfield(method.options,'constant_step_size')
        warning('opticonstant_step_size not specified!!! Setting to default: 1')
        method.options.constant_step_size = 1;
    end

    %check if tao specified
    if ~isfield(method.options,'tao')
        warning('tao not specified!!! Setting to default: 0.5')
        method.options.tao = 0.5;
    end

    %check if c1 specified
    if ~isfield(method.options,'c1')
        warning('c1 not specified!!! Setting to default: 1e-4')
        method.options.c1 = 1e-4;
    end

end