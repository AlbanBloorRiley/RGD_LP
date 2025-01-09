function [CurrentLoop] = Newton(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;      
CurrentLoop.Iterates = x;

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J,X.H] = obj_fun(x, constants); FuncCount = 1;
    X.S = zeros(length(x));
    for i = 1:length(X.R)
        X.S = X.S + X.H(:,:,i)*X.R(i);
    end
    Hess = X.J'*X.J+X.S;
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F, inf, NIter, constants);
% Main Loop
while stop == false

    % p = - Hess\(X.J'*X.R);
    %  p1 = -mldivide(Hess,X.J'*X.R);
    % p = -lsqminnorm(Hess,X.J'*X.R);
    p = -constants.Solver(Hess,X.J'*X.R);
    x = x+p;
    
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J,X.H] = obj_fun(x, constants); FuncCount = 1;
    X.S = zeros(length(x));
    for i = 1:length(X.R)
        X.S = X.S + X.H(:,:,i)*X.R(i);
    end
    Hess = X.J'*X.J+X.S;

     [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, p, NIter, constants);
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end


