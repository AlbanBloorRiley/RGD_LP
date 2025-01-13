function [CurrentLoop] = RGD_LP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;      
CurrentLoop.Iterates = x;

Binv = FormBinv(constants.A);

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F, inf, NIter, constants);
% Main Loop
while stop == false

    p = - Binv*X.J'*X.R;
    if constants.doubled
        p = p*2;
    end
    % xprev = x;
    x = x+p;
    
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;

     [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, p, NIter, constants);
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end
