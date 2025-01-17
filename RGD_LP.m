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
    % alpha = 1;
    % while (obj_fun(xprev+ alpha*p, constants) > X.F + alpha*1e-14*dot(2*X.J'*X.R,p))
    %     alpha = alpha*0.5;
    % 
    % end
    % x = xprev + alpha*p;
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;
    X.F
     [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, p, NIter, constants);
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end
