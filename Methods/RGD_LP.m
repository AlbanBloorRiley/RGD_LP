function [CurrentLoop] = RGD_LP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;
CurrentLoop.Iterates = x;
if ~isfield(constants,'doubled')
    constants.doubled = false;
end
Fprev = inf; pprev=0;
if isfield(constants, 'Verbose')
    Verbose = constants.Verbose;
else
    Verbose = false;
end
usecholesky = false;
if isfield(constants, 'BCholeskyFactor')
    R = constants.BCholeskyFactor;
    usecholesky = true;
elseif isfield(constants,'Binv')
    Binv = constants.Binv;
else
    Binv = FormBinv(constants.A);
end
% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
CurrentLoop.Error = X.F;
if Verbose
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d\n', NIter,X.F,norm(X.J'*X.R));
    fprintf(repmat(' ',1,OutputLineLength))
end
[stop,CurrentLoop.ConvergenceFlag] = isminimum(X.F,x, inf,inf, NIter, constants);
% Main Loop
while stop == false
    if usecholesky
        p = R\(R\X.J'*X.R);
    else
        p = - Binv*X.J'*X.R;
    end
    if constants.doubled
        p = p*2;
    end
    xprev = x;

    x = xprev + p;
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;


    if Fprev<X.F
        X.F;
    end
    Fprev = X.F;
    % [constants,x,Binv] = rescale(constants,x);
    if Verbose
        fprintf(repmat('\b',1,OutputLineLength))
        OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d \n', NIter,X.F,norm(X.J'*X.R));
    end
    [stop, CurrentLoop.ConvergenceFlag] = isminimum(X.F,x, p,pprev, NIter, constants);
    pprev=p;
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
    CurrentLoop.Error=[CurrentLoop.Error,X.F];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end

