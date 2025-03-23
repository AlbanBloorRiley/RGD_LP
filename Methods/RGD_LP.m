function [CurrentLoop] = RGD_LP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;      
CurrentLoop.Iterates = x;
if ~isfield(constants,'doubled')
    constants.doubled = false;
end
alpha = 1; Fprev = inf; pprev=0;
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
OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; alpha = %d \n', NIter,X.F,norm(X.J'*X.R),alpha);
fprintf(repmat(' ',1,OutputLineLength))
[stop,CurrentLoop.ConvergenceFlag] = isminimum(X.F,x, inf,inf, NIter, constants);
% Main Loop
while stop == false
    if usecholesky
     p = R\(R\X.J'*X.R);
    else
     p = - Binv*X.J'*X.R;
    end
    % p = - (X.J'*X.J)\X.J'*X.R;
    % p = - lsqminnorm(X.J'*X.J,X.J'*X.R);
    if constants.doubled
        p = p*2;
    end
    xprev = x;
    x = x+alpha*p;
    
   
    % alpha = 1;
    % while (obj_fun(xprev+ alpha*p, constants) > X.F + alpha*1e-14*dot(2*X.J'*X.R,p))
    %     alpha = alpha*0.5;
    % 
    % end
    % x = xprev + alpha*p;
    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    

    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;

    % if X.F<Fprev
    %     [X1.F,X1.R,X1.J] = obj_fun(xprev+ 1.1*alpha*p, constants); FuncCount = FuncCount +1;
    %     if X1.F<X.F
    %         alpha = alpha*1.1;
    %         x = xprev+alpha*p;
    %         X=X1;
    %     end
    % else
    %     [X1.F,X1.R,X1.J] = obj_fun(xprev+ 0.5*alpha*p, constants); FuncCount = FuncCount +1;
    %     if X1.F<Fprev
    %         alpha = max(alpha*0.5,1);
    %         x = xprev+alpha*p;
    %         X=X1;
    %     else
    %     alpha = 1;
    %     x = xprev+alpha*p;
    %     [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;
    %     end
    % end


    if Fprev<X.F
        X.F;
    end
    Fprev = X.F;
    % [constants,x,Binv] = rescale(constants,x);
    fprintf(repmat('\b',1,OutputLineLength))
OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; alpha = %d \n', NIter,X.F,norm(X.J'*X.R),alpha);
    
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

function [constants,x,Binv] = rescale(constants,x)
for i = 1:length(constants.A)
    constants.A{i} = constants.A{i}*x(i);
    x(i) = 1;
end
Binv =  FormBinv(constants.A);
end