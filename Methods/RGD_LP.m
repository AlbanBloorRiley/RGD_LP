function [CurrentLoop] = RGD_LP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;      
CurrentLoop.Iterates = x;

Binv = FormBinv(constants.A);

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; \n', NIter,X.F,norm(X.J'*X.R));
fprintf(repmat(' ',1,OutputLineLength))
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
    Fprev = X.F;
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;
    if Fprev<X.F
        X.F;
    end
    % [constants,x,Binv] = rescale(constants,x);
    fprintf(repmat('\b',1,OutputLineLength))
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; \n', NIter,X.F,norm(X.J'*X.R));

    
     [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, p, NIter, constants);
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
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