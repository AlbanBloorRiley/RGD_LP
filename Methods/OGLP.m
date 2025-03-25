function [CurrentLoop] = OGLP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0; l = length(x);
N = length(constants.A{1});
CurrentLoop.Iterates = x;
% pprev = 0;
Binv = FormBinv(constants.A);

if isfield(constants, 'Verbose')
    Verbose = constants.Verbose;
else
    Verbose = false;
end

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
[stop,CurrentLoop.ConvergenceFlag] = isminimum(X.F,x, inf,inf, NIter, constants);
if Verbose
    OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d\n', NIter,X.F,norm(X.J'*X.R));
    fprintf(repmat(' ',1,OutputLineLength))
end
% Main Loop
while stop == false

    [QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');

    if length(constants.ev)<N
        C = (DFull'-constants.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        % D = DFull(pairs(:,2));
        Dcomp = DFull;
        % Dmatch = DFull(pairs(:,2));
        Dcomp(pairs(:,2))=[];
        Qmatch = QFull(:,pairs(:,2));
        Qcomp = QFull;
        Qcomp(:,pairs(:,2)) =[];
        
         Z = [Qmatch,Qcomp]*diag([constants.ev;Dcomp])*[Qmatch,Qcomp]';
    else
        Z = QFull*diag(constants.ev)*QFull';
        % Q = QFull;  %D = DFull;
        % Dcomp = [];
    end

    b = zeros(l,1);
    for i = 1:l
        b(i,1) = trace((Z-constants.A0)'*constants.A{i});
    end
    xprev = x;
    x = Binv*b;
    % p1 = x-xprev;
    % p= - Binv*X.J'*X.R;

    % [x,xRGD]

    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;

    if Verbose
        fprintf(repmat('\b',1,OutputLineLength))
        OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d\n', NIter,X.F,norm(X.J'*X.R));
    end

    [stop, CurrentLoop.ConvergenceFlag] = isminimum(X.F,x, x-xprev,Inf, NIter, constants);
    % pprev=p;
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end
