function [CurrentLoop] = OGLP(constants)
obj_fun = constants.obj_fun; NIter = 0;   x = constants.x0;      
CurrentLoop.Iterates = x;

Binv = FormBinv(constants.A);

% Calculate residual, Jacobian and Hessian of R
[X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = 1;
[stop,CurrentLoop.ConvergenceFlag] = ismin(X.F, inf, NIter, constants);
% Main Loop
while stop == false

    [QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
    
    if length(constants.ev)<length(DFull)
        C = (DFull'-constants.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        % D = DFull(pairs(:,2));
        Dcomp = DFull;
        Dmatch = DFull(pairs(:,2));
        Dcomp(pairs(:,2))=[];
        Qmatch = QFull(:,pairs(:,2));
        Qcomp = QFull;
        Qcomp(:,pairs(:,2)) =[];
        % constants.ev = constants.ev((pairs(:,1)));
        Q = [Qmatch,Qcomp];
    else
         Q = QFull;  %D = DFull;
         Dcomp = [];
    end

    % Z = Q*diag(constants.ev)*Q';
    Z = Q*diag([constants.ev;Dcomp])*Q';
    for i = 1:length(constants.A)
        b(i,1) = trace((Z-constants.A0)'*constants.A{i});
    end
    xprev = x;
    x = Binv*b;
    p1 = x-xprev;
     p= - Binv*X.J'*X.R;

    % [x,xRGD]

    NIter = NIter + 1;
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, constants); FuncCount = FuncCount +1;

        [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, x-xprev, NIter, constants);
    % Save iterates for plotting
    CurrentLoop.Iterates = [CurrentLoop.Iterates, x];
end
CurrentLoop.FinalPoint = x;
CurrentLoop.ErrorAtFinalPoint = obj_fun(x, constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount = FuncCount;
end
