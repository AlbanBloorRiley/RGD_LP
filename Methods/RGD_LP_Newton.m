function [LPIterations,NewtonIterations,x,F,ek] =RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,varargin)
if ~isempty(varargin)
problem.doubled = varargin{1};
else
    problem.doubled = false;
end
if LPSteps>0
    problem.MaxIter = LPSteps;  problem.StepTolerance = epsilon; 
    LPIterations = RGD_LP(problem);
    x1= LPIterations.FinalPoint;
else
    LPIterations.Iterates=problem.x0;
    LPIterations.NIter = 0;
    x1=problem.x0;
end
if NewtonSteps>0
    problem.MaxIter = NewtonSteps; problem.StepTolerance = 1e-5; problem.x0=x1;
    [NewtonIterations] = Newton(problem);
    x= NewtonIterations.FinalPoint;
    ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
else
    NewtonIterations.Iterates  =x1;
    NewtonIterations.NIter = 0;
    x = x1;
    ek = norm(LPIterations.Iterates(:,end)-LPIterations.Iterates(:,end-1));
end
[F]=problem.obj_fun(x,problem);
% Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end