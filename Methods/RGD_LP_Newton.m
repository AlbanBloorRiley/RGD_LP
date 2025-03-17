function [LPIterations,NewtonIterations,x,F,ek] =RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,varargin)
if ~isempty(varargin)
problem.doubled = varargin{1};
else
    problem.doubled = false;
end
if LPSteps>0
    LPproblem = problem;
    LPproblem.MaxIter = LPSteps;  LPproblem.StepTolerance = epsilon; 
    LPIterations = RGD_LP(LPproblem);
    x1= LPIterations.FinalPoint;
else
    LPIterations.Iterates=problem.x0;
    LPIterations.NIter = 0;
    x1=problem.x0;
end
if NewtonSteps>0
    Newtonproblem = problem;
    Newtonproblem.MaxIter = NewtonSteps;  Newtonproblem.x0=x1;
    if isfield(problem,'StepTolerance')
        Newtonproblem.StepTolerance =problem.StepTolerance;
    else
        Newtonproblem.StepTolerance = 1e-5;
    end
    NewtonIterations = Newton(Newtonproblem);
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