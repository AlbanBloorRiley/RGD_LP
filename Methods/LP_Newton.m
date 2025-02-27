function [LPIterations,NewtonIterations,x,F,ek,Eigs] =LP_Newton(problem,epsilon,LPSteps,NewtonSteps)
if LPSteps>0
    problem.MaxIter = LPSteps; problem.StepTolerance = epsilon;
    LPIterations = OGLP(problem);
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
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end