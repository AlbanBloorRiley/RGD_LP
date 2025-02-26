function [test,flag] = ismin(f,p,NIter,ConvergenceParams)
flag = "";
test = false;
if norm(p)<ConvergenceParams.StepTolerance
    flag = 'Step Size below tolerance';
    test=true;
elseif isnan(f)||isinf(f)
    test = true ;
    flag = 'NaN/Inf';
elseif NIter>=ConvergenceParams.MaxIter
    flag ='Max Iterations reached';
    test = true  ;
end
end