function [test,flag] = isminimum(f,x,p,pprev,NIter,ConvergenceParams)
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
elseif isfield(ConvergenceParams,'RelativeStepTolerance')&& (norm(p)/norm(x))<ConvergenceParams.RelativeStepTolerance
    flag = 'Relative Step Size below tolerance';
    test=true;
elseif isfield(ConvergenceParams,'StepDifferenceTolerance')&& abs(norm(p)-norm(pprev))<ConvergenceParams.StepDifferenceTolerance
    flag = 'Step Size Difference below tolerance';
    test=true;
end
end