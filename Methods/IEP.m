function [f,varargout] = IEP(x,constants)


if length(constants.A{1})<5000
    [QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
else
    [QFull,DFull] = eigs(full(FormA(x,constants.A,constants.A0)),length(constants.A{1}));
    DFull = diag(DFull);
end

C = (DFull'-constants.ev).^2;
pairs = matchpairs(C,100*max(max(C)));
D = DFull(pairs(:,2));
Q = QFull(:,pairs(:,2));

if any(isnan(D))
    f = nan;
    for i = 1:nargout-1
        varargout{i} = nan;
    end
    return
end

f = sqrt(sum(((D-constants.ev)).^2));
if nargout>1
    varargout{1} = (D-constants.ev);
end
if nargout>2
    varargout{2} =  FormJ_Lambda(Q,constants.A);
end
if nargout>3
    varargout{3} = FormH_Lambda(QFull,constants.A,DFull,constants.ev);
end
end
