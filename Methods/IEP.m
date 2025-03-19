function [f,varargout] = IEP(x,constants,varargin)
if ~isempty(varargin)
    match = varargin{1};
else
    match = false;
end

if length(constants.A{1})<1000
    [QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
elseif length(constants.ev)<0.5*length(constants.A{1}) &&nargout<4
    [QFull,DFull] = eigs(FormA(x,constants.A,constants.A0),length(constants.ev),'smallestreal');
    if  any(isnan(DFull))
        [QFull,DFull] = eigs(FormA(x,constants.A,constants.A0),length(constants.ev),'smallestreal','subspacedimension',2*min([length(constants.ev)+10,2*length(constants.ev),length(constants.A{1})]));
    end
    if any(isnan(DFull))
        [QFull,DFull] = eigs(full(FormA(x,constants.A,constants.A0)),length(constants.A{1}),'smallestreal');
    end
    DFull = diag(DFull);
else
    [QFull,DFull] = eigs(full(FormA(x,constants.A,constants.A0)),length(constants.A{1}),'smallestreal');
    DFull = diag(DFull);
end


if length(constants.ev)<length(DFull)
    if match
        C = (DFull'-constants.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        D = DFull(pairs(:,2));
        Q = QFull(:,pairs(:,2));
        % constants.ev = constants.ev((pairs(:,1)));
    else
        D = DFull(1:length(constants.ev)); Q = QFull(:,1:length(constants.ev));
    end
else
    D = DFull; Q = QFull;
end
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
