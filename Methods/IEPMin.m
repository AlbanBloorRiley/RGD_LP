function [F,R,J,H] = IEPsmallest(x,constants)

Ad = FormA(x,constants.A,constants.A0);
if  length(constants.A{1})<500
    [QFull,DFull] = eig(full(Ad),'vector');
elseif length(constants.A{1})>500 && length(constants.ev)<0.5*length(constants.A{1})&&nargout<4
    [QFull,DFull] = eigs(Ad, length(constants.ev), 'smallestreal');
    if  any(any(isnan(DFull)))
        [QFull,DFull] = eigs(Ad,length(constants.ev),'smallestreal','subspacedimension',min([2*length(constants.ev)+20,4*length(constants.ev),length(constants.A{1})]));
    end
    if  any(any(isnan(DFull)))
        [QFull,DFull] = eigs(Ad,length(constants.ev),'smallestreal','subspacedimension',min([4*length(constants.ev)+50,6*length(constants.ev),length(constants.A{1})]));
    end
    DFull = diag(DFull);
else
    [QFull,DFull] = eigs((Ad),length(constants.A{1}),'smallestreal');
    DFull = diag(DFull);
end

Q = QFull(:,1:length(constants.ev)); D = DFull(1:length(constants.ev));

F = sqrt(sum((D-constants.ev).^2));
if nargout>1
    R = (D-constants.ev);
end
if nargout>2

    J = FormJ_Lambda(Q,constants.A);
end
if nargout>3
    H = FormH_Lambda(Q,constants.A,D,constants.ev);
end
end