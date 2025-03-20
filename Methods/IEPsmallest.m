function [F,R,J,H] = IEPsmallest(x,constants,varargin)
if ~isempty(varargin)
    match = varargin{1};
else
    match = false;
end

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
    elseif match &&DFull(end) < constants.ev(end)
        [QFull,DFull] = eigs(Ad, min(3*length(constants.ev)+10,length(constants.A{1})), 'smallestreal');

        if DFull(end) < constants.ev(end)
            [QFull,DFull] = eigs((Ad),length(constants.A{1}),'smallestreal');
        end
    end
    DFull = diag(DFull);
else
    [QFull,DFull] = eigs((Ad),length(constants.A{1}),'smallestreal');
    DFull = diag(DFull);
end
if match && length(DFull)>length(constants.ev)&&DFull(length(constants.ev)+1) < constants.ev(end)
for i = length(constants.ev)+1:length(constants.A{1})
    if DFull(i) > constants.ev(end)
        break
    end
end
        C = (DFull'-constants.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        D = DFull(pairs(:,2));
        Q = QFull(:,pairs(:,2));
else
Q = QFull(:,1:length(constants.ev)); D = DFull(1:length(constants.ev));
end
F = sqrt(sum((D-constants.ev).^2));
if nargout>1
    R = (D-constants.ev);
end
if nargout>2
    % l = length(A);
    % m = size(Q,2);
    % J = zeros(m,l);
    % for k = 1:l
    %     J(:,k) =real(sum((Q.'*A{k}).*Q',2));
    % end
    J = FormJ_Lambda(Q,constants.A);
end
if nargout>3
    % m = length(constants.ev); l=length(A);
    % H = zeros(l,l,m);
    % QAQ = cell(1,l);
    % for i = 1:l
    %     QAQ{i} = QFull'*A{i}*Q;
    %     % QAQ{i} = Q'*A{i}*Q;
    %     % QAQ{i} =QAQ{i}(1:m,1:m);
    % end
    % DD=D'-DFull;
    % % DD=D'-D;
    % DD(abs(DD)<1e-15) = Inf;
    % for j=1:l
    %     for k = 1:l
    %         H(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
    %     end
    % end
    H = FormH_Lambda(Q,constants.A,D,constants.ev);
end
end