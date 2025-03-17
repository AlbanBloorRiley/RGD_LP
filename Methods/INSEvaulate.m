function [F,R,J,H] = INSEvaulate(x,constants)
A = constants.A;
Ad = constants.A0;
for i = 1:length(x)
    Ad = Ad + x(i)*A{i};
end
if  length(A{1})<500
    [QFull,DFull] = eig(full(Ad),'vector');
elseif length(A{1})>500 && length(constants.ev)<0.5*length(A{1})&&nargout<4
    [QFull,DFull] = eigs(Ad, length(constants.ev)+10, 'smallestreal');
%     if DFull(end) < constants.ev(end)
%         [QFull,DFull] = eigs(Ad, min(3*length(constants.ev)+10,length(A{1})), 'smallestreal');
%         if DFull(end) < constants.ev(end)
%             [QFull,DFull] = eigs((Ad),length(A{1}),'smallestreal');
%         end
%     end
    DFull = diag(DFull);
else
    [QFull,DFull] = eigs((Ad),length(A{1}),'smallestreal');
    DFull = diag(DFull);
end
% if DFull(length(constants.ev)+1) < constants.ev(end)
% for i = length(constants.ev)+1:length(A{1})
%     if DFull(i) > constants.ev(end)
%         break
%     end
% end
%         C = (DFull'-constants.ev).^2;
%         pairs = matchpairs(C,100*max(max(C)));
%         D = DFull(pairs(:,2));
%         Q = QFull(:,pairs(:,2));
% else
Q = QFull(:,1:length(constants.ev)); D = DFull(1:length(constants.ev));
% end
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
    J = FormJ_Lambda(Q,A);
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
    H = FormH_Lambda(Q,A,D,constants.ev);
end
end