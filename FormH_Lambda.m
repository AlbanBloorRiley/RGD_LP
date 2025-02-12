function H_lambda = FormH_Lambda(Q,A,D,ev)
    %Form S
    m = length(ev); l=length(A);    
    H_lambda = zeros(l,l,m);
    QAQ = cell(1,l);
    for i = 1:l
        QAQ{i} = Q'*A{i}*Q;
        QAQ{i} =QAQ{i}(1:m,1:m);
    end
    DD=D'-D;
    DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                H_lambda(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
end
 
% for ii = 1:10
%     m = length(ev); l=length(A);    
%     H_lambda = zeros(l,l,m);
%     for i = 1:m
%         for j=1:l
%             for k = 1:l
%                 s = 0;
%                 for t = 1:length(D) %length(Q)
%                      if abs(D(t)-D(i))<eps
%                         continue 
%                     end
%                     s = s + ((Q(:,t)'*A{k}*Q(:,i))*(Q(:,t)'*A{j}*Q(:,i)))/(D(i)-D(t));
%                 end
%                 %Calculate  H_lambda
%                 H_lambda(k,j,i) = 2*s;
%             end
%         end
%         %S = S+(D(i)-ev(i))*H_lambda(:,:,i);
%     end
%  end
