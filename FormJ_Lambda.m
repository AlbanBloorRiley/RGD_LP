function J_lambda = FormJ_Lambda(Q,A)
%Forms the Jacobian matrix for Lambda
%
%Inputs 
%Q - the eigenvectors of the current Hamiltonian estimate
%A - The basis matrices used
m = size(Q,2);
if iscell(A)
     l = length(A);
    J_lambda = zeros(m,l);
%         QAQ = cell(1,l);
    for k = 1:l
%         J_lambda(:,k)=round(diag(Q'*A{k}*Q),5);
        J_lambda(:,k) =real(sum((Q.'*A{k}).*Q',2)); % Doesn't work for complex A?...
    end
    
% %     for i = 1:m
%         for k = 1:l
% %             J_lambda(i,k) = (Q(:,i)'*A{k}*Q(:,i));
%             J_lambda(:,k)=diag(QAQ{k});
%             
%         end
%         
%     end
else
    l = size(A,3);
    J_lambda = zeros(m,l);
    parfor i = 1:m
        for k = 1:l
            J_lambda(i,k) = Q(:,i)'*A(:,:,k)*Q(:,i);
        end
    end
end
end
