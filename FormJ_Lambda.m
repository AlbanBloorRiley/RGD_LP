function J = FormJ_Lambda(Q,A)
%Forms the Jacobian matrix for Lambda
%
%Inputs 
%Q - the eigenvectors of the current Hamiltonian estimate
%A - The basis matrices used

l = length(A);
J = zeros(size(Q,2),l);
for k = 1:l
    J(:,k) =real(sum((Q.'*A{k}).*Q',2));
end

