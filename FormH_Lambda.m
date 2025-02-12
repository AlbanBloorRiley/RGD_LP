function H = FormH_Lambda(Q,A,D,evORm)
%Forms the Jacobian matrix for Lambda
%
%Inputs 
%Q - the eigenvectors of A(x^k), ordered so the first eigenvectors are assosiated with the prescribed eigenvalues
%A - The basis matrices 
%D - The current eigenvalues of A(x^k), ordered so the first eigenvectors are assosiated with the prescribed eigenvalues
%evORm - Either the vector of prescribed eigenvalues or the length of said vector
  
if isscalar(evORm)
    m = evORm;
else
    m = length(evORm);
end
l=length(A);
H = zeros(l,l,m);
QAQ = cell(1,l);
for i = 1:l
    QAQ{i} = Q'*A{i}*Q(:,1:m);
end
if ~isvector(D)
    D  =diag(D);
end
DD=D(1:m)'-D;
DD(abs(DD)<1e-15) = Inf;
for j=1:l
    for k = 1:l
        H(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
    end
end