function [Binv,B,CholR]= FormBinv(A)
l = length(A);
B = zeros(l);
for i = 1:(l)
    B(i,i) = sum(sum(A{i}'.*A{i}));
    for j = i+1:(length(A))
        % B(i,j) = sum(sum(A{j}'.*A{i}));
        B(i,j) = A{j}(:)'*A{i}(:);
        B(j,i) = B(i,j);
    end
end
B = sparse(B);
Binv = inv(B);
if nargout>2
    CholR = chol(B);
end
end