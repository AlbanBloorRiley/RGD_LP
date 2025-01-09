function [Binv,B]= FormBinv(A)
B = zeros(length(A));
for i = 1:(length(A))
    for j = 1:(length(A))
        B(i,j) = sum(sum(A{j}'.*A{i}));
    end
end
Binv = inv(B);
end