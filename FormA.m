function Ad = FormA(d,A,A0)
%Forms Current estimate for the Desired Matrix
%
%Inputs
%d - current estimate for the parameters
%A and A0 - The basis matrices used

Ad = A0+d(1)*A{1};
if iscell(A)
    for i = 2:length(d)
        Ad = Ad + d(i)*A{i};
    end
else
    for i = 2:length(d)
        Ad = Ad + d(i)*A(:,:,i);
    end
end

end
