function Ad = FormA(d,A,varargin)
%Forms Current estimate for the Desired Matrix
%
%Inputs
%d - current estimate for the parameters
%A and A0 - The basis matrices used (A0 is optional)
if nargin < 3
    Ad = d(1)*A{1};
elseif nargin == 3
    Ad = varargin{1}+d(1)*A{1};
else
    error("incorect number of inputs")
end

for i = 2:length(d)
    Ad = Ad + d(i)*A{i};
end

end
