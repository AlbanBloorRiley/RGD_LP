function Ad = FormA(x,A,varargin)
%Forms Current estimate for the Desired Matrix
%
%Inputs
%d - current estimate for the parameters
%A and A0 - The basis matrices used (A0 is optional)
if nargin>3
    error("incorect number of inputs")
end
if iscell(A)
I = [];
J = [];
V = [];
n = 0;
for k = 1:length(A)
    [i,j,v] = find(A{k});
    p = n + numel(i);
    m = numel(I);
    if p > m
        m = max(p,2*m);
        I(m) = 0;
        J(m) = 0;
        V(m) = 0;
    end
    idx = (n+1:p);
    I(idx) = i;
    J(idx) = j;
    V(idx) = v*x(k);
    n = p;
end
idx = (n+1:numel(I));
I(idx) = [];
J(idx) = [];
V(idx) = [];
[m,n] = size(A{1});
Ad = sparse(I,J,V,m,n);

if nargin == 3
    Ad = varargin{1}+Ad;
end
elseif isstruct(A)

lidx =1; v = zeros(sum(A.NumV),1);
    for k = 1:length(x)
        uidx = lidx+A.NumV(k)-1;
        v(lidx:uidx) = A.v(lidx:uidx)*x(k);
        lidx =uidx+1;
    end
Ad = sparse(A.i,A.j,v,A.m,A.n);

else
    error("A must be a cell array of basis matrices, or a structure used to construct a sparse matrix")
end

end
