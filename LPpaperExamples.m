warning('off','MATLAB:nearlySingularMatrix');
%% Example 1
feature('jit', 'off')
clear all
clf
N = 5;
vec = zeros(1,N);
vec(2) = -1;
% problem.A0 =[ 0 -1 0 0 0; -1 0 -1 0 0;0 -1 0 -1 0;0 0 -1 0 -1;0 0 0 -1 0];
problem.A0 = toeplitz(vec);
for k = 1:N
     D= zeros(N,1); D(k) = 4;
    problem.A{k} =  diag(D);
end
problem.x0 = [0.63160, 0.23780, 0.90920, 0.98660, 0.50070]';
% problem.x0 = randn([N,1]);
problem.ev =  [1, 1, 2, 3, 4]';
% problem.ev = 1:N;
problem.obj_fun = @IEP;
%% Example 1 bigger
feature('jit', 'off')
clear all
clf
N = 25;
vec = zeros(1,N);
vec(2) = -1;
% problem.A0 =[ 0 -1 0 0 0; -1 0 -1 0 0;0 -1 0 -1 0;0 0 -1 0 -1;0 0 0 -1 0];
problem.A0 = toeplitz(vec);
for k = 1:N
     D= zeros(N,1); D(k) = 4;
    problem.A{k} =  diag(D);
end
problem.x0 = ones(N,1);
% % problem.ev =  [1, 1, 2, 3, 4]';
problem.ev = (1:N)';
problem.obj_fun = @IEP;

%%
epsilon = 0.00; doubled = false; Opt = [];
% [x,ek,Eigs,F,LPIterations,NewtonIterations] =LP_Newton(problem,epsilon,LPSteps,Opt,doubled);
% F,ek 
% Opt.Regularisation = @(NIter,X)max(max(X.J))/max(max(diag(diag(X.J'*X.J))))/(NIter);
N = 200;
repeats = 2;
times = zeros(1,N);
timesOG = zeros(repeats,N);
timesNew = zeros(repeats,N);

NewtonIterOG = inf;
NewtonIterNew = inf;
NIterChangeOG = [];
NIterChangeNew = [];
for  i = 1:repeats
    for LPSteps = 1:N
        tic
        [LPIterations,NewtonIterations]=Old_LP_Newton(problem,0,LPSteps,Inf);
        timesOG(i,LPSteps) = toc;
        if i==repeats
            if NewtonIterOG > NewtonIterations.NIter
                NewtonIterOG = NewtonIterations.NIter;
                NIterChangeOG(end+1) = LPSteps;
            end
        end
    end
    for  LPSteps = 1:N
        tic
        [LPIterations,NewtonIterations]=RGD_LP_Newton(problem,0,LPSteps,Inf,false);
        timesNew(i,LPSteps) =  toc;
        if i==repeats
            if NewtonIterNew > NewtonIterations.NIter
                NewtonIterNew = NewtonIterations.NIter;
                NIterChangeNew(end+1) = LPSteps;
            end
        end
    end
    timesNew(LPSteps) = timesNew(LPSteps)+ toc;
end
%%

clf
% sum(sum(timesOG))
% TFOG = ~isoutlier(timesOG);
TFOG = ~isoutlier(timesOG',ThresholdFactor=2)';

TFOGrows = sum(TFOG,1);
OGdata = timesOG;
OGdata(~TFOG) = 0;
% sum(sum(OGdata))
OGdata = sum(OGdata,1)./TFOGrows;
OGdata = smoothdata(OGdata,'movmean',1);
%  sum(sum(timesNew))
% TFNew = ~isoutlier(timesNew);
TFNew = ~isoutlier(timesNew',ThresholdFactor=2)';
TFNewrows =  sum(TFNew,1);
Newdata = timesNew;
Newdata(~TFNew) = 0;
% sum(sum(Newdata))
Newdata = sum(Newdata,1)./TFNewrows;
Newdata = smoothdata(Newdata,'movmean',1);
hold on
plot(1:N,Newdata,LineWidth=2,Color=	"#0072BD")
plot(NIterChangeNew,Newdata(NIterChangeNew),'x',MarkerSize = 5,Color="black")
plot(1:N,OGdata,LineWidth=2,Color="#D95319")
plot(NIterChangeOG,OGdata(NIterChangeOG),'x',MarkerSize = 5,Color="black")
hold off
xlabel("The number of LP iterations")
ylabel("CPU time (in seconds)")

%
clf

hold on
plot(1:N,sum(timesNew,1),LineWidth=2,Color=	"#0072BD")
% plot(NIterChangeNew,timesNew(NIterChangeNew),'x',MarkerSize = 3,Color=	"#0072BD")
plot(1:N,sum(timesOG,1),LineWidth=2,Color="#D95319")
% plot(NIterChangeOG,timesOG(NIterChangeOG),'x',MarkerSize = 3,Color="#D95319")
hold off

%% Example 2
clear all
clf
problem.A0 =zeros(20);
for k = 1:20
     D= zeros(20,1); D(k) = 1;
    problem.A{k} =  toeplitz(D);
end
problem.x0 = [1.1650, 0.6268, 0.0751, 0.3516, -0.6965, 1.6961, 0.0591, 1.7971, ...
    0.2641, 0.8717, -1.4462, -0.7012, 1.2460, -0.6390, 0.5773, -0.3600,...
    -0.1356, -1.3493, -1.2704, 0.9845]';
problem.ev =  [-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5]';

epsilon = 0.00; LPSteps = 0;  doubled = false;  problem.obj_fun = @IEP;
problem.StepTolerence = 1e-8;
repeats = 1;
problem.Solver = @lsqminnorm;
tic
for i = 1:repeats
[~,OnlyNewtonIterations] =RGD_LP_Newton(problem,epsilon,LPSteps,500,doubled);
end
NewtonTime = toc/repeats;
ek1 = vecnorm(OnlyNewtonIterations.Iterates(:,2:end)-OnlyNewtonIterations.Iterates(:,1:end-1),2,1);
epsilon = 0.01; LPSteps = inf; 
tic
for i = 1:repeats
[LPIterations01,NewtonIterations01] = RGD_LP_Newton(problem,epsilon,LPSteps,Inf,doubled);
end
New01Time = toc/repeats;
Iterates01 = [LPIterations01.Iterates,NewtonIterations01.Iterates(:,2:end)];
ek2 = vecnorm(Iterates01(:,2:end)-Iterates01(:,1:end-1),2,1);
epsilon = 0.001; 
tic
for i = 1:repeats
[LPIterations001,NewtonIterations001] =RGD_LP_Newton(problem,epsilon,LPSteps,Inf,doubled);
end
New001Time  = toc/repeats;
Iterates001 = [LPIterations001.Iterates,NewtonIterations001.Iterates(:,2:end)];
ek3 = vecnorm(Iterates001(:,2:end)-Iterates001(:,1:end-1),2,1);

epsilon = 0.01;  
tic
for i = 1:repeats
[oldLPIterations01,oldNewtonIterations01] =Old_LP_Newton(problem,epsilon,LPSteps,Inf);
end
Old01Time  = toc/repeats;
oldIterates01 = [oldLPIterations01.Iterates,oldNewtonIterations01.Iterates(:,2:end)];
ek4 = vecnorm(oldIterates01(:,2:end)-oldIterates01(:,1:end-1),2,1);
epsilon = 0.001; 
tic
for i = 1:repeats
[oldLPIterations001,oldNewtonIterations001,x,F,ek,Eigs] =Old_LP_Newton(problem,epsilon,LPSteps,Inf);
end
Old001Time  = toc/repeats;
oldIterates001 = [oldLPIterations001.Iterates,oldNewtonIterations001.Iterates(:,2:end)];
ek5 = vecnorm(oldIterates001(:,2:end)-oldIterates001(:,1:end-1),2,1);



["Algorithm", "No. Iterations","CPU Time (seconds)";...
  "Newton", OnlyNewtonIterations.NIter, NewtonTime ;...    
  "Old LP-Newton, \epsilon = 0.01", string(oldLPIterations01.NIter+"+"+oldNewtonIterations01.NIter),Old01Time;...
  "Old LP-Newton, \epsilon = 0.001", string(oldLPIterations001.NIter+"+"+oldNewtonIterations001.NIter),Old001Time;...
  "New LP-Newton, \epsilon = 0.01", string(LPIterations01.NIter+"+"+NewtonIterations01.NIter),New01Time;...
  "New LP-Newton, \epsilon = 0.001", string(LPIterations001.NIter+"+"+NewtonIterations001.NIter),New001Time
]


f = figure(1);
semilogy(1:length(ek1),ek1,1:length(ek4),ek4,1:length(ek5),ek5,1:length(ek2),ek2,1:length(ek3),ek3, LineWidth=1.5)
xlabel("The number of iterations")
ylabel("log_{10}(error)")
legend({"Newton Method","Original LP-Newton"+newline+ "(\epsilon = 0.01)",...
    "Original LP-Newton"+newline+ "(\epsilon = 0.001)","New LP-Newton"+newline+...
    "(\epsilon = 0.01)","New LP-Newton"+newline+ "(\epsilon = 0.001)"},Location='ne')
ylim([1e-10,1e10])
f.Units = 'centimeters';
f.Position = [-50 10 20 14];
linestyleorder('default')
print(f, 'Example2.eps', '-depsc')
%



%% Example 3 
clear problem
problem.A0 =zeros(16);
T = [4 -1 0 0 ;-1 4 -1 0; 0 -1 4 -1; 0 0 -1 4];
E = -eye(4);
A = [T E zeros(4) zeros(4); E T E zeros(4); zeros(4) E T E; zeros(4) zeros(4) E T];
L = chol(A,'lower');
for k = 1:16
     D= zeros(16,1); D(k) = 1;
    problem.A{k} =  L'*diag(D)*L;
end
problem.x0 = [1.5578, -2.4443, -1.0982, 1.1226, 0.5817,-0.2714, 0.4142, -0.9778,...
    -1.0215, 0.3177,1.5161, 0.7494, -0.5077, 0.8853, -0.2481, -0.7262]';
problem.ev =  [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]';

%RGD_LP_Newton(problem,epsilon,LPSteps,Inf,false);

repeats = 100;
LPSteps = 0; doubled = false; problem.obj_fun = @IEP; problem.Solver = @mldivide;
tic
for i = 1:repeats
[~,OnlyNewtonIterations] =RGD_LP_Newton(problem,Inf,LPSteps,40,doubled);
end
Newtontime = toc/repeats;
ek1 = vecnorm(OnlyNewtonIterations.Iterates(:,2:end)-OnlyNewtonIterations.Iterates(:,1:end-1),2,1);
epsilon = 0.01; LPSteps = 99;%Opt.MaxIter =Opt.MaxIter-LPSteps;
tic
for i = 1:repeats
[LPIterations,NewtonIterations] =RGD_LP_Newton(problem,epsilon,LPSteps,Inf,doubled);
end
LPNewtontime = toc/repeats;
Iterates = [LPIterations.Iterates,NewtonIterations.Iterates(:,2:end)];
ek2 = vecnorm(Iterates(:,2:end)-Iterates(:,1:end-1),2,1);

tic
for i = 1:repeats
[OGLPIterations,OGNewtonIterations] =Old_LP_Newton(problem,epsilon,LPSteps,Inf);
end
OGLPNewtontime = toc/repeats;
OGIterates = [OGLPIterations.Iterates,OGNewtonIterations.Iterates(:,2:end)];
ek3 = vecnorm(OGIterates(:,2:end)-OGIterates(:,1:end-1),2,1);

% epsilon=0.001;
% tic
% [LPIterations001,NewtonIterations001] =RGD_LP_Newton(problem,epsilon,LPSteps,Inf,doubled);
% LPNew001tontime = toc;
% Iterates = [LPIterations001.Iterates,NewtonIterations001.Iterates(:,2:end)];
% ek4 = vecnorm(Iterates(:,2:end)-Iterates(:,1:end-1),2,1);


["Algorithm", "No. Iterations","CPU Time (seconds)";...
  "Newton", OnlyNewtonIterations.NIter, Newtontime ;...
  "New LP-Newton", string(LPIterations.NIter+"+"+NewtonIterations.NIter),LPNewtontime;...
  "OG LP-Newton", string(OGLPIterations.NIter+"+"+OGNewtonIterations.NIter),OGLPNewtontime]
  % "LP-Newton", string(LPIterations001.NIter+"+"+NewtonIterations001.NIter),LPNew001tontime]
%
clf
f=figure(1);
semilogy(1:length(ek1),ek1,1:length(ek2),ek2,1:length(ek3),ek3, LineWidth=1.5)
legend({"Newton Method","LP-Newton New","LP-Newton Original"},Location='se')
% ylim([1e-10,1e40])
xlabel("The number of iterations")
ylabel("log_{10}(error)")

f.Units = 'centimeters';
f.Position = [-50 10 20 14];
linestyleorder('default')
% print(f, 'Example3.eps', '-depsc')
%% Mn12 Example
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;  
B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;
B22 = (6.8*10^-4)*meV;
x0 = round([B20;B40;B44;B22;-9.1192e+05],2,'significant'); % 
Sys.S = 10; 
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = ham(Sys,[0,0,0]);  [~,E]=eig(H);
EE = diag(E);  Exp.ev=EE-EE(1);

% A{5} = speye(length(EE));
%The Stevens Operators
A0 = sparse(21,21);
A{1} = stev(10,[2,0]);
A{2} = stev(10,[4,0]);
A{3} = stev(10,[4,4]);
A{4} = stev(10,[2,2]);
A{5} = eye(21);
constants.A = A;
constants.ev = EE;
constants.A0 = A0;
obj_fun = @INSEvaulate;
problem = constants; problem.obj_fun = obj_fun; problem.x0 = x0;
[B,Binv] = FormB(problem);
Opt = struct('NDeflations',1,'Verbose',false,'constants',problem,'Method',"UGradientDescent","ScalingMatrix",Binv,...
    'ObjectiveTolerance',0,'StepTolerance',1e-10,'MaxIter',1000000,'Linesearch','No');
tic
newIterations = DMin(obj_fun,problem.x0,Opt);
LPnewTime = toc;
ekLPnew = vecnorm(newIterations.Iterates(:,2:end)-newIterations.Iterates(:,1:end-1),2,1);
%
Opt = struct('NDeflations',1,'Verbose',false,'constants',problem,'Method',"LP","ScalingMatrix",Binv,...
    'ObjectiveTolerance',0,'StepTolerance',1e-10,'MaxIter',1000,'Linesearch','No');
tic
oldIterations = DMin(obj_fun,problem.x0,Opt);
LPoldTime = toc;
ekLPold = vecnorm(oldIterations.Iterates(:,2:end)-oldIterations.Iterates(:,1:end-1),2,1);

epsilon = 0.01; LPSteps = inf; doubled = false;
tic
[~,~,~,~,newLPIterations,newNewtonIterations] =LP_Newton(problem,epsilon,LPSteps,Opt,doubled,0);
LPNnewTime =toc;
newIterates = [newLPIterations.Iterates,newNewtonIterations.Iterates(:,2:end)];
ekLPNnew = vecnorm(newIterates(:,2:end)-newIterates(:,1:end-1),2,1);

tic
[~,~,~,~,oldLPIterations,oldNewtonIterations] =LP_Newton(problem,epsilon,LPSteps,Opt,doubled,1);
LPNoldTime =toc;
oldIterates = [oldLPIterations.Iterates,oldNewtonIterations.Iterates(:,2:end)];
ekLPNold = vecnorm(oldIterates(:,2:end)-oldIterates(:,1:end-1),2,1);
%

[LPnewTime;LPoldTime;LPNnewTime;LPNoldTime]
clf
% hold on
% semilogy(ekLPnew,LineWidth=1.2)
% semilogy(ekLPold,LineWidth=1.2)
% semilogy(ekLPNnew,LineWidth=1.2)
% semilogy(ekLPNold,LineWidth=1.2)
% set(gca, 'YScale', 'log')
% hold off

semilogy(vecnorm(oldIterations.Iterates-oldIterations.DeflatedPoint,2,1),LineWidth=1.2)
hold on
semilogy(vecnorm(newIterations.Iterates-newIterations.DeflatedPoint,2,1),LineWidth=1.2)
semilogy(vecnorm([newLPIterations.Iterates,newNewtonIterations.Iterates(:,2:end)]-newNewtonIterations.DeflatedPoint,2,1),LineWidth=1.2)
semilogy(vecnorm([oldLPIterations.Iterates,oldNewtonIterations.Iterates(:,2:end)]-oldNewtonIterations.DeflatedPoint,2,1),LineWidth=1.2)

set(gca, 'YScale', 'log')
hold off
% semilogy(1:length(ek1),ek1,1:length(ek2),ek2)





%% Aux Funcs



function [LPIterations,NewtonIterations,x,F,ek,Eigs] =RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,doubled)
if LPSteps>0
 problem.MaxIter = LPSteps;  problem.StepTolerance = epsilon; problem.doubled = doubled;
LPIterations = RGD_LP(problem);
x1= LPIterations.FinalPoint;
else
    LPIterations=[];
    x1=problem.x0;
end
problem.MaxIter = NewtonSteps; problem.StepTolerance = 1e-8; problem.x0=x1;
[NewtonIterations] = Newton(problem);

x= NewtonIterations.FinalPoint;
ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
[F]=problem.obj_fun(x,problem);
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end


function [LPIterations,NewtonIterations,x,F,ek,Eigs] =Old_LP_Newton(problem,epsilon,LPSteps,NewtonSteps)
if LPSteps>0
 problem.MaxIter = LPSteps; problem.StepTolerance = epsilon;
LPIterations = OGLP(problem);
else
    LPIterations=[];
    x1=problem.x0;
end
x1= LPIterations.FinalPoint;
problem.MaxIter = NewtonSteps; problem.StepTolerance = 1e-8; problem.x0=x1;
[NewtonIterations] = Newton(problem);

x= NewtonIterations.FinalPoint;
ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
[F]=problem.obj_fun(x,problem);
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end





function [x,ek,Eigs,F,LPIterations,NewtonIterations] =LP_Newton(problem,epsilon,LPSteps,Opt,doubled,original)
[~,Binv] = FormB(problem);
% if nargin>2 &&  varargin{1}
if doubled 
    ScalingMatrix = 2*Binv;
    else
    ScalingMatrix = Binv;
end

if isempty(Opt)
    Opt = struct;
end
if isempty(LPSteps)
    LPSteps = inf;
end
if isfield(Opt,'MaxIter')
    NewtonMaxIter = Opt.MaxIter;
else
    NewtonMaxIter = 1e3;
end
if ~isfield(Opt,'Linesearch')
    Opt.Linesearch = 'No';
end

Opt.SupressConvergedFlag = true; 
Opt.NDeflations = 1; Opt.constants = problem;
if original 
    Opt.Method = 'LP'; 
else
    Opt.Method = 'UGradientDescent'; 
end
Opt.ScalingMatrix = ScalingMatrix; Opt.StepTolerance = epsilon; Opt.MaxIter = LPSteps;
[~,LPIterations] = evalc('DMin(problem.obj_fun,problem.x0,Opt)'); 
x1= LPIterations.DeflatedPoint;

Opt.Method = 'UNewton'; Opt.ScalingMatrix = []; Opt.StepTolerance = 1e-10; 
Opt.MaxIter =NewtonMaxIter;
for i = 1:1
[NewtonIterations,res,params] = eval('DMin(problem.obj_fun,x1,Opt)');  
end
x= NewtonIterations.DeflatedPoint;
ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
[F]=problem.obj_fun(x,problem);
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end


function [B,Binv]= FormB(constants)
B = zeros(length(constants.A));
for i = 1:(length(constants.A))
    for j = 1:(length(constants.A))
        B(i,j) = sum(sum(constants.A{j}'.*constants.A{i}));
    end
end
Binv = inv(B);
end


function [f,varargout] = IEP(x,constants)
[QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');
if any(isnan(DFull))
    f = nan;
    for i = 1:nargout-1
        varargout{i} = nan;
    end
    return
end
if length(constants.ev)<length(DFull)
    C = (DFull'-constants.ev).^2;
    pairs = matchpairs(C,1e10);
    D = DFull(pairs(:,2));
    Q = QFull(:,pairs(:,2));
    constants.ev = constants.ev((pairs(:,1)));
else
    D = DFull; Q = QFull;
end

f = sqrt(sum(((D-constants.ev)).^2));
if nargout>1
    varargout{1} = (D-constants.ev);
end
if nargout>2
    varargout{2} =  FormJ_Lambda(Q,constants.A);
end
if nargout>3
    varargout{3} = FormH_Lambda(QFull,constants.A,D,constants.ev);
end
end

function [F,R,J,H] = INSEvaulateDifference(x,constants)
A = constants.A; 
    Ad = x(1)*A{1};
    for i = 2:length(x)
        Ad = Ad + x(i)*A{i};
    end
    [Q,D] = eig(full(Ad),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
    R = ((D(2:end) - D(1:end-1)) - (constants.ev(2:end) - constants.ev(1:end-1)));
    F = sqrt(sum((R).^2));
    if nargout>2
        l = length(A);
        m = size(Q,2);
        LJ = zeros(m,l);
        for k = 1:l
            LJ(:,k) =real(sum((Q.'*A{k}).*Q',2)); 
        end
        J = LJ(2:end,:) - LJ(1:end-1,:);
    end
    if nargout>3
        m = length(constants.ev); l=length(A);
        LH = zeros(l,l,m);
        QAQ = cell(1,l);
        for i = 1:l
            QAQ{i} = Q'*A{i}*Q;
            QAQ{i} =QAQ{i}(1:m,1:m);
        end
        DD=D'-D;
        DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                LH(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
        H = LH(:,:,2:end) - LH(:,:,1:end-1);
    end
end

function [F,R,J,H] = INSEvaulate(x,constants)
A = constants.A; 
    Ad = x(1)*A{1};
    for i = 2:length(x)
        Ad = Ad + x(i)*A{i};
    end
    [Q,D] = eig(full(Ad),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
    F = sqrt(sum((D-constants.ev).^2));
    if nargout>1
        R = (D-constants.ev);
    end
    if nargout>2
        l = length(A);
        m = size(Q,2);
        J = zeros(m,l);
        for k = 1:l
            J(:,k) =real(sum((Q.'*A{k}).*Q',2)); 
        end
    end
    if nargout>3
        m = length(constants.ev); l=length(A);
        H = zeros(l,l,m);
        QAQ = cell(1,l);
        for i = 1:l
            QAQ{i} = Q'*A{i}*Q;
            QAQ{i} =QAQ{i}(1:m,1:m);
        end
        DD=D'-D;
        DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                H(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
    end
end
% dlp=[0.8486, 0.8424, -0.0050, 0.3076, -0.5089, 1.6325, -0.0659,1.72764, -0.00038, 1.1018, -1.5155, -0.8286, 1.1952, -0.7433,0.0336, -0.0737, 0.0356, -1.5870, -0.1220, -0.2275]';


