%% Example 1
clear all
N =1000;
problem.A0 =sparse(zeros(N));
for k = 1:floor(N/3)
    D= zeros(N,1); D(k) = 1;
    problem.A{k} =  sparse(toeplitz(D));
end
problem.x0 = [1:floor(N/3)]';
problem.ev =  [-floor(N/4):-1,1:ceil(N/4)]';
% problem.ev =[-N/2:-1,1:N/2]';
   problem.obj_fun = @IEP;
repeats = 1;
problem.Solver = @mldivide;

% tic
% for i = 1:repeats
%     [~,OnlyNewtonIterations] =RGD_LP_Newton(problem,epsilon,0,500,doubled);
% end
% NewtonTime = toc/repeats;
% ek1 = vecnorm(OnlyNewtonIterations.Iterates(:,2:end)-OnlyNewtonIterations.Iterates(:,1:end-1),2,1);
% %
epsilon = 0.00001; LPSteps = 1000;NSteps =0;
tic
for i = 1:repeats
    [RGDLPIterations01,RGDNewtonIterations01] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps);
end
RGDN01Time = toc/repeats;
%
RGDNIterates01 = [RGDLPIterations01.Iterates,RGDNewtonIterations01.Iterates(:,2:end)];
ek2 = vecnorm(RGDNIterates01(:,2:end)-RGDNIterates01(:,1:end-1),2,1);
% %
% epsilon = 0.001;
% tic
% for i = 1:repeats
%     [RGDLPIterations001,RGDNewtonIterations001] =RGD_LP_Newton(problem,epsilon,LPSteps,NSteps);
% end
% RGDN001Time  = toc/repeats;
% RGDNIterates001 = [RGDLPIterations001.Iterates,RGDNewtonIterations001.Iterates(:,2:end)];
% ek3 = vecnorm(RGDNIterates001(:,2:end)-RGDNIterates001(:,1:end-1),2,1);
% 
% epsilon = 0.01;
%
tic
for i = 1:repeats
    [LPIterations01,NewtonIterations01] =LP_Newton(problem,epsilon,LPSteps,NSteps);
end
Old01Time  = toc/repeats;
oldIterates01 = [LPIterations01.Iterates,NewtonIterations01.Iterates(:,2:end)];
ek4 = vecnorm(oldIterates01(:,2:end)-oldIterates01(:,1:end-1),2,1);

% epsilon = 0.001;
% tic
% for i = 1:repeats
%     [LPIterations001,NewtonIterations001,x,F,ek,Eigs] =Old_LP_Newton(problem,epsilon,LPSteps,NSteps);
% end
% Old001Time  = toc/repeats;
% oldIterates001 = [LPIterations001.Iterates,NewtonIterations001.Iterates(:,2:end)];
% ek5 = vecnorm(oldIterates001(:,2:end)-oldIterates001(:,1:end-1),2,1);
%%

difference01 = sum(abs(LPIterations01.Iterates(:,1:min(length(LPIterations01.Iterates),length(RGDLPIterations01.Iterates)))-RGDLPIterations01.Iterates(:,1:min(length(LPIterations01.Iterates),length(RGDLPIterations01.Iterates)))),1);
% difference001 = sum(abs(LPIterations001.Iterates(:,1:min(length(LPIterations001.Iterates),length(RGDLPIterations001.Iterates)))-RGDLPIterations001.Iterates(:,1:min(length(LPIterations001.Iterates),length(RGDLPIterations001.Iterates)))),1);

plot(1:length(difference01),difference01)


%%



["Algorithm", "No. Iterations","CPU Time (seconds)";...
    "Newton", OnlyNewtonIterations.NIter, NewtonTime ;...
    "Old LP-Newton, \epsilon = 0.01", string(LPIterations01.NIter+"+"+NewtonIterations01.NIter),Old01Time;...
    "Old LP-Newton, \epsilon = 0.001", string(LPIterations001.NIter+"+"+NewtonIterations001.NIter),Old001Time;...
    "New LP-Newton, \epsilon = 0.01", string(RGDLPIterations01.NIter+"+"+RGDNewtonIterations01.NIter),RGDN01Time;...
    "New LP-Newton, \epsilon = 0.001", string(RGDLPIterations001.NIter+"+"+RGDNewtonIterations001.NIter),RGDN001Time
    ]
clf

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

%% Example 2
N = 10;
problem.A0 = toeplitz([2,-1,zeros(1,N-2)])';
for k = 1:N
    problem.A{k} =  sparse(diag([zeros(1,k-1),1,zeros(1,N-k)]));
end
problem.x0 = [1:N]';
problem.ev =  [-floor(N/4):-1,1:ceil(N/4)]';
% problem.ev =[-N/2:-1,1:N/2]';

   problem.obj_fun = @IEP;
repeats = 1;
problem.Solver = @mldivide;
LPSteps = 1000; NSteps = 0;epsilon = 0.001;

[RGDLPIterations01,RGDNewtonIterations01] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps);


[LPIterations01,NewtonIterations01] = LP_Newton(problem,epsilon,LPSteps,NSteps);
%% Example 3 - Mn12 
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;
B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;
B22 = (6.8*10^-4)*meV;

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
problem.A = A; problem.ev = EE; problem.A0 = A0;
problem.x0 = round([B20;B40;B44;B22;-1e6],2,'significant');

% for i = 1:length(A)
%     problem.A{i} = problem.A{i}*problem.x0(i);
% problem.x0(i) = 1;
% end
%%

problem.Solver = @mldivide;  problem.obj_fun = @INSEvaulate;
repeats = 50;

epsilon = 0.001; LPSteps = Inf; doubled = false; NewtonSteps = 0; doubled = false;
tic
for i = 1:repeats
    [newIterations,~] = RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,doubled);
end; RGDLP01Time = toc/repeats;
ekLPnew = vecnorm(newIterations.Iterates(:,2:end)-newIterations.Iterates(:,1:end-1),2,1);

tic
for i = 1:repeats
    oldIterations = LP_Newton(problem,epsilon,LPSteps,NewtonSteps);
end;  LP01Time = toc/repeats;
ekLPold = vecnorm(oldIterations.Iterates(:,2:end)-oldIterations.Iterates(:,1:end-1),2,1);


epsilon = 0.00001; LPSteps = Inf; doubled = false; NewtonSteps = 0; doubled = false;
tic
for i = 1:repeats
    [newIterations001,~] = RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,doubled);
end; RGDLP001Time = toc/repeats;
ekLPnew = vecnorm(newIterations001.Iterates(:,2:end)-newIterations001.Iterates(:,1:end-1),2,1);

tic
for i = 1:repeats
    oldIterations001 = LP_Newton(problem,epsilon,LPSteps,NewtonSteps);
end;  LP001Time = toc/repeats;
ekLPold = vecnorm(oldIterations001.Iterates(:,2:end)-oldIterations001.Iterates(:,1:end-1),2,1);



["Algorithm", "No. Iterations","CPU Time (seconds)";...
    " RGD-LP, \epsilon = 0.01", newIterations.NIter ,RGDLP01Time;...
    "RGD-LP, \epsilon = 0.001", newIterations001.NIter,RGDLP001Time;...
    "LP, \epsilon = 0.01", oldIterations.NIter,LP01Time;...
    "LP, \epsilon = 0.001", oldIterations001.NIter,LP001Time
    ]




% epsilon = 0.01; LPSteps = Inf; doubled = false; NewtonSteps = Inf;
% tic
% for i = 1:repeats
%     [newLPIterations,newNewtonIterations] =RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,doubled);
% end
% LPNnewTime =toc/repeats;
% newIterates = [newLPIterations.Iterates,newNewtonIterations.Iterates(:,2:end)];
% ekLPNnew = vecnorm(newIterates(:,2:end)-newIterates(:,1:end-1),2,1);

% 
% f=figure(1);
% clf
% 
% semilogy(vecnorm(oldIterations.Iterates-oldIterations.FinalPoint,2,1),LineWidth=1.5)
% hold on
% semilogy(vecnorm(newIterations.Iterates-newIterations.FinalPoint,2,1),LineWidth=1.5)
% % semilogy(vecnorm([newLPIterations.Iterates,newNewtonIterations.Iterates(:,2:end)]-newNewtonIterations.FinalPoint,2,1),LineWidth=1.5)
% set(gca, 'YScale', 'log')
% hold off
% legend({"RGD LP","LP","LP-Newton"},Location='ne')
% % ylim([1e-10,1e40])
% xlabel("The number of iterations")
% ylabel("log_{10}(error)")
% %%
% f.Units = 'centimeters';
% f.Position = [-50 10 20 14];
% setMarkerNumber(f.Children(2),20)
% linestyleorder('default')
% print(f, 'ExampleMn12.eps', '-depsc')


%% Example 4 - Cr6 
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;

N_electrons = 6;
B20 = (-0.041/3).*meV; % convert from meV to MHz
B22 = 0.007.*meV; % convert from meV to MHz
% BValues = [-3304.4;1692.5;353000];
% We have 6 S=3/2 spins that are coupled together
S=3/2;
Sys.S = S;
for i = 2:N_electrons;Sys.S = [Sys.S,S]; end
% n = prod(2*Sys.S+1);
% each of the spins has a non zero B22 and B20 
B2 = [B22 0 B20 0 0]; % B(k=2,q) with q = +2,+1,0,-1,-2
B2s = [1 0 1 0 0];
Sys1.B2 = B2s;
for i = 2:N_electrons; Sys1.B2 = [Sys1.B2;B2s]; end
Jval = 1.46*meV;
J =[]; %ee = [];
for i = 2:N_electrons
    J = [1,zeros(1,i-2),J];
    % ee = [eye(3);zeros(3*(i-2),3);ee];
end
% Sys.ee = Jval.*ee;
Sys1.J = J;
Sys1.S = Sys.S;

Sys.B2 = Sys1.B2.*B2;
Sys.J = Jval.*J;

A{1} = stev(Sys1.S,[2,2,1],'sparse');
A{2} = stev(Sys1.S,[2,0,1],'sparse');
for i = 2:N_electrons       
    A{1} = A{1} + stev(Sys1.S,[2,2,i],'sparse');
    A{2} = A{2} + stev(Sys1.S,[2,0,i],'sparse');
end   
A{3} = ham_ee(Sys1,1:N_electrons,'sparse');
A{4} = speye(size(A{1}));
A0 = sparse(length(A{1}),length(A{1}));
H=ham(Sys,[0,0,0],'sparse');
[Vecs,EE] = eig(full(H),'vector');
Eigs=EE(1:24)-EE(1);

problem.A = A; problem.ev = Eigs; problem.A0 = A0;
problem.Solver = @mldivide;  problem.obj_fun = @INSEvaulate;

f = figure(1);
subplot(1,3,1)
spy(problem.A{1})
subplot(1,3,2)
spy(problem.A{2})
subplot(1,3,3)
spy(problem.A{3})
f.Units = 'centimeters';
f.Position = [-50 10 20 14];
print(f, 'SpyPlots.eps', '-depsc')

problem.x0 = [1e3,-1e3,1e5,1]';
problem.x0 = [100,-100,1000,100000]';
problem.x0 = [1,1,1,1]';
% for i = 1:length(A)
%     problem.A{i} = problem.A{i}*problem.x0(i);
% problem.x0(i) = 1;
% end

%%
repeats = 1;
epsilon = 1e0; LPSteps = 100; doubled = false; NewtonSteps = 0;
tic
for i = 1:repeats
    % [LPIterations,NewtonIterations] = Old_LP_Newton(problem,epsilon,LPSteps,NewtonSteps);
    [LPIterations,NewtonIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps, doubled);
end
LPnewTime = toc/repeats;
NewtonIterations
%%

epsilon = 1e-10; LPSteps = 10; doubled = false; NewtonSteps = 0;
tic
for i = 1:repeats
    % [LPIterations,NewtonIterations] = O ld_LP_Newton(problem,epsilon,LPSteps,NewtonSteps);
    [LPIterations,NewtonIterations] = LP_Newton(problem,epsilon,LPSteps,NewtonSteps);
end
LPnewTime = toc/repeats;



% ekLPnew = vecnorm(newIterations.Iterates(:,2:end)-newIterations.Iterates(:,1:end-1),2,1);




%% Aux Funcs



function [LPIterations,NewtonIterations,x,F,ek,Eigs] =RGD_LP_Newton(problem,epsilon,LPSteps,NewtonSteps,varargin)
if ~isempty(varargin)
problem.doubled = varargin{1};
else
    problem.doubled = false;
end
if LPSteps>0
    problem.MaxIter = LPSteps;  problem.StepTolerance = epsilon; 
    LPIterations = RGD_LP(problem);
    x1= LPIterations.FinalPoint;
else
    LPIterations.Iterates=problem.x0;
    LPIterations.NIter = 0;
    x1=problem.x0;
end
if NewtonSteps>0
    problem.MaxIter = NewtonSteps; problem.StepTolerance = 1e-8; problem.x0=x1;
    [NewtonIterations] = Newton(problem);
    x= NewtonIterations.FinalPoint;
    ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
else
    NewtonIterations.Iterates  =x1;
    NewtonIterations.NIter = 0;
    x = x1;
    ek = norm(LPIterations.Iterates(:,end)-LPIterations.Iterates(:,end-1));
end
[F]=problem.obj_fun(x,problem);
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end


function [LPIterations,NewtonIterations,x,F,ek,Eigs] =LP_Newton(problem,epsilon,LPSteps,NewtonSteps)
if LPSteps>0
    problem.MaxIter = LPSteps; problem.StepTolerance = epsilon;
    LPIterations = OGLP(problem);
    x1= LPIterations.FinalPoint;
else
    LPIterations.Iterates=problem.x0;
    LPIterations.NIter = 0;
    x1=problem.x0;
end
if NewtonSteps>0
    problem.MaxIter = NewtonSteps; problem.StepTolerance = 1e-8; problem.x0=x1;
    [NewtonIterations] = Newton(problem);
    x= NewtonIterations.FinalPoint;
    ek = norm(NewtonIterations.Iterates(:,end)-NewtonIterations.Iterates(:,end-1));
else
     NewtonIterations.Iterates  =x1;
     NewtonIterations.NIter = 0;
    x = x1;
    ek = norm(LPIterations.Iterates(:,end)-LPIterations.Iterates(:,end-1));
end
[F]=problem.obj_fun(x,problem);
Eigs = eig(full(FormA(x,problem.A,problem.A0)),'vector');
end





function [f,varargout] = IEP(x,constants)
% 
if length(constants.A{1})<500 && length(constants.ev)<0.5*length(constants.A{1}) &&nargout<4
    [Q,D] = eigs(FormA(x,constants.A,constants.A0),length(constants.ev),constants.ev);
else
    [QFull,DFull] = eig(full(FormA(x,constants.A,constants.A0)),'vector');


    if length(constants.ev)<length(DFull)
        C = (DFull'-constants.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        D = DFull(pairs(:,2));
        Q = QFull(:,pairs(:,2));
        % constants.ev = constants.ev((pairs(:,1)));
    else
        D = DFull; Q = QFull;
    end
end
if any(isnan(D))
    f = nan;
    for i = 1:nargout-1
        varargout{i} = nan;
    end
    return
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




function [F,R,J,H] = INSEvaulate(x,constants)
A = constants.A;
Ad = x(1)*A{1};
for i = 2:length(x)
    Ad = Ad + x(i)*A{i};
end
if length(A{1})>500 && length(constants.ev)<0.5*length(A{1})&&nargout<4
    [Q,D] = eigs(Ad, length(constants.ev), 'smallestreal');
    D = diag(D);
else
    [Q,D] = eig(full(Ad),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
end
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


function setMarkerNumber(f,n)
    for i = 1:length(f.Children)
        f.Children(i).MarkerIndices = 1:(floor(length(f.Children(i).MarkerIndices)/n)):f.Children(i).MarkerIndices(end);
    end
end
