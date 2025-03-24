%% Example 1
clear problem A Binv 
N =5000;
l = 40;
problem.A0 =sparse(zeros(N));

A{1} = speye(N);
Aindexed.i = [1:N]';  Aindexed.j = [1:N]';  Aindexed.v = ones(N,1);
Aindexed.m = N;Aindexed.n = N;Aindexed.NumV = N;
for k = 2:l
    A{k} = (spdiags(1,[-k+1,k-1],N,N));
    nz = nnz(A{k});
    [Aindexed.i(end+1:end+nz),Aindexed.j(end+1:end+nz),Aindexed.v(end+1:end+nz)] =find(A{k});
    % H{k,k,k} =find(A{k});
    Aindexed.NumV(k) = nz;
end
problem.A = A;

step = 10;
problem.x0 = 10*[1:step:step*l]';
% problem.x0 =([-60,-3,0,3,60]')*1e1;
problem.x0 = ones(l,1);
% problem.x0 = [-N/2:-1,1:N/2]';
% problem.x0 = [1.1650, 0.6268, 0.0751, 0.3516, -0.6965, 1.6961, 0.0591, 1.7971, ...
%     0.2641, 0.8717, -1.4462, -0.7012, 1.2460, -0.6390, 0.5773, -0.3600,...
%     -0.1356, -1.3493, -1.2704, 0.9845]';      %For N = 20 (original example 2)
problem.ev =  [-floor(N/4):-1,0,1:ceil(N/4)]';
problem.ev = [-5,-4,-3,-2,-1]'*100;
problem.ev = -[10,9,8,7,6,5,4,3,2,1]'*1e1;
problem.ev = (-1.1:0.001:-1.081)'*5e2;

% eigs(FormA(problem.x0,problem.A,problem. A0),21,'smallestreal')
% eigs(FormA(problem.x0,problem.A,problem.A0),30,'smallestreal')


   problem.obj_fun = @(a,b)IEPsmallest(a,b,false);
repeats = 1;
problem.Solver = @mldivide;
tic;[Binv]= FormBinv(problem.A);toc
problem.Binv = Binv;
%%
problem.StepTolerance= 1e-6; problem.MaxIter = 500;
%
problem.obj_fun = @(a,b)IEPsmallest(a,b,false);
tic
for i = 1:repeats
    RGDLPMinIterations = RGD_LP(problem)
    % [RGDLPMinIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps,false)
end
RGDLPMinTime = toc/repeats


%%

problem.obj_fun = @IEP;
tic
for i = 1:repeats
    RGDLPMinIterations = RGD_LP(problem)
    % [RGDLPIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps,false)
end
RGDLPTime = toc/repeats
%%
for i = 1:repeats
    [LPIterations,] = LP_Newton(problem,epsilon,LPSteps,NSteps)
end
LPN01Time = toc/repeats


%%
% ["Algorithm", "No. Iterations","CPU Time (seconds)";...
%     "Newton", OnlyNewtonIterations.NIter, NewtonTime ;...
%     "Old LP-Newton, \epsilon = 0.01", string(LPIterations01.NIter+"+"+NewtonIterations01.NIter),Old01Time;...
%     "Old LP-Newton, \epsilon = 0.001", string(LPIterations001.NIter+"+"+NewtonIterations001.NIter),Old001Time;...
%     "New LP-Newton, \epsilon = 0.01", string(RGDLPIterations01.NIter+"+"+RGDNewtonIterations01.NIter),RGDN01Time;...
%     "New LP-Newton, \epsilon = 0.001", string(RGDLPIterations001.NIter+"+"+RGDNewtonIterations001.NIter),RGDN001Time
%     ]
colours = colororder;
f = figure(1);
clf
semilogy(sum((RGDLPIterations.Iterates(:,2:end)-RGDLPIterations.Iterates(:,1:end-1)).^2), LineWidth=1.5)
hold on
semilogy([nan(1,RGDLPIterations.NIter-1),sum((RGDNewtonIterations.Iterates(:,2:end)-RGDNewtonIterations.Iterates(:,1:end-1)).^2)], 'color',colours(1,:),LineWidth=1.5)

hold off


%% Example 2 - Mn12 
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
problem.x0 = round([B20;B40;B44;B22;-1e6],1,'significant');
problem.x0 = [-1000,1,1,1,1]';


problem.Solver = @mldivide;  problem.obj_fun = @IEP;
repeats = 1;

%
epsilon = 1e-10; LPSteps = 500;NSteps =0;

tic
for i = 1:repeats
    [RGDLPIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps,false)
end
RGDLPTime = toc/repeats
%
for i = 1:repeats
    [LPIterations,] = LP_Newton(problem,epsilon,LPSteps,NSteps)
end
LPN01Time = toc/repeats

%

problem.obj_fun = @(a,b)IEPsmallest(a,b,false);
tic
for i = 1:repeats
    [RGDLPMinIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps,false)
end
RGDLPMinTime = toc/repeats


%%

["Algorithm", "No. Iterations","CPU Time (seconds)";...
    " RGD-LP, \epsilon = 0.01", newIterations.NIter ,RGDLP01Time;...
    "RGD-LP, \epsilon = 0.001", newIterations001.NIter,RGDLP001Time;...
    "LP, \epsilon = 0.01", oldIterations.NIter,LP01Time;...
    "LP, \epsilon = 0.001", oldIterations001.NIter,LP001Time
    ]

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
A{3} = ham_ee(Sys1,1:N_electrons,'sparse')*1;


A0 = sparse(length(A{1}),length(A{1}));
% A0 = speye(size(A{1}))*5.21e6;
H=ham(Sys,[0,0,0],'sparse');
[Vecs,EE] = eig(full(H),'vector');
Eigs=EE(1:24)-EE(1);

problem.A = A; problem.ev = Eigs; problem.A0 = A0;
problem.Solver = @mldivide;  problem.obj_fun = @IEP;
problem.A{4} = speye(size(A{1}));
tic;[Binv]= FormBinv(problem.A);toc
problem.Binv = Binv;
% f = spyplots(1,3,A{1},A{2},A{3});
% f.Units = 'centimeters';
% f.Position = [-50 10 20 14];
% print(f, 'Figures/SpyPlots.eps', '-depsc')
%%
problem.x0 = [1692.5,-3304,3.53e5,5.2117e6]';
problem.x0 = [1000,-3000,353000,5.21e6]';

problem.x0 = [1,-1,1,1]'*1e4;
% x0 = [1e3,-1e3,3e5,5e6]';

repeats = 1;
Npoints = 2;

X1 = linspace(0,5e3,Npoints);
X2 = linspace(0,-5e3,Npoints);
X3 = linspace(0,1e6,Npoints);
X4 = linspace(0,1e7,Npoints);
[x1,x2,x3,x4]=ndgrid(X1,X2,X3,X4);

%%
tic
for i = 1:numel(x1)
    [F,R,J] = IEPsmallest([x1(i);x2(i);x3(i);x4(i)],problem);
    p = - Binv*J'*R;
    RGDMinx(:,:,i) = [x1(i);x2(i);x3(i);x4(i)] +p;
end
toc
%%
tic
for i = 1:numel(x1)
    [F,R,J] = IEP([x1(i);x2(i);x3(i);x4(i)],problem);
    p = - Binv*J'*R;
    RGDx(:,:,i) = [x1(i);x2(i);x3(i);x4(i)] +p;
end
toc
%%
tic
for i = 1:numel(x1)
    [QFull,DFull] = eig(full(FormA([x1(i);x2(i);x3(i);x4(i)],problem.A,problem.A0)),'vector');
    if length(problem.ev)<length(DFull)
        C = (DFull'-problem.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        % D = DFull(pairs(:,2));
        Dcomp = DFull;
        Dmatch = DFull(pairs(:,2));
        Dcomp(pairs(:,2))=[];
        Qmatch = QFull(:,pairs(:,2));
        Qcomp = QFull;
        Qcomp(:,pairs(:,2)) =[];
        % constants.ev = constants.ev((pairs(:,1)));
        Q = [Qmatch,Qcomp];
        % if any(pairs(:,2)~=[1:length(pairs)]')
        %     pairs
        % end
    else
         Q = QFull;  %D = DFull;
         Dcomp = [];
    end

    % Z = Q*diag(constants.ev)*Q';
    Z = Q*diag([problem.ev;Dcomp])*Q';
    for j = 1:length(problem.A)
        b(j,1) = trace((Z-problem.A0)'*problem.A{j});
    end
    x(:,:,i) = Binv*b;
end
toc



%% Example 5 - Mn6 
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;

problem.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;  
problem.ev =  [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;
Sys1.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

clear A
 Sys1.J = [10*meV,0 0 0 0 0 0 0 0 0 0 0 0 0 0];
 problem.A0 = ham_ee(Sys1,1:length(Sys1.S),'sparse');% problem.A0 =sparse(32400,32400);
A{1} = stev(Sys1.S,[2,0,1],'sparse')+ stev(Sys1.S,[2,0,2],'sparse');
% A{end+1} = stev(Sys1.S,[2,2,1],'sparse')+ stev(Sys1.S,[2,2,1],'sparse');
Sys1.J = [0,1 0 1 0 0 1 0 1 0 0 0 0 0 0];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');
Sys1.J = [0,0 1 0 1 1 0 1 0 0 0 0 0 0 0];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');
% A{end} = A{end}-ham_ee(Sys1,1:length(Sys1.S),'sparse');


Sys1.J = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');

% Sys1.J = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse'); 


problem.Solver = @mldivide;  problem.obj_fun = @IEPsmallest;
% f = spyplots(2,2,A{1},A{2},A{3},A{4});

%%
problem.x0 = [-21035 -1.9826e+05 1.9343e+05  48357 ]';
%  problem.x0 = [-21035 -1.9826e+05  48357 ]';
% problem.x0 = [-2e4 -2e5 2e5   5e4 ]';
% problem.x0 = [-20000 -200000   200000 50000 2000000]';
problem.x0(end+1) =  -eigs(FormA(problem.x0,A,problem.A0),1,"smallestreal");
problem.A = A; problem.A{end+1} = speye(32400,32400);
problem.Binv = FormBinv(problem.A);
%%
epsilon = 1e-2; LPSteps = 10000;NSteps =0;

tic
[RGDLPIterations] = RGD_LP_Newton(problem,epsilon,LPSteps,NSteps,true)
RGDLPTime = toc

%%


colours = colororder;
f = figure(1);
clf
semilogy(sum((LPIterations.Iterates(:,2:end)-LPIterations.Iterates(:,1:end-1)).^2), LineWidth=1.5)
hold on
semilogy([nan(1,LPIterations.NIter-1),sum((NewtonIterations.Iterates(:,2:end)-NewtonIterations.Iterates(:,1:end-1)).^2)], 'color',colours(1,:),LineWidth=1.5)
% 
%%
[Q,D] = eigs(FormA(LPIterations.FinalPoint, problem.A,problem.A0),length(problem.ev),'smallestreal');
 [Q0,D0] = eigs(FormA(problem.x0, problem.A,problem.A0),length(problem.ev),'smallestreal');
D0 = diag(D0); D0 =  (D0-D0(1))./meV;
D = diag(D); D = (D-D(1))./meV;
MintOpt.Eigs = D;
MintOpt.Vecs = Q;
%%
LPIterations.FinalPoint = problem.x0;
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2
% J_S4_S4   = LPIterations.FinalPoint(5)./(-2); %
J_S4_S4   = -5*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
J_S4_S5_1 = LPIterations.FinalPoint(2)./(-2);  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = LPIterations.FinalPoint(3)./(-2); % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = LPIterations.FinalPoint(4)./(-2); 
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = -0.00.*meV;         J_23 = J_14; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = -0.00.*meV;         J_24 = J_13; 
Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2JS.S formalism
DIII = LPIterations.FinalPoint(1)./3; % MnIII anisotropy. Free Value 
% DIII = -0
% .02*meV;
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism
DII = -0.00*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism
Sys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];
MintSys = Sys;

% Sim INS powder spectrum
MintExp.SpectrumType = 'SE'; %INS intensity vs. energy
Ei = 4; %Incident neutron energy in meV
MintExp.lwfwhm = 0.02*Ei/2.355; %calculated line width full-width-at-half-max, 2 % of incident energy
% MintExp.lwfwhm = 0.05;
MintExp.Energy = linspace(-Ei*0.8,Ei*0.8,1000); %calculating the spectrum in the interval -0.8*Ei to 0.8*Ei
MintExp.Q = 0.1:0.01:2.5; %Q-range which the simulation integrates over.
MintExp.Temperature = [1.5 5 10 30];
MintOpt.NumEigs =length(problem.ev)+1; %100 eigenvalues gives a good INS sim

b =[0 0.4470 0.7410];
r=[0.8500 0.3250 0.0980];
y=[0.9290 0.6940 0.1250];
g=[0.4660 0.6740 0.1880];
colours = [b;y;g;r];

FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};
% Metal ion coordinates for INS. Only relative positions are important.
Coords = [24.60550  13.06674   7.07523; % A
    22.27025  11.49336   6.95626; % B
    21.99544  13.57222   9.30268; % 1
    24.62949  10.95050   9.41655; % 2
    24.24486  10.55088   4.68547; % 3
    22.84040  14.03029   4.66904; % 4
    ];
MintSys.FormFactor = FormFactor;
MintSys.Coords = Coords;

    [cross_sect,Eigs,Vecs,I_nm] = mint(MintSys,MintExp,MintOpt);
    plotmint(cross_sect, MintExp, colours)

%%

%% Aux Funcs




% function setMarkerNumber(f,n)
%     for i = 1:length(f.Children)
%         f.Children(i).MarkerIndices = 1:(floor(length(f.Children(i).MarkerIndices)/n)):f.Children(i).MarkerIndices(end);
%     end
% end
% 
function f = spyplots(x,y,varargin)
l = length(varargin);
if l~= x*y
    error("must be the same dimension")
end
f = figure(1);
for i = 1:l
    subplot(x,y,i)
    spy(varargin{i})
end
end

function A = scaleAx(x,A)
for i = 1:length(x)
    A{i} = A{i}*x(i);
end
end

function plotmint(cross_sect, MintExp, colours)
subplot(2,1,1)
plot(MintExp.Energy,cross_sect*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('New Sim')
ax = gca;;ax.ColorOrder = colours;
%
 subplot(2,1,2)
load("Sys0_Sim.mat")
plot(MintExp.Energy,cross_sect_Sys0*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('Old Sim')
ax = gca;;ax.ColorOrder = colours;
f.Units = 'centimeters';
f.Position = [10 10 20 25];
end

function plotek(Iterations,varargin)
colours = colororder;
f = figure(1);
clf

% for i = 1:length(varargin)
semilogy(sum((Iterations.Iterates(:,2:end)-Iterations.Iterates(:,1:end-1)).^2), LineWidth=1.5)
hold on
if nargin>1
semilogy([nan(1,varargin{1}.NIter-1),sum((RGDNewtonIterations.Iterates(:,2:end)-varargin{1}.Iterates(:,1:end-1)).^2)], 'color',colours(1,:),LineWidth=1.5)
end
hold off
end