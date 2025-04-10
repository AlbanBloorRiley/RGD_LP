%% Example 1
clear problem A Binv 
N =5000;  %Dimension of matrices 
l = 40;   %Number of free parameters
problem.A0 =sparse(zeros(N));
%Form basis matrices
A{1} = speye(N);
for k = 2:l
    A{k} = (spdiags(1,[-k+1,k-1],N,N));
end
%Set up problem parameters
problem.A = A;   problem.x0 = -ones(l,1);
problem.ev = (-1.1:0.002:-1.061)'*1e2;
problem.StepTolerance= 1e-3; problem.MaxIter = 1000; 
problem.obj_fun = @IEPMin;

%Run timings for all three methods
tic
RGDLPMinIterations = RGD_LP(problem);
RGDLPMinTime = toc

problem.obj_fun = @IEP;
tic
RGDLPIterations = RGD_LP(problem);
RGDLPTime = toc

tic
LPIterations = LP(problem);
LPTime = toc


%% Example 2 - Mn12 
clear all
%Set up units
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;
%Define solution
B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;
B22 = (6.8*10^-4)*meV;
%Simulate eigenvalues
Sys.S = 10;
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = ham(Sys,[0,0,0]);  [~,E]=eig(H);
EE = diag(E);  Exp.ev=EE-EE(1);

%The Stevens Operators
A0 = sparse(21,21);
A{1} = stev(10,[2,0]);
A{2} = stev(10,[4,0]);
A{3} = stev(10,[4,4]);
A{4} = stev(10,[2,2]);
%Ground State operator
A{5} = eye(21);
%Set up problem parameters
problem.A = A; problem.ev = EE; problem.A0 = A0;
problem.x0 = round([B20;B40;B44;B22;-1e6],1,'significant');
problem.x0 = [-1000,1,1,1,0]';
problem.StepTolerance= 1e-8; problem.MaxIter = 500; repeats = 100;

%Run timings for all three methods
problem.obj_fun = @IEPMin;
tic
for i = 1:repeats
    RGDLPMinIterations = RGD_LP(problem);
end
RGDLPMinTime = toc/repeats

problem.obj_fun = @IEP;
tic
for i = 1:repeats
    RGDLPIterations = RGD_LP(problem);
end
RGDLPTime = toc/repeats

tic
for i = 1:repeats
    LPIterations = LP(problem);
end
LPTime = toc/repeats


%% Example 4 - Cr6 
clear all
%Set up units
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;

%Simulate eigenvalues
N_electrons = 6;
B20 = (-0.041/3).*meV; % convert from meV to MHz
B22 = 0.007.*meV; % convert from meV to MHz
S=3/2;
Sys.S = S;
for i = 2:N_electrons;Sys.S = [Sys.S,S]; end
B2 = [B22 0 B20 0 0]; % B(k=2,q) with q = +2,+1,0,-1,-2
B2s = [1 0 1 0 0];
Sys1.B2 = B2s;
for i = 2:N_electrons; Sys1.B2 = [Sys1.B2;B2s]; end
Jval = 1.46*meV;
J =[]; %ee = [];
for i = 2:N_electrons
    J = [1,zeros(1,i-2),J];
end
Sys1.J = J;
Sys1.S = Sys.S;

Sys.B2 = Sys1.B2.*B2;
Sys.J = Jval.*J;
H=ham(Sys,[0,0,0],'sparse');
[Vecs,EE] = eig(full(H),'vector');
Eigs=EE(1:24)-EE(1);

%Form basis matrices
A{1} = stev(Sys1.S,[2,2,1],'sparse');
A{2} = stev(Sys1.S,[2,0,1],'sparse');
for i = 2:N_electrons       
    A{1} = A{1} + stev(Sys1.S,[2,2,i],'sparse');
    A{2} = A{2} + stev(Sys1.S,[2,0,i],'sparse');
end   
A{3} = ham_ee(Sys1,1:N_electrons,'sparse')*1;

%Set up problem parameters
A0 = sparse(length(A{1}),length(A{1}));
problem.A = A; problem.ev = Eigs; problem.A0 = A0;
problem.Solver = @mldivide;  problem.obj_fun = @IEP;
problem.A{4} = speye(size(A{1}));

%Form B and Binverse matrices 
tic;[Binv]= FormBinv(problem.A);toc
problem.Binv = Binv;
%Form spy plots for Example 3
f = spyplots(1,3,A{1},A{2},A{3});
f.Units = 'centimeters';
f.Position = [-50 10 20 14];
print(f, 'Figures/SpyPlots.eps', '-depsc')
%%
%Set up grid points
Npoints = 5;
X1 = logspace(3,7,Npoints);
X2 = -logspace(3,7,Npoints);
X3 = logspace(3,7,Npoints);
X4 = logspace(3,7,Npoints);
[x1,x2,x3,x4]=ndgrid(X1,X2,X3,X4);
l = 4; 
NIterates = Npoints^l;


%Run timings for all three methods
RGDMinx = zeros(l,NIterates);
tic
for i = 1:NIterates
    [F,R,J] = IEPMin([x1(i);x2(i);x3(i);x4(i)],problem);
    p(:,i) = - Binv*J'*R;
    RGDMinx(:,i) = [x1(i);x2(i);x3(i);x4(i)] +p(:,i);
end
RGDLPMinTime = toc

RGDx = zeros(l,NIterates);
tic
for i = 1:NIterates
    [F,R,J] = IEP([x1(i);x2(i);x3(i);x4(i)],problem);
    p = - Binv*J'*R;
    RGDx(:,i) = [x1(i);x2(i);x3(i);x4(i)] +p;
end
RGDLPTime = toc

tic
x = zeros(l,NIterates);
for i = 1:NIterates
    [QFull,DFull] = eig(full(FormA([x1(i);x2(i);x3(i);x4(i)],problem.A,problem.A0)),'vector');
    if length(problem.ev)<length(DFull)
        C = (DFull'-problem.ev).^2;
        pairs = matchpairs(C,100*max(max(C)));
        Dcomp = DFull;
        Dcomp(pairs(:,2))=[];
        Qmatch = QFull(:,pairs(:,2));
        Qcomp = QFull;
        Qcomp(:,pairs(:,2)) =[];
        
         Z = [Qmatch,Qcomp]*diag([problem.ev;Dcomp])*[Qmatch,Qcomp]';
    else
        Z = QFull*diag(problem.ev)*QFull';
    end
    b = zeros(l,1);
    for j = 1:l
        b(j,1) = trace((Z-problem.A0)'*problem.A{j});
    end
    x(:,i) = Binv*b;
end
LPTime = toc



%% Example 5 - Mn6 
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;

%Define Eigenvalues
problem.ev = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]'.*meV;  
problem.ev =  [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]'.*meV;
Sys1.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

%Form basis matrices
 Sys1.J = [10*meV,0 0 0 0 0 0 0 0 0 0 0 0 0 0];
 problem.A0 = ham_ee(Sys1,1:length(Sys1.S),'sparse'); problem.A0 =sparse(32400,32400);
A{1} = stev(Sys1.S,[2,0,1],'sparse')+ stev(Sys1.S,[2,0,2],'sparse');
Sys1.J = [0,1 0 1 0 0 1 0 1 0 0 0 0 0 0];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');
Sys1.J = [0,0 1 0 1 1 0 1 0 0 0 0 0 0 0];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');
Sys1.J = [0 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse');
Sys1.J = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
A{end+1} = ham_ee(Sys1,1:length(Sys1.S),'sparse'); 

%Set up problem parameters
problem.Solver = @mldivide;  problem.obj_fun = @IEPMin;
% problem.x0 = [-21035 -1.9826e+05 1.9343e+05  48357 10*meV]';
problem.x0 = [1,1,1,1,1]'*1e2;
problem.x0(end+1) = 2e7;
problem.A = A; problem.A{end+1} = speye(32400,32400);
problem.Binv = FormBinv(problem.A);
problem.Verbose = true;
problem.StepTolerance= 1e-1; problem.MaxIter = 600; 
problem.StepDifferenceTolerance = 1e-5;problem.StepRelativeTolerance = 1e-7;

%Run timings for RGDLPMin
tic
    RGDLPMinIterations = RGD_LP(problem);
Mn6Time = toc

%plot convergence of method
f = figure(1);
clf
% semilogy(sum((RGDLPMinIterations.Iterates(:,2:end)-RGDLPMinIterations.Iterates(:,1:end-1)).^2), LineWidth=2)
[Binv,B] = FormBinv(problem.A);

p =(RGDLPMinIterations.Iterates(:,2:end)-RGDLPMinIterations.Iterates(:,1:end-1));
% semilogy(sqrt(sum(p.^2)./(sum(RGDLPMinIterations.Iterates(:,end-1).^2))));ylabel("||p^k||/||x^k||")
% semilogy(sqrt(sum(p.^2)));ylabel("||p^k||")
% semilogy(sqrt(sum(B*p.*p,1)));ylabel("||p^k||_g")
semilogy(sqrt(sum(B*p.*p,1))./sqrt(sum(RGDLPMinIterations.Iterates(:,end-1)'*B*RGDLPMinIterations.Iterates(:,end-1),1)),LineWidth=2);       ylabel("||p^k||_g/||x^k||_g")
yscale('log')
xlabel("Iteration")
grid("on")


f.Units = 'centimeters';
f.Position = [-50 10 20 14];
% linestyleorder('mixedstyles')

% print(f, ["Figures/Mn6Convergence"+num2str(problem.x0')+".png"], '-dpng')
print(f, "Figures/Mn6Convergence100.eps", '-depsc')

%%
%Can check validity of INS simulation using Mint (https://mlbakerlab.co.uk/mint/)
% Iterates = problem.x0;
Iterates = RGDLPMinIterations.FinalPoint;
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2
% J_S4_S4   = RGDLPMinIterations.FinalPoint(5)./(-2); %
% J_S4_S4   = -5*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
J_S4_S5_1 = Iterates(2)./(-2);  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = Iterates(3)./(-2); % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = Iterates(4)./(-2); 
J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = -0.00.*meV;         J_23 = J_14; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = -0.00.*meV;         J_24 = J_13; 
Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2JS.S formalism
DIII = Iterates(1)./3; % MnIII anisotropy. Free Value 
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


%% Aux Funcs


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


function plotmint(cross_sect, MintExp, colours)
subplot(2,1,1)
plot(MintExp.Energy,cross_sect*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('New Sim')
ax = gca;ax.ColorOrder = colours;
%
subplot(2,1,2)
load("Sys0_Sim.mat")
plot(MintExp.Energy,cross_sect_Sys0*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('Old Sim')
ax = gca;ax.ColorOrder = colours;
f.Units = 'centimeters';
f.Position = [10 10 20 25];
end

