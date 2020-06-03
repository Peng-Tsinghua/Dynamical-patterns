% This program simulates the dynamical trajectories of the Illinois power
% grid under different voltage damings in the nominal operating condition
% It automatically outputs 100 data files (ODE_down_x.mat, ODE_up_x.ma)
% that store the system trajectories. Based on these trajectories, results
% in Fig.6 can be obtained.
%% Load basic data and set parameters
clear,clc
load('0Illinois200');
% get nominal operation point by AC-pf solover
define_constants;
mpc0=loadcase('case_illinois200');
LDbus=find(mpc0.bus(:,PD)>0); % find load bus index
nLD=length(LDbus); % total number of loads
indexSG=find(mpc0.gen(:,GEN_STATUS)==1); % index of turn-on SG
indexGen=find(mpc0.gen(:,PG)>0); % index of SG in matrix mpc.gen
SGbus=mpc0.gen(indexSG,1); 
nSG=length(SGbus); % total number of generators
n0=size(mpc0.bus,1); % total number of buses
result0=runpf(mpc0);
% estimate xd', xd, and H for generators
S0=sqrt(result0.gen(indexSG,QG).^2+result0.gen(indexSG,PG).^2);
Q0=mpc0.gen(indexSG,QMAX); % QMAX is provided in the case
V0=mpc0.gen(indexSG,VM);
QV0=abs(Q0)./V0;
xdd=17.91./QV0;
xd=7.398*xdd-0.1333;
H=3.312*S0+3.308;
%% Set more parameters
% Initial network parameters
Bbus0=Bbus-diag(diag(Bbus));
Bbus0=Bbus0-diag(sum(Bbus0));% Bbus0: whiteout bi
b=diag(Bbus0)-diag(Bbus); % shunt,b<0 for inductance
n=length(b); % number of nodes
n2=nSG; % number of second order dynamics (SG)
n1=n-n2; % number of first order dynamics (loads)
U=Vd(:,1);V0=abs(U);the0=angle(U);
I=1i*Bbus0*U;S=U.*conj(I);
Pg=real(S);Qg=imag(S);
% set dynamical parameter
Td0=5+5.*rand(n2,1); % random from [5,10]
T2=Td0./(xd-xdd);K2=1./(xd-xdd)-b(n1+1:end)+1;
T1=ones(n1,1)*20; % set T1=20
D=0.1*ones(n,1); % set frequency damping=0.1 for loads
D(end-n2+1:end)=1; % set frequency damping=1 for generations
K1=-b(1:n1)+1; % set K1=1
M=H.*2; % SG inertia M, M=2H
T=[T1;T2];
K0=[K1;K2];
%% theoretical phase transition points
Bzz=func_Bzz(the0,V0,Bbus0);
Schur=Bzz(n+1:end,n+1:end)-Bzz(1:n,n+1:end)'*pinv(Bzz(1:n,1:n))*Bzz(1:n,n+1:end);
min(eig(diag(K0)+Schur))
min(eig(diag(K0)-Bbus0))
%% set initial values by flawed direction
x0=[the0;V0;zeros(n2,1)];
[Ve,De]=eig(Bzz);
find(real(diag(De))<0)
temp=diag(De); temp([38,220])
directBzz=Ve(:,220);
direct0=[directBzz;zeros(n2,1)];
direct=0.01*direct0;
x00=x0+direct;x000=x0-direct;
%% solve system trajectories under different K
m=50;
kplus=linspace(-1.15,-0.9,m);
R=zeros(2,m); v=R;
options = odeset('Events',@nbusEventsFcn,'NonNegative',[1+n:2*n]...
    ,'Jacobian',@Jacnbus); % nonnegative constrain of voltages
fields = {'extdata','xe','ye','ie','stats','idata'};

for i=1:m
    clear sol
    K=K0+kplus(i);
     u=K.*V0+Qg./V0;
    % downwards instability, use solver ode23
    [t1,x1]=ode15s(@ode_nbus,[0 10000],x00,options,Bbus0,M,D,T,K,Pg,u);
    if t1(end)<10000
        [t2,x2]=ode23(@ode_nbus,[t1(end-1) 10000],x1(end-1,:)',options,Bbus0,M,D,T,K,Pg,u);
        t=[t1(1:end-1);t2];x=[x1(1:end-1,:);x2]';
    else
        t=t1;x=x1';
    end
    % compress data if neccessary
%     index=1:floor(length(t)/4);
%     sol.x=t(4*index);sol.y=x(:,4*index);
    sol.x=t;sol.y=x;
    save(['ODE_down_i=',num2str(i)],'sol')
    % upwards instability, stiff, use solver ode15s
    [t,x]=ode15s(@ode_nbus,[0 10000],x000,options,Bbus0,M,D,T,K,Pg,u);
    sol.x=t;sol.y=x';
    save(['ODE_up_i=',num2str(i)],'sol')
end