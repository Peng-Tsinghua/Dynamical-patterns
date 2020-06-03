function out=ode_nbus(t,x,B,M,D,T,K,Pg,u)
% x=[the,V,omega] in the order of bus index
% B=[suspectance] without shunts
% M=[inertia]
% D=[frequency damping] 
% T=[voltage time constant]
% K=[voltage damping] including shunts
% Pg,u=[input]
n2=length(M);% number of second order dynamics (SG)
n1=length(D)-n2; % number of first order (load)
n=n1+n2;
the=mod(x(1:n),2*pi);
% the=x(1:n);
V=x(n+1:2*n);w=x(2*n+1:end);
% first order variables & parameters
V1=V(1:n1);
D1=D(1:n1);T1=T(1:n1);K1=K(1:n1);
% % second order variables & parameters
V2=V(n1+1:end);
D2=D(n1+1:end);T2=T(n1+1:end);K2=K(n1+1:end);
%% power flow
% power flow complex
U=V.*exp(1i*the);
I=1i*B*U;S=diag(conj(I))*U;
P=real(S);
Bc=zeros(n,n);
for i=1:n
    for j=1:n
    Bc(i,j)=B(i,j)*cos(the(i)-the(j));
    end
end
QV=-Bc*V;
%% bus dynamics
% ddel=zeros(10,1);dw=ddel;dEqq=dw;
% First order dynamics
for i=1:n1
    dthe1(i)=(-P(i)+Pg(i))/D1(i);
    dV1(i)=(-K1(i)*V1(i)-QV(i)+u(i))/T1(i);
end
% Second order dynamics
for i=1:n2
    dthe2(i)=w(i);
    dw(i)=(-D2(i)*w(i)-P(n1+i)+Pg(n1+i))/M(i);%\omega i
    dV2(i)=(-K2(i)*V2(i)-QV(n1+i)+u(n1+i))/T2(i);
end

dx=[dthe1,dthe2,dV1,dV2,dw];
out=dx';

