function J=Jacnbus(t,x,Bbus0,M,d,T,K,~,~)
% input: columns
n2=length(M);n=length(d);n1=n-n2;
the=x(1:n);V=x(n+1:2*n);w=x(2*n+1:end);
Bzz=func_Bzz(the,V,Bbus0);
A=Bzz(1:n,1:n);D=Bzz(1:n,n+1:2*n);C=Bzz(n+1:2*n,n+1:2*n);
temp=-[A(1:n1,1:n1),A(1:n1,n1+1:n),D(1:n1,:),zeros(n1,n2);...
       zeros(n2,2*n),-eye(n2);...
       D(1:n1,:)',D(n1+1:n,:)',diag(K)+C,zeros(n,n2);...
       A(n1+1:n,1:n1),A(n1+1:n,n1+1:n),D(n1+1:n,:),diag(d(n1+1:n))];
J=diag([1./d(1:n1);ones(n2,1);1./T;1./M])*temp;
end