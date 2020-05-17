%function [Fx,Fy]=FDMcircular(X,Y,X_dot,Y_dot)
clc
clear
X=0.0001;
Y=0.0002;
X_dot=0.3;
Y_dot=0.4;
Fx=0;
Fy=0;
epsilon=0.6;%???=0.6
L=0.03;
d=0.06;%????d
lambda=L/d;
m=40;
n=30;
phi1=0;
phi2=2*pi;
delta_phi=(phi2-phi1)/m;
delta_lambda=2/n;
k=1;
P0=zeros(m+1,n+1);% zeros???????????0??
A=zeros(m+1,n+1);
B=zeros(m+1,n+1);
C=zeros(m+1,n+1);
D=zeros(m+1,n+1);
E=zeros(m+1,n+1);
F=zeros(m+1,n+1);
H=zeros(m+1,n+1);
for i=1:1:m+1
    % theta(i)=(i-1)*delta_phi;
       % end
        for j=1:1:n+1
            H(i,j)=1+epsilon*cos((i-1)*delta_phi);
        end 
end
S=0;
T=0;
ERR=1e-3;
GAP=1;
while GAP>ERR
    k=k+1;
    for i=1:1:m+1
        for j=1:1:n+1
            if (i==1)|(j==1)|(i==m+1)|(j==n+1)
                P(i,j)=0;
            else
                A(i,j)=(1+epsilon*cos((i+1/2-1)*delta_phi))^3;
                B(i,j)=(1+epsilon*cos((i-1/2-1)*delta_phi))^3;
                C(i,j)=(d/L)^2*(delta_phi/delta_lambda)^2*(1+epsilon*cos((i-1)*delta_phi))^3;
                D(i,j)=(d/L)^2*(delta_phi/delta_lambda)^2*(1+epsilon*cos((i-1)*delta_phi))^3;
                E(i,j)=A(i,j)+B(i,j)+C(i,j)+D(i,j);
                F(i,j)=6*delta_phi*((1+epsilon*cos((i+1/2-1)*delta_phi))-(1+epsilon*cos((i-1/2-1)*delta_phi)))...
                +12*(delta_phi)^2*(X_dot*cos(phi1+(i-1)*delta_phi)+Y_dot*sin(phi1+(i-1)*delta_phi));
                %P(i,j)=(A(i,j)*P0(i+1,j)+B(i,j)*P0(i-1,j)+C(i,j)*P0(i,j+1)+D(i,j)*P0(i,j-1)-F(i,j))/E(i,j);
                P(i,j)=(A(i,j)*P0(i+1,j)+B(i,j)*P0(i-1,j)+C(i,j)*P0(i,j+1)+D(i,j)*P0(i,j-1)-F(i,j))/E(i,j); 
                if P(i,j)<0
                    P(i,j)=0;
                else
                end
            end
        end
    end
% for i=2:1:m
 %       for j=2:1:n
            S=sum(sum(abs(P-P0)));
            T=sum(sum(abs(P)));
%         end
%     end
    GAP=S/T;
    P0=P;
end
for i=1:1:m+1
    for j=1:1:n+1
        Fx=Fx+P(i,j)*cos(phi1+(i-1)*delta_phi)*delta_phi*delta_lambda;
        Fy=Fy+P(i,j)*sin(phi1+(i-1)*delta_phi)*delta_phi*delta_lambda;
    end
end