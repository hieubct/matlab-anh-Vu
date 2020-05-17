clc;
clear;
M=500;% mot
N=1024;
stdnoise = 1.2;
f=1024;% mot giay lay dc 1024 mau.
% t?o chu?i tín hi?u theo th?i gian
t = ((0:N-1)/f)';
X1= sin(40*pi*t+20);% tín hi?u th? nh?t
X2=0.5*sin(100*pi*t+40);
X3=0.7*sin(98*pi*t+50);
X4=0.8*sin(160*pi*t+60);
X5=X1+X2+X3+X4;
stdnoise = 0.2;
noise = stdnoise*randn(N,5);% nhi?u
X1 = X1 + noise(:,1);% tín hi?u có nhi?u
X2 = X2 + noise(:,2);
X3 = X3 + noise(:,3);
X4 = X4 + noise(:,4);
X5 = X5 + noise(:,5);
% tính các ph??ng sai
X1 = X1 - mean(X1);
X1 = X2 - mean(X2);
X3 = X3 - mean(X3);
X4 = X4 - mean(X4);
X5 = X5 - mean(X5);
% chia cho ?? l?ch chu?n
X1 = X1/std(X1);
X2 = X2/std(X2);
X3 = X3/std(X3);
X4 = X4/std(X4);
X5 = X5/std(X5);
X = [X1 X2 X3 X4];
% v? ?? th? các tín hi?u 
figure(1);
    clf;
    subplot(5,1,1);
    plot(t,X1);
   subplot(5,1,2);
    plot(t,X2);
    subplot(5,1,3);
    plot(t,X3);
    subplot(5,1,4);
    plot(t,X4);
    subplot(5,1,5);
    plot(t,X5);
 % tính ma tr?n hi?p ph??ng sai
 Y1=zeros(N-M+1,M);
 Y2=zeros(N-M+1,M);
 Y3=zeros(N-M+1,M);
 Y4=zeros(N-M+1,M);
 for i=1:M
     Y1(:,i)=X1((1:N-M+1)+i-1);
     Y2(:,i)=X2((1:N-M+1)+i-1);
     Y3(:,i)=X3((1:N-M+1)+i-1);
     Y4(:,i)=X4((1:N-M+1)+i-1);
 end
 Y=[Y1 Y2 Y3 Y4];
 C=Y'*Y/(N-M+1);
 % ma tran rieng vec to rieng
 [vtr gtr]=eig(C);
 gtr=diag(gtr);
 [gtr cs]=sort(gtr,'descend');
 gtrd=gtr(1:100);
 vtr = vtr(:,cs);
 % ve do thi vecto rieng
 figure(3)
plot(gtrd,'*');
%figure(4);
%subplot(3,1,2);
%plot(vtr(:,1:2),'.');
%legend('hang 1','hang 2');
%subplot(3,1,3);
%plot(vtr(:,3:4));
%legend('hang 3','hang 4');