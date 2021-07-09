%--------PROGRAM DTR-----------
%'---chuong trinh giai bai toan TPT cua DTR M21-4/10/07--'
clear all;
close all;
clf;
clc;
%' Cac thong so ban dau- cua file so lieu M21'
buocTP  =0.0001;
spt     =5;% ch? so ma tran ket qua
Dn     =0.0048;% ?uòng kính ngoai thuoc phong
dtr    =0.0014;% d??ng kinh trong thuoc phong
Lt     =0.085; %chieu dai thu?c phong
Ro     =1540;% kh?i l??ng riêng c?a TP
n      =19;% so luong thanh thuoc phong
kn      =1.225;
f0      =960000;
u1      =0.6*10^(-3); %0.6*10^(-6);
nuy     =0.8;
slf     =1;%7;
dth     =0.003;
Dk      =0.025;
kw      =0.6*10^(-5);% he so trong ham xoi mon
kt      =0.011;% he so nhiet do
Tbd     =293; % nhiet do ban dau
Tbdc    =273;   % nhiet do ban dau chuan
a       =0.3; % hs thuc nghiem tinh he so ton that nhiet
b       =5.0; % hs thuc nghiem tinh he so ton that nhiet
phi2    =0.97;
pmoi    =10*10^6;
%'Ket thuc cac tham so mau ban dau'
%Bat dau tinh toan trong chuong trinh
Lk=Lt*1.1;%chieu dai buong dot
Fth=slf*pi*dth^2/4;% dien tich toi han
Smd=n*pi/4*(Dn^2-dtr^2);%di?n tích m?t ??u c?a thu?c phóng
Om=Smd*Lt*Ro;%khoi luong thuoc phong
Fk=pi/4*Dk^2;%dien tich mat cat ngang buong dot
Vk=pi*Dk^2*Lk/4;% the tich dong co
ft=1/(1-kt*(Tbd-Tbdc));%ham fu thuoc nhiet do ban dau
k0=sqrt(2*kn/(kn+1))*(2/(kn+1))^(1/(kn-1));%hàm cua chi so mu doan nhiet
tgchay=0.0;
x(1)=0.00001;
x(2)=0.00001;
x(3)=0.00001;
x(4)=0.7105;
x(5)=pmoi;
tg(1)=tgchay;
si(1)=x(3);
khi(1)=x(4);
p(1)=x(5);
k=2;
pmax=0.0;


while x(3)<1.0 % Of vong si
    %---Procedure runkut-------
    w1=[0 0.5 0.5 1];
    for ir=1:spt
        sr(1,ir)=0;
    end;
    for jr=1:4
        nr=jr+1;
        for kr=1:spt
            y(kr)=x(kr)+w1(jr)*sr(jr,kr);
        end;
        %'-khai bao ve phai HPT-'
        %y(1)=si1; y(2)=si2; y(3)=si; y(4)=khi; y(5)=p;
        S=n*(pi*(Dn+dtr)*Lt);%dien dich be matchay
        fp=u1*y(5)^nuy;% toc do chay cua thuoc phong
        et=1/4*(Dn-dtr)*y(1);% 1/2 b? dày thuoc da cháy
        ftd=Fk-pi/4*((Dn-2*et)^2-(dtr+2*et)^2);% di?n tich thoatkhi tai thoi diem t
        w=phi2*k0*Fth*sqrt(y(4)*f0)/ftd;%
        fw=1+kw*w^2;
        uc=fp*ft*fw;
        if y(1)>=1
            y(1)=1;
            h(1)=0;
        else h(1)=S*uc*Ro/Om
        end
        h(2)=0;
        h(3)=(Om*h(1))/Om;
        h(4)=a*b*h(3)/(1+b*y(3))^2;
        khi1=h(4)/y(4);
        V=Vk-Om*(1-y(1))/Ro;%th? tích t? do
        ht1=(phi2*k0*Fth*sqrt(y(4)*f0)+S*uc-V*khi1)*y(5);
        ht2=(S*uc*Ro)*y(4)*f0;
        h(5)=-(ht1-ht2)/V;
        %'- ket thuc khai bao VP-'
        for lr=1:spt
            sr(nr,lr)=buocTP*h(lr);
        end;
    end
    for lr=1:spt
        x(lr)=x(lr)+(sr(2,lr)+2*sr(3,lr)+2*sr(4,lr)+sr(5,lr))/6;
    end;
    %'-ket thuc Procedure runkut-'
    tg(k)=tgchay;
    si(k)=x(3);
    khi(k)=x(4);
    p(k)=x(5);
    if pmax<=x(5);
        pmax=x(5);%tinh pmax
    end
    tgchay=tgchay+buocTP;
    k=k+1;
end %End cua vong si







hold on;
plot(tg,si); %Ve do thi P(t)
title('QUAN HE Si=g(t) DAN 9M22Y');
xlabel('t[s]');
ylabel('Si');
grid on;



%   tk=tgchay-buocTP; %thoi gian chay cua tich phan
%   kk=k-1;
%   k=kk;
%   pk=x(5);
%   mk=pk*Vk/f0;
%   vkk=Vk/mk;
%   bnk=(phi2*k0*kn*Fth*sqrt(pk/vkk))/mk;
%   btk=(kn-1)*bnk/2/kn;
%   k=kk+1;
%   tchk=buocTP;
%   while x(5)>=1*10^5 %cua thoi ki phut khi tu do
%     x(5)=pk*(1+btk*tchk)^(-(2*kn/(kn-1)));
%     p(k)=x(5);
% 	tchk=tchk+buocTP;
%     tg(k)=tk+tchk;
%     k=k+1;
%   end %end cua thoi ki phut khi tu do
%
%
%
%
%
%
%   tc=tk+tchk-buocTP; %thoi gian lam viec cua DC
%   hold on;
%   plot(tg,p); %Ve do thi p(t)
%   title('QUAN HE p=f(t) DAN 9M22Y');
%   xlabel('t[s]');
%   ylabel('p[Pa]');
%   grid on;
%   P=1.5*Fth*p;
%   Pmax=1.5*Fth*pmax;
%   disp('----------------------------------------------');
%   disp('  Ket qua chay chuong trinh thu duoc:');
%   fprintf('   + Luc day lon nhat Pmax= %g ? \n',Pmax);
%   fprintf('   + Ap suat lon nhat pmax= %g (MPa) \n',pmax/1e6);
%   fprintf('   + Thoi gian lam viec cua dong co %g (s) \n',tc);
%   fprintf('   + Thoi gian thuoc phong chay het  %g (s) \n',tk);
%   hold on;
%   plot(tg,P); %Ve do thi P(t)
%   title('QUAN HE P=g(t) DAN 9M22Y');
%   xlabel('t[s]');
%   ylabel('P[N]');
%   grid on;
