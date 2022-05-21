clear all;
clc
%'---Chuong trinh giai bai toan TPT cua DTR 9M-22M-10/12/05--'
%'Cac tham so ban dau - cua file solieu 9m22m'
  buocTP=0.01; spt=3;
  dn=0.104;dt=0.020;ltp=0.902; delt=1580; 
  kn=1.25;f0=850000; u1=0.12*10^(-6);nuy=0.708;
  slf=7; dth=0.013; dk=0.114;lk=0.934;
  kt=0.002; Tbd=293; Tbdc=293;kvkhi=0.000006;
  a=0.3; b=5.0;phi2=0.92;
  %'ket thuc cac tham so ban dau'
    s=pi*(dn+dt)*ltp;
    fth=slf*pi*dth^2/4;
    smd=pi/4*(dn^2-dt^2);
    om=smd*ltp*delt;
    fk=pi/4*dk^2;
    f1tbd=1/(1-kt*(Tbd-Tbdc));
    k0k=sqrt(2*kn/(kn+1))*(2/(kn+1))^(1/(kn-1));
    w0=pi/4*dk^2*lk; 
    tgchay=0.0;
  x(1)=0.0001;  %x(1)=si
  x(2)=0.7015; %x(2)=khi
  x(3)=40*10^5;%x(3)=p
%for k=1:100 
k=1;
while x(1)<=1.0     
   %' --Procedure runkut-- '   
    w1=[0 0.5 0.5 1];
for ir=1:spt     %1
    sr(1,ir)=0;
end;              %1
for jr=1:4   %2
    nr=jr+1;
    for kr=1:spt  %3
        y(kr)=x(kr)+w1(jr)*sr(jr,kr);
    end;          %3    
   %' -khai bao ve phai HPT- '
    %y(1)=si; y(2)=khi; y(3)=p;
    ec=1/4*(dn-dt)*y(1);
    ftd=fk-pi/4*((dn-2*ec)^2-(dt+2*ec)^2);
    fp=u1*y(3)^nuy;
    vkhi=phi2*k0k*fth*sqrt(y(2)*f0)/ftd;
    fvkhi=1+kvkhi*vkhi^2;   
    w=w0-om*(1-y(1))/delt;
    u=fp*f1tbd*fvkhi;
      h(1)=s*u*delt/om;  
      h(2)=a*b*h(1)/(1+b*y(1))^2;
      khi1=h(2)/y(2);
      ht1=(phi2*k0k*fth*sqrt(y(2)*f0)+s*u-w*khi1)*y(3);
      ht2=s*u*y(2)*f0*delt;
      h(3)=-(ht1-ht2)/w;
    %' -ket thuc khai bao VP- '
    for lr=1:spt  %4
        sr(nr,lr)=buocTP*h(lr);
    end;          %4
end          %2
    for lr=1:spt  %5
        x(lr)=x(lr)+(sr(2,lr)+2*sr(3,lr)+2*sr(4,lr)+sr(5,lr))/6;
    end;          %5    
    %'--Ket thuc Procedure runkut--'
      tg(k)=tgchay;
      si(k)=x(1);
      khi(k)=x(2);
      p(k)=x(3);    
    tgchay=tgchay+buocTP;
     k=k+1;
end %of while    
plot(tg,p);
title('QUAN HE p=f(t): DTR 9M-22M');
xlabel('t[s]');
ylabel('p[Pa]');
grid;

