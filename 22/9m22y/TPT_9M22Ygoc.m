%--------PROGRAM DTR-----------
%'---chuong trinh giai bai toan TPT cua DTR M21-4/10/07--'
clear all;
close all;
clf;
' Cac thong so ban dau- cua file so lieu M21'
    buocTP  =0.0001;
    spt     =5;
    Dn1     =0.1; 
    Dn2     =0.098; 
    dtr1    =0.024;
    dtr2    =0.016;
    Lt1     =0.9;
    Lt2     =0.893;
    Ro1     =1580;
    Ro2     =1580;
    kn      =1.25; 
    f0      =850000; 
    u1      =0.12*10^(-6); 
    nuy     =0.69;
    slf     =7; 
    dth     =0.018; 
    Dk      =0.112; 
    kw      =0.6*10^(-5);
    kt      =0.011; 
    Tbd     =288.9; 
    Tbdc    =273; 
    a       =0.3; 
    b       =5.0; 
    phi2    =0.92;
    pmoi    =40*10^5;
%'Ket thuc cac tham so mau ban dau'
 
%Bat dau tinh toan trong chuong trinh
    Lk=(Lt1+Lt2)*1.1;
    Fth=slf*pi*dth^2/4;
    Smd1=pi/4*(Dn1^2-dtr1^2);
    Smd2=pi/4*(Dn2^2-dtr2^2);
    Om1=Smd1*Lt1*Ro1;
    Om2=Smd2*Lt2*Ro2;
    Om=Om1+Om2;
    Fk=pi/4*Dk^2;
    Vk=pi*Dk^2*Lk/4;
    ft=1/(1-kt*(Tbd-Tbdc));
    k0=sqrt(2*kn/(kn+1))*(2/(kn+1))^(1/(kn-1));
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
    
  while x(3)<=1.0 % Of vong si
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
     
      S1=(pi*(Dn1+dtr1)*Lt1);
      S2=(pi*(Dn2+dtr2)*Lt2);
      fp1=u1*y(5)^nuy;
      fp2=u1*y(5)^nuy;
      et1=1/4*(Dn1-dtr1)*y(1);
      et2=1/4*(Dn2-dtr2)*y(2);
      ftd=Fk-pi/4*((Dn1-2*et1)^2-(dtr1+2*et1)^2+(Dn2-2*et2)^2-(dtr2+2*et2)^2)/2;
      w=phi2*k0*Fth*sqrt(y(4)*f0)/ftd;
      fw=1+kw*w^2;
      uc1=fp1*ft*fw;
      uc2=fp2*ft*fw;
      if y(1)>=1 
          y(1)=1;
          h(1)=0; 
      else h(1)=S1*uc1*Ro1/Om1;
      end
      if y(2)>=1
          y(2)=1;
          h(2)=0; 
      else h(2)=S2*uc2*Ro2/Om2;
      end          
      h(3)=(Om1*h(1)+Om2*h(2))/Om;
      h(4)=a*b*h(3)/(1+b*y(3))^2;
      khi1=h(4)/y(4);
      V=Vk-(Om1*(1-y(1))/Ro1+Om2*(1-y(2))/Ro2);
      ht1=(phi2*k0*Fth*sqrt(y(4)*f0)+S1*uc1+S2*uc2-V*khi1)*y(5);
      ht2=(S1*uc1*Ro1+S2*uc2*Ro2)*y(4)*f0;
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
  
  tk=tgchay-buocTP; %thoi gian chay cua tich phan
  kk=k-1;
  k=kk;
  pk=x(5);
  mk=pk*Vk/f0; 
  vkk=Vk/mk;
  bnk=(phi2*k0*kn*Fth*sqrt(pk/vkk))/mk;
  btk=(kn-1)*bnk/2/kn;
  k=kk+1;
  tchk=buocTP;
  while x(5)>=1*10^5 %cua thoi ki phut khi tu do
    x(5)=pk*(1+btk*tchk)^(-(2*kn/(kn-1)));
    p(k)=x(5);
	tchk=tchk+buocTP;
    tg(k)=tk+tchk;
    k=k+1;
  end %end cua thoi ki phut khi tu do
  tc=tk+tchk-buocTP; %thoi gian lam viec cua DC
  hold on;
  plot(tg,p); %Ve do thi p(t)
  title('QUAN HE p=f(t) DAN 9M22Y');
  xlabel('t[s]');
  ylabel('p[Pa]');
  grid on;
  P=1.5*Fth*p; 
  Pmax=1.5*Fth*pmax;
  disp('----------------------------------------------');
  disp('  Ket qua chay chuong trinh thu duoc:');
  fprintf('   + Luc day lon nhat Pmax= %g (N) \n',Pmax);
  fprintf('   + Ap suat lon nhat pmax= %g (MPa) \n',pmax/1e6);
  fprintf('   + Thoi gian lam viec cua dong co %g (s) \n',tc);
  fprintf('   + Thoi gian thuoc phong chay het  %g (s) \n',tk);
  hold on;
  plot(tg,P); %Ve do thi P(t)
  title('QUAN HE P=g(t) DAN 9M22Y');
  xlabel('t[s]');
  ylabel('P[N]');
  grid on;

  
  