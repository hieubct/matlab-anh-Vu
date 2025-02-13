% Chuong trinh giai bai toan thuat phong trong, xac dinh ap suat lon nhat
% Ha Noi, 13/4/2011

clc;
clear all;
%--------------------------------------------------------------------------
%du lieu ban dau
  %mCd=18;
  n=1;                 %so thanh thuoc phong (n).
  %cac dac trung cua buong dot
  Dk=0.112;            %duong kinh trong buong dot.
  %simaB=6300*10^5;
  %roK=7800;
  
  %cac dac trung vat ly cua thuoc phong
  roT=1580;            %khoi luong rieng thuoc phong.
  u1=0.12*10.^-6;
  vChay=0.69;
  a=0.3;
  b=5;
  k=1.25;
  KT=0.011;
  f0=850000;
  
  %cac thong so ban dau.
  %pGh=5*10^6;
  pMoi=4*10^6;
  khin0=0.7015;
  phi2=0.92;
  kV=6*10^-6;
  TbdMin=223;
  Tbd=323;
  Tbdc=288.9;
  h=0.01;              %buoc tich phan (chu y: khong duoc de buoc tich phan
                        %qua nho ket qua ra dang bang se tran man hinh).
  pa=1*10.^5;           %ap suat moi truong .
  t=[];
  p=[];
  psi=[];
  khi=[];
  t(1)=0;
  p(1)=pMoi;                                %ap suat ban dau = ap suat moi.
  khi(1)=khin0;

  
  %cac dai luong tinh toan
  H=zeros(1,3);         %ve phai he phuong trinh vi phan.
  DH=zeros(5,3);        %cac he so k.
  KQtg=zeros(1,3);      %ket qua trung gian tuong ung psi, khi va ap suat.
  hs=[0,0.5,0.5,1];     %he so.
  KQ=zeros(1,3);        %ket qua tuong ung psi, khi va ap suat.
  
  %cac dac trung ket cau lieu phong
  mT=7.57;                 %khoi luong thuoc phong wT(kg).
  Fk=pi/4*Dk^2;
  K0=(2/(1+k))^(1/(k-1))*sqrt(2*k/(1+k));   %ham so mu doan nhiet.
  f1=1/(1-KT*(Tbd-Tbdc));
  f1min=1/(1-KT*(TbdMin-Tbdc));
  C=phi2*K0/sqrt(f0);
  
%tinh khoi luong ket cau nho nhat khi duong kinh ngoai thanh thuoc(Dn) thay doi
  Dn=0.101;
  dt=0.02;
  L=0.902;           %chieu dai thanh thuoc phong.
  Lk=0.99;                                %chieu dai buong dot.
  Vk=Lk*Fk;                                %the tich buong dot.
  Vtp=mT/roT;
  S=n*pi*(Dn*L+L*dt);                      %dien tich chay ban dau.
  %Fth=roT*u1*f1min*S/(C*pGh^(1-vChay));    %dien tich toi han cua loa phut.
  Fth=0.018;
  psi0=pMoi*(Vk-Vtp)/(mT*f0-pMoi*Vtp);
  psi(1)=psi0;

  dem=2;
  KQ(1)=psi(1);
  KQ(2)=khi(1);
  KQ(3)=p(1);
  %giai he TPT bang phuong phap runge-kutty
  while psi(dem-1)<=1
  for i=1:3
    DH(1,i)=0;
  end;
  for i=1:4
    for j=1:3
      KQtg(j)=KQ(j)+hs(i)*DH(i,j); 
    end;
    ec=1/4*(Dn-dt)*KQtg(1);
    Ftd=Fk-pi/4*((Dn-2*ec)^2-(dt+2*ec)^2);
    %Ftd=Fk-mT*(1-KQtg(1))/roT/L;          %LA CACH TINH Ftd KHAC
        
    vKhi=phi2*K0*Fth*sqrt(KQtg(2)*f0)/Ftd; %van toc khi thuoc
    phiV=1+kV*vKhi^2;                      %ham anh huong boi van toc khi
    %van toc chay
    u=u1*(KQtg(3))^vChay*f1*phiV;
    %cac he so k
    DH(i+1,1)=h*S*u*roT/mT;
    DH(i+1,2)=a*b*DH(i+1,1)/(1+b*KQtg(1))^2;
    khi1=DH(i+1,2)/(h*KQtg(2));
    V=Vk-mT*(1-KQtg(1))/roT;
    
    DH(i+1,3)=-h*( ( phi2*K0*Fth*sqrt(f0*KQtg(2))+S*u-V*khi1 )*KQtg(3)-S*u*KQtg(2)*roT*f0 )/V;               
  end;
  %thiet lap ket qua
  for i=1:3               
      KQ(i)=KQ(i)+( DH(2,i)+2*DH(3,i)+2*DH(4,i)+DH(5,i) )/6;
  end;                 
  psi(dem)=KQ(1);
  khi(dem)=KQ(2);
  p(dem)=KQ(3);
  t(dem)=t(dem-1)+h;
  
  dem=dem+1;
end;
  dem=dem-1;
  demTk=dem;
  
%-------------------------------------------------
% Code giai bai toan thuat phong trong giai doan phut khi tu do
%giai doan chay het thuoc phong
    %ap suat tai thoi diem chay het
    pk=p(dem);
    %khoi luong khi thuoc trong buong dot
    mk=Vk*pk/f0;
    %the tich rieng cua khi trong buong dot
    vk=f0/pk;
    %he so b1
    b1=phi2*k*K0*Fth*sqrt(pk/vk)/mk;
    B=(k-1)*b1/(2*k);
%ap suat y
while p(dem)>1.75*10.^5
    dem=dem+1;
    t(dem)=t(dem-1)+h; 
    p(dem)=pk*( 1+B*t(dem) )^(-2*k/(k-1));
end;
dem=dem-1;
%--------------------------------------------------
%tim ap suat lon nhat
  imax=2;
  for i=3:dem
    if p(imax)<p(i)
        imax=i;
    end;
end;

tk=t(demTk);                              %thoi gian chay cua thuoc phong la
pMax=p(imax);                           %ap suat lon nhat
ptt=1.1*pMax;                           %ap suat tinh toan (Pa).
%--------------------------------------------------
%Ket qua giai bai toan thuat phong trong
tk
ptt
plot(t,p);
xlabel('t(s)');
ylabel('p(Pa)');
grid;
