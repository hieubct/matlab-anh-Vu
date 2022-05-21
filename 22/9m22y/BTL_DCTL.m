%du lieu ban dau
  %cac dac trung ket cau
  Dk=0.112;
  Lk=0.99;
  dth=7*0.018;
  D=0.101;
  Di=[];        %mang duong kinh ngoai D
  Di(1)=D;
  d2=0.02;
  d1=d2;
  d1i=[];        %mang duong kinh trong d1
  d1i(1)=d1;
 
  d2i=[];        %mang duong kinh trong d2
  d2i(1)=d2;
  L=0.902;
  L1=300;

  %cac dac trung vat ly
  roT=1580;     %mat do thuoc phong
  u1=0.12*10.^-6;
  v=0.69;
  k=1.25;
  KT=0.011;
  f0=850000;
  Tbd=288.9;
  
  %cac thong so khac
  p=[];
  p(1)=3*10.^6;           %ap suat ban dau = ap suat moi
  pa=1*10.^5;         %ap suat moi truong 
  a=0.3;
  b=5;
  psi0=0.001;
  psi=[];
  psi(1)=psi0;
  khin0=0.7;
  khi=[];
  khi(1)=khin0;
  phi2=0.98;
  h=0.01;               %buoc tich phan (chu y: khong duoc de buoc tich phan 
                        %qua nho ket qua ra dang bang se tran man hinh)
  
  %cac dai luong tinh toan
  H=zeros(1,3);         %ve phai he phuong trinh vi phan
  DH=zeros(5,3);        %cac he so k
  KQtg=zeros(1,3);      %ket qua trung gian tuong ung psi, khi va ap suat
  hs=[0,0.5,0.5,1];     %he so
  KQ=zeros(1,3);        %ket qua tuong ung psi, khi va ap suat
  
%xu ly so lieu ban dau
  Vk=Lk*(pi*Dk*Dk/4); %the tich buong dot
  f1=1/(1-KT*(Tbd-288));
  K0=(2/(1+k)).^(1/(k-1))*sqrt(2*k/(1+k));  %ham so mu doan nhiet
  Fth=pi*dth*dth/4;                         %dien tich toi han
  S0=pi*(D*L+L1*d1+(L-L1)*d2+(d1*d1-d2*d2)/4);     %dien tich chay ban dau
  Vtp=L1*(pi*(D*D-d1*d1)/4)+(L-L1)*(pi*(D*D-d2*d2)/4);         %the tich thuoc phong
  m0=roT*Vtp;                               %khoi luong thuoc phong ban dau
   
%buoc 1 t0=0;  
  KQ(1)=psi0;
  KQ(2)=khin0;
  KQ(3)=p(1);
  
dem=2;
while D>=d1
  for i=1:3
    DH(1,i)=0;
  end;
  
  for i=1:4
    for j=1:3
      KQtg(j)=KQ(j)+hs(i)*DH(i,j); 
    end;
    %van toc chay
    u=u1*f1*(KQtg(3)).^v;
    %dien tich chay
    S=pi*(D*L+L1*d1)+(L-L1)*d2+(d1*d1-d2*d2)/4;     %dien tich chay 
    %cac he so k
    i1=i+1;
    DH(i1,1)=h*S*u*roT/m0;
    DH(i1,2)=a*b*DH(i1,1)/(1+b*KQtg(1)).^2;
    khi1=DH(i1,2)/(h*KQtg(2));
    V=Vk-L1*(pi*(D*D-d1*d1)/4)-(L-L1)*(pi*(D*D-d2*d2)/4);
    DH(i1,3)=-h*( ( phi2*K0*Fth*sqrt(f0*KQtg(2))+S*u-V*khi1 )*KQtg(3)-S*u*KQtg(2)*roT*f0 )/V;               
  end;
  %thiet lap ket qua
  for i=1:3               
      KQ(i)=KQ(i)+( DH(2,i)+2*DH(3,i)+2*DH(4,i)+DH(5,i) )/6;
  end;  
                
  u2=u1*f1*(KQ(3)).^v;
  D=D-2*h*u2;
  d1=d1+2*h*u2;
  d2=d2+2*h*u2; 
  %L=L-2*h*u2;
  
  psi(dem)=KQ(1);
  khi(dem)=KQ(2);
  p(dem)=KQ(3);
  Di(dem)=D;
  d1i(dem)=d1;
  d2i(dem)=d2;
  
  dem=dem+1;
end;         

L=L-L1;

while D>=d2
  for i=1:3
    DH(1,i)=0;
  end;
  
  for i=1:4
    for j=1:3
      KQtg(j)=KQ(j)+hs(i)*DH(i,j); 
    end;
    %van toc chay
    u=u1*f1*(KQtg(3)).^v;
    %dien tich chay
    S=pi*L*(D+d2);
    %cac he so k
    i1=i+1;
    DH(i1,1)=h*S*u*roT/m0;
    DH(i1,2)=a*b*DH(i1,1)/(1+b*KQtg(1)).^2;
    khi1=DH(i1,2)/(h*KQtg(2));
    V=Vk-L*(pi*(D*D-d2*d2)/4);
    DH(i1,3)=-h*( ( phi2*K0*Fth*sqrt(f0*KQtg(2))+S*u-V*khi1 )*KQtg(3)-S*u*KQtg(2)*roT*f0 )/V;
    
  end;
  %thiet lap ket qua
  for i=1:3
      KQ(i)=KQ(i)+( DH(1,i)+2*DH(2,i)+2*DH(3,i)+DH(4,i) )/6;
  end;                                  
  u2=u1*f1*(KQ(3)).^v;
  D=D-2*h*u2;
  d2=d2+2*h*u2;  
  %L=L-2*h*u2;
  
  psi(dem)=KQ(1);
  khi(dem)=KQ(2);
  p(dem)=KQ(3);  
  Di(dem)=D;
  d1i(dem)=0;
  d2i(dem)=d2;
  
  dem=dem+1;
end;

dem=dem-1;
dem2=dem;
pk=p(dem);

%giai doan chay het thuoc phong
    %khoi luong khi thuoc trong buong dot
    mk=Vk*pk/f0;
    %the tich rieng cua khi trong buong dot
    vk=f0/pk;
    %he so b1
    b1=phi2*k*K0*Fth*sqrt(pk/vk)/mk;
    B=(k-1)*b1/(2*k);
    
%ap suat y
x(1)=0;
y(1)=p(1);
P(1)=phi2*K0*Fth*y(1)*sqrt(2*k*(1-(pa/y(1)).^((k-1)/k))/(k-1));
j=2;
while y(j-1)>1.75*10.^5
    if j<dem
        x(j)=(j-1)*h;
        y(j)=p(j);
    else
       %giai doan chay het thuoc phong
       x(j)=(j-1)*h; 
       y(j)=pk*( 1+B*(x(j)-(dem-1)*h ) ).^(-2*k/(k-1));
    end;
    %luc day o che do tinh toan
    P(j)=phi2*K0*Fth*y(j)*sqrt(2*k*(1-(pa/y(j)).^((k-1)/k))/(k-1));
    j=j+1;
end;
j=j-1;

plot(x,y);
xlabel('Thoi gian: t(s)');
ylabel('Ap suat buong dot:p(pa)');
title('Do Thi Moi Quan He Ap Suat Buong Dot Va Thoi Gian');
grid on;

plot(x,P);
xlabel('Thoi gian: t(s)');
ylabel('Luc day dong co:P(N)');
title('Do Thi Moi Quan He Luc Day Dong Co Va Thoi Gian');
grid on;

plotyy(x,y,x,P);
xlabel('Thoi gian: t(s)');
ylabel('Ap suat buong dot:p(pa)');
title('Do Thi Moi Quan He p(pa),P(N) Voi Thoi Gian');
text(0.3,8.1*10.^6,'Ap suat:p(pa)');
text(0.5,6.4*10.^6,'Luc day:P(N)');
grid on;

imax=1;
for i=2:j
    if y(imax)<y(i)
        imax=i;
    end;
end;

%thoi gian chay la
tk=dem2*h;
%ap suat tai thoi diem thuoc phong chay het
pk;
imax;
%ap suat lon nhat
y(imax);
%luc day lon nhat
P(imax);

