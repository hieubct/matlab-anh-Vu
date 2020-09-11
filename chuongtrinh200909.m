clc
clear
format long g
% cac thong so ban dau
R = 100e-3; % ban kinh truc
l = 150e-3; % chieu dai o
Rt=100.21e-3;% ban kinh trong cua pad
Rn=115e-3; % ban kinh ngoai cua pad
Rlr=100.09e-3; % ban kinh lap rap
cl = Rt-R; % khe ho cua o truc ( tam truc va tam pad trung nhau)
cb = Rlr-R; % khe ho cai dat ( tam truc va tam o trung nhau)
nguy = 0.063; % do nhot cua dau
d1 = 2*R; % duong kinh truc
N=2000;
W=4000;
lrz=2*Rn/l;
%posi = c1/R; % ti le khe ho
exilonn1=0;%1.5e-5;% do nang cua pad1
%cbtt1=cb-exilonn1;
exilonn2=0;%7e-6;% do nang cua pad1
%cbtt2=cb-exilonn2;
exilonn3=0;%3e-6;% do nang cua pad1
%cbtt3=cb-exilonn3;
M =0.5714;%
M1 = 0.5714;%1-(cbtt1/cl);
M2 = 0.5714;%1-(cbtt2/cl);
M3 = 0.5714;%1-(cbtt3/cl);
anpha1 = 0;%0.035*pi/180;%0.009*pi/180;% goc nghieng cua pad1
anpha2 = 0.0266*pi/180;% goc nghieng cua pad2
Precess2=3.5;%0.53;%0.5;%0.551;% áp su?t vùng lõm2
anpha3 = 0;% goc nghieng cua pad3
capa1=Rt*anpha1/cl;
capa2=Rt*anpha2/cl;
capa3=Rt*anpha3/cl;
ld = l/d1; % ti le chieu dai tren duong kinh
m = 160; % so luoi theo chu vi
n = 90; % so luoi theo phuong doc truc
%e1 =0.75e-04;% 0.1e-3; % do lech tam ban dau
exilon = 0.15;%0.05;%e1/c1; % ti le lech tam
phibd = 1.0455; % goc lech ban dau cua truc
% thong so pad1
phi11 = 3*pi/180; % goc dau cua pad
phi12 = 117*pi/180;% goc sau cua pad
beta1 = (phi11+phi12)/2; % toa do cua diem xoay
% thong so pad2
phi21 = 123*pi/180; % goc dau cua pad
phi22 = 237*pi/180;% goc sau cua pad
beta2 = (phi22+phi21)/2; % toa do cua diem xoay
% thong so pad3
phi31 = 243*pi/180; % goc dau cua pad
phi32 = 357*pi/180;% goc sau cua pad
beta3 = (phi32+phi31)/2; % toa do cua diem xoay
% tinh cac thong so ban dau
deltalanda = 2/n; % so luoi theo phuong doc truc
deltaphi = (phi22-phi21)/m; % chia luoi theo chu vi
% ma tran he so ban dau cho phan hydrodynamic
p01 = zeros(m+1,n+1);
p1 = zeros(m+1,n+1);
a1 = zeros(m+1,n+1);
b1 = zeros(m+1,n+1);
c1 = zeros(m+1,n+1);
d1 = zeros(m+1,n+1);
e1 = zeros(m+1,n+1);
f1 = zeros(m+1,n+1);
p02 = zeros(m+1,n+1);
p2 = zeros(m+1,n+1);
a2 = zeros(m+1,n+1);
b2 = zeros(m+1,n+1);
c2 = zeros(m+1,n+1);
d2 = zeros(m+1,n+1);
e2 = zeros(m+1,n+1);
f2 = zeros(m+1,n+1);
p03 = zeros(m+1,n+1);
p3 = zeros(m+1,n+1);
a3 = zeros(m+1,n+1);
b3 = zeros(m+1,n+1);
c3 = zeros(m+1,n+1);
d3 = zeros(m+1,n+1);
e3 = zeros(m+1,n+1);
f3 = zeros(m+1,n+1);
% ma tran he so cho phan hydrostatic
p01s = zeros(m+1,n+1);
p1s  = zeros(m+1,n+1);
a1s  = zeros(m+1,n+1);
b1s  = zeros(m+1,n+1);
c1s  = zeros(m+1,n+1);
d1s  = zeros(m+1,n+1);
e1s  = zeros(m+1,n+1);
f1s  = zeros(m+1,n+1);
p02s = zeros(m+1,n+1);
p2s  = zeros(m+1,n+1);
a2s  = zeros(m+1,n+1);
b2s  = zeros(m+1,n+1);
c2s  = zeros(m+1,n+1);
d2s  = zeros(m+1,n+1);
e2s  = zeros(m+1,n+1);
f2s  = zeros(m+1,n+1);
p03s = zeros(m+1,n+1);
p3s  = zeros(m+1,n+1);
a3s  = zeros(m+1,n+1);
b3s  = zeros(m+1,n+1);
c3s  = zeros(m+1,n+1);
d3s  = zeros(m+1,n+1);
e3s  = zeros(m+1,n+1);
f3s  = zeros(m+1,n+1);
% khai bao cac bien ban dau
S1 = 0;
T1 = 0;
ERR1 = 10e-4;
GAP1 =1;
k1 = 1; % he so lap
S2 = 0;
T2 = 0;
ERR2 = 10e-4;
GAP2 =1;
k2 = 1; % he so lap
S3 = 0;
T3= 0;
ERR3 = 10e-4;
GAP3 =1;
k3 = 1; % he so lap
S1s = 0;
T1s = 0;
ERR1s = 10e-4;
GAP1s =1;
k1s = 1; % he so lap
S2s = 0;
T2s = 0;
ERR2s = 10e-4;
GAP2s =1;
k2s = 1; % he so lap
S3s = 0;
T3s = 0;
ERR3s = 10e-4;
GAP3s =1;
k3s = 1; % he so lap
%% TINH CHO PAD1
%  do day mang dau
for i=1:m+1
    for j=1:n+1
        % phan dong
        phi1(i,j) = phi11+(i-1)*deltaphi;
        phic1(i,j)= phi11+(i-1+1/2)*deltaphi;
        phit1(i,j)= phi11+(i-1-1/2)*deltaphi;
        h1(i,j) = 1+exilon*cos(phi1(i,j)-phibd)-M1*cos(phi1(i,j)-beta1)+capa1*sin(phi11-phi1(i,j));
        hc1(i,j) = 1+exilon*cos(phic1(i,j)-phibd)-M1*cos(phic1(i,j)-beta1)+capa1*sin(phi11-phic1(i,j));
        ht1(i,j) = 1+exilon*cos(phit1(i,j)-phibd)-M1*cos(phit1(i,j)-beta1)+capa1*sin(phi11-phit1(i,j));
        % phan tinh
        h1s(i,j)  = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phi1(i,j));
        hc1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phic1(i,j));
        ht1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phit1(i,j));
    end
end
% ch??ng trinh tinh ap suat dong
while GAP1>ERR1
    k1 = k1+1;
    fx1=0;
    fy1=0;
    for i =1:m+1
        for j =1:n+1
            if i==1||i==m+1||j==1||j==n+1
                p1(i,j)=0;
            else
                a1(i,j)= hc1(i,j)^3;
                b1(i,j)= ht1(i,j)^3;
                c1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                d1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                e1(i,j)= a1(i,j)+b1(i,j)+c1(i,j)+d1(i,j);
                f1(i,j)= 3*deltaphi*(hc1(i,j)-ht1(i,j));
                p1(i,j)= (a1(i,j)*p01(i+1,j)+b1(i,j)*p01(i-1,j)+c1(i,j)*p01(i,j+1)... % he phuong trinh tinh ap suat
                    +d1(i,j)*p01(i,j-1)-f1(i,j))/e1(i,j);
                if p1(i,j)<=0
                    p1(i,j)=0;
                end
            end
        end
    end
    % Kiem tra dieu kien hoi tu
    for i= 1:m+1
        for j= 1:n+1
            S1=S1+abs(p1(i,j)-p01(i,j));
            T1=T1+abs(p1(i,j));
            GAP1=S1/T1;
        end
    end
    p01=p1;
end
% chuong trinh tinh ap suat tinh
while GAP1s>ERR1s
    k1s = k1s+1;
    for i =1:m+1
        for j =1:n+1
            if j==1||j==n+1||i==1||i==m+1
                p1s(i,j)=0;
            else
                a1s(i,j)= hc1s(i,j)^3;
                b1s(i,j)= ht1s(i,j)^3;
                c1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
                d1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
                e1s(i,j)= a1s(i,j)+b1s(i,j)+c1s(i,j)+d1s(i,j);
                f1s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                p1s(i,j)= (a1s(i,j)*p01s(i+1,j)+b1s(i,j)*p01s(i-1,j)+c1s(i,j)*p01s(i,j+1)...  % he phuong trinh tinh ap suat
                    +d1s(i,j)*p01s(i,j-1)-f1s(i,j))/e1s(i,j);
                %if p(i,j)<=0
                %p(i,j)=0;
                %end
            end
        end
    end
    for i= 1:m+1
        for j= 1:n+1
            S1s=S1s+abs(p1s(i,j)-p01s(i,j));
            T1s=T1s+abs(p1s(i,j));
            GAP1s=S1s/T1s;
        end
    end
    p01s=p1s;
end
% Can bang luu luong


%% TINH CHO PAD2
%  do day mang dau
for i=1:m+1
    for j=1:n+1
        % phan dong
        phi2(i,j) = phi21+(i-1)*deltaphi;
        phic2(i,j)= phi21+(i-1+1/2)*deltaphi;
        phit2(i,j)= phi21+(i-1-1/2)*deltaphi;
        h2(i,j) = 1+exilon*cos(phi2(i,j)-phibd)-M2*cos(phi2(i,j)-beta2)+capa2*sin(phi21-phi2(i,j));
        hc2(i,j) = 1+exilon*cos(phic2(i,j)-phibd)-M2*cos(phic2(i,j)-beta2)+capa2*sin(phi21-phic2(i,j));
        ht2(i,j) = 1+exilon*cos(phit2(i,j)-phibd)-M2*cos(phit2(i,j)-beta2)+capa2*sin(phi21-phit2(i,j));
        % phan tinh
        h2s(i,j)  = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phi2(i,j));
        hc2s(i,j) = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phic2(i,j));
        ht2s(i,j) = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phit2(i,j));
    end
end
% ch??ng trinh tinh ap suat dong
while GAP2>ERR2
    k2 = k2+1;
    fx1=0;
    fy1=0;
    for i =1:m+1
        for j =1:n+1
            if i==1||i==m+1||j==1||j==n+1
                p2(i,j)=0;
            else
                a2(i,j)= hc2(i,j)^3;
                b2(i,j)= ht2(i,j)^3;
                c2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
                d2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
                e2(i,j)= a2(i,j)+b2(i,j)+c2(i,j)+d2(i,j);
                f2(i,j)= 3*deltaphi*(hc2(i,j)-ht2(i,j));
                p2(i,j)= (a2(i,j)*p02(i+1,j)+b2(i,j)*p02(i-1,j)+c2(i,j)*p02(i,j+1)... % he phuong trinh tinh ap suat
                    +d2(i,j)*p02(i,j-1)-f2(i,j))/e2(i,j);
                if p2(i,j)<=0
                    p2(i,j)=0;
                end
            end
        end
    end
    % Kiem tra dieu kien hoi tu
    for i= 1:m+1
        for j= 1:n+1
            S2=S2+abs(p2(i,j)-p02(i,j));
            T2=T2+abs(p2(i,j));
            GAP2=S2/T2;
        end
    end
    p02=p2;
end
% chuong trinh tinh ap suat tinh
while GAP2s>ERR2s
    k2s = k2s+1;
    for i =1:m+1
        for j =1:n+1
            if j==1||j==n+1||i==1||i==m+1
                p2s(i,j)=0;
                %p(i,j)=0;%p(i,j+1);
            elseif i>=46&&i<=116&&j>=19&&j<=73
                p2s(i,j)=Precess2;
            else
                a2s(i,j)= hc2s(i,j)^3;
                b2s(i,j)= ht2s(i,j)^3;
                c2s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h2s(i,j)^3;
                d2s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h2s(i,j)^3;
                e2s(i,j)= a2s(i,j)+b2s(i,j)+c2s(i,j)+d2s(i,j);
                f2s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                p2s(i,j)= (a2s(i,j)*p02s(i+1,j)+b2s(i,j)*p02s(i-1,j)+c2s(i,j)*p02s(i,j+1)...  % he phuong trinh tinh ap suat
                    +d2s(i,j)*p02s(i,j-1)-f2s(i,j))/e2s(i,j);
            end
        end
    end
    for i= 1:m+1
        for j= 1:n+1
            S2s=S2s+abs(p2s(i,j)-p02s(i,j));
            T2s=T2s+abs(p2s(i,j));
            GAP2s=S2s/T2s;
        end
    end
    p02s=p2s;
end
Qd21=0;
Qd22=0;
Qd23=0;
Qs21=0;
Qs22=0;
Qs23=0;
for j=44:48
    %Qd21=Qd21+(1/2)*ld^2*((h2(79,44)-(1/3)*(h2(79,44))^3*((-3*p2(78,j)+4*p2(79,j)-p2(80,j)))/(2*deltaphi))*(deltalanda));
    Qd21=Qd21+(1/2)*ld^2*(h2(79,43)-((1/3)*(h2(79,43)^3)*(3*Precess2-4*p2(78,j)+p2(77,j))/(2*deltaphi)))*deltalanda;
    %Qd22=Qd22+(1/4)*ld^2*(h2(83,44)-(1/3)*(h2(83,44))^3*((3*p2(83,j)-4*p2(84,j)+p2(85,j))/(2*deltaphi))*(deltalanda));
    Qd22=Qd22+(1/2)*ld^2*(-h2(83,43)-((1/3)*(h2(83,43)^3)*(-3*Precess2+4*p2(84,j)-p2(85,j))/(2*deltaphi)))*deltalanda;
end
for j=43
    for i=79:83
        Qd23=Qd23-(1/6)*(h2(i,43))^3*(3*Precess2-4*p2(i,43)+p2(i,42))*(deltaphi/deltalanda);
    end
end
Q2d=Qd21+Qd22+Qd23;
%for i=51
% for j=16:76
%Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(51,16))^3*((-3*p2s(51,j)+4*p2s(50,j)-p2s(49,j))/(2*deltaphi))*(deltalanda));
%Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(111,16))^3*((3*p2s(111,j)-4*p2s(112,j)+p2s(113,j))/(2*deltaphi))*(deltalanda));
%end
%end
for i=46
    for j=19:73
        Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(46,19))^3*((-3*p2s(46,j)+4*p2s(45,j)-p2s(44,j))/(2*deltaphi))*(deltalanda));
        Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(116,19))^3*((3*p2s(116,j)-4*p2s(117,j)+p2s(118,j))/(2*deltaphi))*(deltalanda));
    end
end
%for j=16
%   for i=51:111
%   Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
%Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
% end
%end
for j=19
    for i=46:116
        Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
        %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
    end
end
Q2s=Qs21+Qs22+Qs23;
% can bang momen
Md2=0;
Msr2=0;
Mso12=0;
Mso22=0;
Mso32=0;
Mso42=0;
for i=1:m+1
    for j=1:n+1
        Md2=Md2-R*p2(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=19:73
        Msr2=Msr2-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
for i=1:45
    for j=1:n+1
        Mso12=Mso12-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
for i=117:m+1
    for j=1:n+1
        Mso22=Mso22-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=1:18
        Mso32=Mso32-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=73:n+1
        Mso42=Mso42-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
    end
end
Ms2=Msr2+Mso12+Mso22+Mso32+Mso42;
%% TINH CHO PAD3
%  do day mang dau
for i=1:m+1
    for j=1:n+1
        % phan dong
        phi3(i,j) = phi31+(i-1)*deltaphi;
        phic3(i,j)= phi31+(i-1+1/2)*deltaphi;
        phit3(i,j)= phi31+(i-1-1/2)*deltaphi;
        h3(i,j) = 1+exilon*cos(phi3(i,j)-phibd)-M3*cos(phi3(i,j)-beta3)+capa3*sin(phi31-phi3(i,j));
        hc3(i,j) = 1+exilon*cos(phic3(i,j)-phibd)-M3*cos(phic3(i,j)-beta3)+capa3*sin(phi31-phic3(i,j));
        ht3(i,j) = 1+exilon*cos(phit3(i,j)-phibd)-M3*cos(phit3(i,j)-beta3)+capa3*sin(phi31-phit3(i,j));
        % phan tinh
        h3s(i,j)  = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phi3(i,j));
        hc3s(i,j) = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phic3(i,j));
        ht3s(i,j) = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phit3(i,j));
    end
end
% ch??ng trinh tinh ap suat dong
while GAP3>ERR3
    k3 = k3+1;
    fx3=0;
    fy3=0;
    for i =1:m+1
        for j =1:n+1
            if i==1||i==m+1||j==1||j==n+1
                p3(i,j)=0;
            else
                a3(i,j)= hc3(i,j)^3;
                b3(i,j)= ht3(i,j)^3;
                c3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
                d3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
                e3(i,j)= a3(i,j)+b3(i,j)+c3(i,j)+d3(i,j);
                f3(i,j)= 3*deltaphi*(hc3(i,j)-ht3(i,j));
                p3(i,j)= (a3(i,j)*p03(i+1,j)+b3(i,j)*p03(i-1,j)+c3(i,j)*p03(i,j+1)... % he phuong trinh tinh ap suat
                    +d3(i,j)*p03(i,j-1)-f3(i,j))/e3(i,j);
                if p3(i,j)<=0
                    p3(i,j)=0;
                end
            end
        end
    end
    % Kiem tra dieu kien hoi tu
    for i= 1:m+1
        for j= 1:n+1
            S3=S3+abs(p3(i,j)-p03(i,j));
            T3=T3+abs(p3(i,j));
            GAP3=S3/T3;
        end
    end
    p03=p3;
end
% chuong trinh tinh ap suat tinh
fx1=0;
fy1=0;
for i=1:m+1
    for j=1:n+1
        fx1=fx1-(p1(i,j))*sin(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
        fy1=fy1-(p1(i,j))*cos(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
    end
end
% tren pad2
fx2=0;
fy2=0;
for i=1:m+1
    for j=1:n+1
        fx2=fx2-(p2(i,j))*sin(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
        fy2=fy2-(p2(i,j))*cos(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
    end
end
% trên pad 3
fx3=0;
fy3=0;
for i=1:m+1
    for j=1:n+1
        fx3=fx3-(p3(i,j))*sin(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
        fy3=fy3-(p3(i,j))*cos(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
    end
end
fx=fx1+fx2+fx3;
fy=fy1+fy2+fy3;