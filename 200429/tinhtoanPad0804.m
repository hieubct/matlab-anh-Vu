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
cbtt1=cb-exilonn1;
exilonn2=0;%7e-6;% do nang cua pad1
cbtt2=cb-exilonn2;
exilonn3=0;%3e-6;% do nang cua pad1
cbtt3=cb-exilonn3;
%M =0.57;%
M1 = 0.5714;%1-(cbtt1/cl);
M2 = 1-(cbtt2/cl);
M3 = 1-(cbtt3/cl);
anpha1 = 0.0361*pi/180;%0.009*pi/180;% goc nghieng cua pad1
anpha2 = 0.0159*pi/180;% goc nghieng cua pad2
anpha3 = 0.0434*pi/180;% goc nghieng cua pad3
Precess1=0.87;% áp su?t vùng lõm2
Precess2=8.4;%3.43;% áp su?t vùng lõm2
Precess3=1.31;%3.43;% áp su?t vùng lõm2
capa1=Rt*anpha1/cl;1
capa2=Rt*anpha2/cl;
capa3=Rt*anpha3/cl;
ld = l/d1; % ti le chieu dai tren duong kinh
m = 160; % so luoi theo chu vi
n = 90; % so luoi theo phuong doc truc
%e1 =0.75e-04;% 0.1e-3; % do lech tam ban dau
exilon = 0.25;%0.05;%e1/c1; % ti le lech tam
phibd = 0.35; % goc lech ban dau cua truc
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
        h1(i,j)   = 1+exilon*cos(phi1(i,j)-phibd)-M1*cos(phi1(i,j)-beta1)+capa1*sin(phi11-phi1(i,j));
        hc1(i,j)  = 1+exilon*cos(phic1(i,j)-phibd)-M1*cos(phic1(i,j)-beta1)+capa1*sin(phi11-phic1(i,j));
        ht1(i,j)  = 1+exilon*cos(phit1(i,j)-phibd)-M1*cos(phit1(i,j)-beta1)+capa1*sin(phi11-phit1(i,j));
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
            elseif i>=56&&i<=58&&j>=45&&j<=47
                p1(i,j)=Precess1;
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
            elseif i>=46&&i<=116&&j>=19&&j<=73
                p1s(i,j)=Precess1;
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
Qd11=0;
Qd12=0;
Qd13=0;
Qs11=0;
Qs12=0;
Qs13=0;
for j=44:48
    Qd11=Qd11+(1/2)*ld^2*(h1(55,43)-((1/3)*(h1(55,43)^3)*(3*Precess1-4*p1(54,j)+p1(53,j))/(2*deltaphi)))*deltalanda;
    Qd12=Qd12+(1/2)*ld^2*(-h1(59,43)-((1/3)*(h1(59,43)^3)*(-3*Precess1+4*p1(60,j)-p1(61,j))/(2*deltaphi)))*deltalanda;
end
for j=43
    for i=55:59
        Qd13=Qd13-(1/6)*(h1(i,43))^3*(3*Precess1-4*p1(i,43)+p1(i,42))*(deltaphi/deltalanda);
    end
end
Q1d=Qd11+Qd12+Qd13;


for j=19:73
    Qs11=Qs11-(1/2)*ld^2*((1/3)*(h1s(46,19))^3*((-3*p1s(46,j)+4*p1s(45,j)-p1s(44,j))/(2*deltaphi))*(deltalanda));
    Qs12=Qs12-(1/2)*ld^2*(-(1/3)*(h1s(116,19))^3*((3*p1s(116,j)-4*p1s(117,j)+p1s(118,j))/(2*deltaphi))*(deltalanda));
end
for j=19
    for i=46:116
        Qs13=Qs13+(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j-1)+p1s(i,j-2))*(deltaphi/deltalanda);
    end
end
Q1s=Qs11+Qs12+Qs13;
% can bang momen
Md1=0;
Msr1=0;
Mso11=0;
Mso21=0;
Mso31=0;
Mso41=0;
for i=1:m+1
    for j=1:n+1
        Md1=Md1-R*p1(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=19:73
        Msr1=Msr1-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
for i=1:45
    for j=1:n+1
        Mso11=Mso11-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
for i=117:m+1
    for j=1:n+1
        Mso21=Mso21-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=1:18
        Mso31=Mso31-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=73:n+1
        Mso41=Mso41-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
    end
end
Ms1=Msr1+Mso11+Mso21+Mso31+Mso41;
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
            elseif i>=80&&i<=82&&j>=45&&j<=47
                p2(i,j)=Precess2;
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
    Qd21=Qd21+(1/2)*ld^2*(h2(79,43)-((1/3)*(h2(79,43)^3)*(3*Precess2-4*p2(78,j)+p2(77,j))/(2*deltaphi)))*deltalanda;
    Qd22=Qd22+(1/2)*ld^2*(-h2(83,43)-((1/3)*(h2(83,43)^3)*(-3*Precess2+4*p2(84,j)-p2(85,j))/(2*deltaphi)))*deltalanda;
end
for j=43
    for i=79:83
        Qd23=Qd23-(1/6)*(h2(i,43))^3*(3*Precess2-4*p2(i,43)+p2(i,42))*(deltaphi/deltalanda);
    end
end
Q2d=Qd21+Qd22+Qd23;
for i=46
    for j=19:73
        Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(46,19))^3*((-3*p2s(46,j)+4*p2s(45,j)-p2s(44,j))/(2*deltaphi))*(deltalanda));
        Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(116,19))^3*((3*p2s(116,j)-4*p2s(117,j)+p2s(118,j))/(2*deltaphi))*(deltalanda));
    end
end
for j=19
    for i=46:116
        Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
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
            elseif i>=38&&i<=40&&j>=45&&j<=47
                p3(i,j)=Precess3;
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
while GAP3s>ERR3s
    k3s = k3s+1;
    for i =1:m+1
        for j =1:n+1
            if j==1||j==n+1||i==1||i==m+1
                p3s(i,j)=0;
            elseif i>=46&&i<=116&&j>=19&&j<=73
                p3s(i,j)=Precess3;
            else
                a3s(i,j)= hc3s(i,j)^3;
                b3s(i,j)= ht3s(i,j)^3;
                c3s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h3s(i,j)^3;
                d3s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h3s(i,j)^3;
                e3s(i,j)= a3s(i,j)+b3s(i,j)+c3s(i,j)+d3s(i,j);
                f3s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                p3s(i,j)= (a3s(i,j)*p03s(i+1,j)+b3s(i,j)*p03s(i-1,j)+c3s(i,j)*p03s(i,j+1)...  % he phuong trinh tinh ap suat
                    +d3s(i,j)*p03s(i,j-1)-f3s(i,j))/e3s(i,j);
                %if p(i,j)<=0
                %p(i,j)=0;
                %end
            end
        end
    end
    for i= 1:m+1
        for j= 1:n+1
            S3s=S3s+abs(p3s(i,j)-p03s(i,j));
            T3s=T3s+abs(p3s(i,j));
            GAP3s=S3s/T3s;
        end
    end
    p03s=p3s;
end
Qd31=0;
Qd32=0;
Qd33=0;
Qs31=0;
Qs32=0;
Qs33=0;
for j=44:48
    Qd31=Qd31+(1/2)*ld^2*(h3(37,43)-((1/3)*(h3(37,43)^3)*(3*Precess3-4*p3(36,j)+p3(35,j))/(2*deltaphi)))*deltalanda;
    Qd32=Qd32+(1/2)*ld^2*(-h3(41,43)-((1/3)*(h3(41,43)^3)*(-3*Precess3+4*p3(42,j)-p3(43,j))/(2*deltaphi)))*deltalanda;
end
for j=43
    for i=37:41
        Qd33=Qd33-(1/6)*(h3(i,43))^3*(3*Precess3-4*p3(i,43)+p3(i,42))*(deltaphi/deltalanda);
    end
end
Q3d=Qd31+Qd32+Qd33;


for j=19:73
    Qs31=Qs31-(1/2)*ld^2*((1/3)*(h3s(46,19))^3*((-3*p3s(46,j)+4*p3s(45,j)-p3s(44,j))/(2*deltaphi))*(deltalanda));
    Qs32=Qs32-(1/2)*ld^2*(-(1/3)*(h3s(116,19))^3*((3*p3s(116,j)-4*p3s(117,j)+p3s(118,j))/(2*deltaphi))*(deltalanda));
end
for j=19
    for i=46:116
        Qs33=Qs33+(1/6)*((h3s(i,j))^3)*(3*p3s(i,j)-4*p3s(i,j-1)+p3s(i,j-2))*(deltaphi/deltalanda);
    end
end
Q3s=Qs31+Qs32+Qs33;
% can bang momen
Md3=0;
Msr3=0;
Mso13=0;
Mso23=0;
Mso33=0;
Mso43=0;
for i=1:m+1
    for j=1:n+1
        Md3=Md3-R*p3(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=19:73
        Msr3=Msr3-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
for i=1:45
    for j=1:n+1
        Mso13=Mso13-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
for i=117:m+1
    for j=1:n+1
        Mso23=Mso23-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=1:18
        Mso33=Mso33-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
for i=46:116
    for j=73:n+1
        Mso43=Mso43-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
    end
end
Ms3=Msr3+Mso13+Mso23+Mso33+Mso43;
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



% Ve do thi
for i=1:m+1
    phi11(i,1) = phi1(i,1)*180/pi;
    phi22(i,1) = phi2(i,1)*180/pi;
    phi33(i,1) = phi3(i,1)*180/pi;
end

figure
[X,Y]=meshgrid(1:1:91,phi11(:,1));
Z=p1;
mesh(X,Y,Z)
hold on
[X2,Y2]=meshgrid(1:1:91,phi22(:,1));
Z2=p2;
mesh(X2,Y2,Z2)
mesh(X,Y,Z)
[X3,Y3]=meshgrid(1:1:91,phi33(:,1));
Z3=p3;
mesh(X3,Y3,Z3)
hold off


figure
[Xs,Ys]=meshgrid(1:1:91,phi11(:,1));
Zs=p1s;
mesh(Xs,Ys,Zs)
hold on
[X2s,Y2s]=meshgrid(1:1:91,phi22(:,1));
Z2s=p2s;
mesh(X2s,Y2s,Z2s)
[X3s,Y3s]=meshgrid(1:1:91,phi33(:,1));
Z3s=p3s;
mesh(X3s,Y3s,Z3s)
hold off
colorbar

figure
plot(max(transpose(p1)),transpose(phi11))
hold on
plot(max(transpose(p2)),transpose(phi22))
plot(max(transpose(p3)),transpose(phi33))





