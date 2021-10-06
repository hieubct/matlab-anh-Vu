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


%%
% trong ch??ng tr�nh n�y, m�nh m?i ch? t�m ???c Precess2 theo ?i?u ki?n c�n b?ng l?u l??ng (Qd2=Qs2), sau khi thu ???c Precess2, m�nh s? ti?n h�nh ki?m tra ?i?u ki?n c�n b?ng m� men, n?u Md2#Ms2 th� l?i thay ??i anpha2 v� ch?y l?i ch??ng tr�nh, c? nh? v?y ?? t�m ra anpha2 th?a m�n ?k Md2=Ms2 (L�c n�y Precess2 c?ng s? ??t gi� tr? m?i).
% B�y gi? a mu?n t�m anpha2 t? ??ng (th?a m�n ?k c�n b?ng momen) b?ng c�ch s? d?ng thu?t to�n t�m nghi?m gi?ng nh? t�m Precess2(theo ?k c�n b?ng l?u l??ng).

%%

exilonn1=0;%1.5e-5;% do nang cua pad1
exilonn2=0;%7e-6;% do nang cua pad1
exilonn3=0;%3e-6;% do nang cua pad1
M =0.57;%
M1 = 0.57;%1-(cbtt1/cl);
M2 = 0.57;%1-(cbtt2/cl);
M3 = 0.57;%1-(cbtt3/cl);
anpha11 = 0;%0.035*pi/180;%0.009*pi/180;% goc nghieng cua pad1
anpha12 = 0.1*pi/180;%0.009*pi/180;% goc nghieng cua pad1
anpha1 = (anpha11+anpha12)/2;  %0.0167*pi/180;%0.0165*pi/180;%0.0168*pi/180;%0.017*pi/180;%0.016*pi/180;%0.015*pi/180;%0.018*pi/180;%goc nghieng cua pad2

anpha21 = 0.002*pi/180;
anpha22 = 0.05*pi/180;
anpha2 = (anpha21+anpha22)/2;  %0.0167*pi/180;%0.0165*pi/180;%0.0168*pi/180;%0.017*pi/180;%0.016*pi/180;%0.015*pi/180;%0.018*pi/180;%goc nghieng cua pad2
%anpha2 =0.0117*pi/180;%0.0167*pi/180;%0.0165*pi/180;%0.0168*pi/180;%0.017*pi/180;%0.016*pi/180;%0.015*pi/180;%0.018*pi/180;%goc nghieng cua pad2




% hs_anpha2 = (anpha2*180)/pi;
% hs_anpha21 = (anpha21*180)/pi;
% hs_anpha22 = (anpha22*180)/pi;



anpha31 = 0*pi/180;
anpha32 = 0.02*pi/180;
anpha3 = (anpha31+anpha32)/2;

exilon = 0.2;% ti le lech tam
phibd =0.985;%0.9;%0.95;%1;%1.23;% goc lech ban dau cua truc


Precess11=0;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
Precess12=2.5;%0.9;%0.7;%1;%
Precess21=1;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
Precess22=2;%0.9;%0.7;%1;%
Precess31=0;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
Precess32=3;%0.9;%0.7;%1;%


ld = l/d1; % ti le chieu dai tren duong kinh
m = 204; % so luoi theo chu vi
n = 150; % so luoi theo phuong doc truc
% thong so pad1
phi11 = 3*pi/180;%3*pi/180; % goc dau cua pad
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

err = 10e-3;


Precess1 = (Precess11+Precess12)/2;
Precess2 = (Precess21+Precess22)/2;
Precess3 = (Precess31+Precess32)/2;
kk1=0;
cc1=0;
kk2=0;
cc2=0;
kk3=0;
cc3=0;
k=0;
c=0;
while true
    
    capa1=Rt*anpha1/cl;
    % ma tran he so ban dau cho phan hydrodynamic
    p01 = zeros(m+1,n+1);
    p1 = zeros(m+1,n+1);
    a1 = zeros(m+1,n+1);
    b1 = zeros(m+1,n+1);
    c1 = zeros(m+1,n+1);
    d1 = zeros(m+1,n+1);
    e1 = zeros(m+1,n+1);
    f1 = zeros(m+1,n+1);
    % ma tran he so cho phan hydrostatic
    p01s = zeros(m+1,n+1);
    p1s  = zeros(m+1,n+1);
    a1s  = zeros(m+1,n+1);
    b1s  = zeros(m+1,n+1);
    c1s  = zeros(m+1,n+1);
    d1s  = zeros(m+1,n+1);
    e1s  = zeros(m+1,n+1);
    f1s  = zeros(m+1,n+1);
    
    % khai bao cac bien ban dau
    S1 = 0;
    T1 = 0;
    ERR1 = 10e-4;
    GAP1 =1;
    k1 = 1; % he so lap
    
    S1s = 0;
    T1s = 0;
    ERR1s = 10e-4;
    GAP1s =1;
    k1s = 1; % he so lap
    
    
    
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
                elseif i>=72&&i<=74&&j>=75&&j<=77
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
                elseif i>=62&&i<=144&&j>=31&&j<=121
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
    for j=75:77%102:104   73:79%100:106
        %Qd21=Qd21+(1/2)*ld^2*(h2(100,72)-((1/3)*(h2(100,72)^3)*(3*Precess2-4*p2(99,j)+p2(98,j))/(2*deltaphi)))*deltalanda;
        Qd11=Qd11+(1/2)*ld^2*(h1(72,74)-((1/3)*(h1(72,74)^3)*(3*Precess1-4*p1(71,j)+p1(70,j))/(2*deltaphi)))*deltalanda;
        %Qd22=Qd22+(1/2)*ld^2*(-h2(106,72)-((1/3)*(h2(106,72)^3)*(-3*Precess2+4*p2(107,j)-p2(108,j))/(2*deltaphi)))*deltalanda;
        Qd12=Qd12+(1/2)*ld^2*(-h1(74,74)-((1/3)*(h1(74,74)^3)*(-3*Precess1+4*p1(77,j)-p1(76,j))/(2*deltaphi)))*deltalanda;
    end
    for j=74%j=72
        for i=72:74%
            %Qd23=Qd23-(1/6)*(h2(i,44))^3*(3*Precess2-4*p2(i,44)+p2(i,43))*(deltaphi/deltalanda);
            Qd13=Qd13-(1/6)*(h1(i,74))^3*(3*Precess1-4*p1(i,74)+p1(i,73))*(deltaphi/deltalanda);
        end
    end
    Q1d=Qd11+Qd12+Qd13;
    %for i=51
    % for j=16:76
    %Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(51,16))^3*((-3*p2s(51,j)+4*p2s(50,j)-p2s(49,j))/(2*deltaphi))*(deltalanda));
    %Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(111,16))^3*((3*p2s(111,j)-4*p2s(112,j)+p2s(113,j))/(2*deltaphi))*(deltalanda));
    %end
    %end
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62
        for j=31:121
            Qs11=Qs11-(1/2)*ld^2*((1/3)*(h1s(62,31))^3*((-3*p1s(62,j)+4*p1s(61,j)-p1s(60,j))/(2*deltaphi))*(deltalanda));
            Qs12=Qs12-(1/2)*ld^2*(-(1/3)*(h1s(144,31))^3*((3*p1s(144,j)-4*p1s(145,j)+p1s(146,j))/(2*deltaphi))*(deltalanda));
        end
    end
    %for j=16
    %   for i=51:111
    %   Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
    %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
    % end
    %end
    for j=31
        for i=62:144
            Qs13=Qs13+(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j-1)+p1s(i,j-2))*(deltaphi/deltalanda);
            %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
        end
    end
    Q1s=Qs11+Qs12+Qs13;
    
    
    
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
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62:144
        for j=31:121
            Msr1=Msr1-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
        end
    end
    for i=1:61
        for j=1:n+1
            Mso11=Mso11-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
        end
    end
    for i=145:m+1
        for j=1:n+1
            Mso21=Mso21-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=1:30
            Mso31=Mso31-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=122:n+1
            Mso41=Mso41-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
        end
    end
    Ms1=Msr1+Mso11+Mso21+Mso31+Mso41;
    
    
    
    kk1=kk1+1;
    condition11 = abs((Q1d-Q1s)/Q1d);
    condition12 = abs((Q1d-Q1s)/Q1s);
    if(condition11 > err && condition12 > err && kk1<10)
        %Kiểm tra nếu
        condition13 = Q1d-Q1s;
        if(condition13 > 0)
            Precess11 = Precess1;
        else
            Precess12 = Precess1;
        end
        Precess1 = (Precess11+Precess12)/2;
        continue
    end
    
    conditionM11 = abs((Md1-Ms1)/Md1);
    conditionM12 = abs((Md1-Ms1)/Ms1);
    if(conditionM11 > err && conditionM12 > err)
        conditionM13 = Md1-Ms1;
        if(conditionM13 < 0)
            anpha11 = anpha1;
        else
            anpha12 = anpha1;
        end
        anpha1 = (anpha11+anpha12)/2;
        
        Precess11=0;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
        Precess12=2.5;%0.9;%0.7;%1;%
        Precess1 = (Precess11+Precess12)/2;
        kk1=0;
        
        cc1=cc1+1;
        
        
        Md1;
        Ms1;
        hs_anpha1 = (anpha1*180)/pi;
        hs_anpha11 = (anpha11*180)/pi;
        hs_anpha12 = (anpha12*180)/pi;
        
        disp('---------------');
        
        
        
        if(cc1==10)
            break
        end
        
        continue
    end
    
end


%% TINH CHO PAD2
%  do day mang dau
while true
    
    %-----------------------------------
    
    capa2=Rt*anpha2/cl;
    
    p02 = zeros(m+1,n+1);
    p2 = zeros(m+1,n+1);
    a2 = zeros(m+1,n+1);
    b2 = zeros(m+1,n+1);
    c2 = zeros(m+1,n+1);
    d2 = zeros(m+1,n+1);
    e2 = zeros(m+1,n+1);
    f2 = zeros(m+1,n+1);
    
    p02s = zeros(m+1,n+1);
    p2s  = zeros(m+1,n+1);
    a2s  = zeros(m+1,n+1);
    b2s  = zeros(m+1,n+1);
    c2s  = zeros(m+1,n+1);
    d2s  = zeros(m+1,n+1);
    e2s  = zeros(m+1,n+1);
    f2s  = zeros(m+1,n+1);
    
    S2 = 0;
    T2 = 0;
    ERR2 = 10e-4;
    GAP2 =1;
    k2 = 1; % he so lap
    
    S2s = 0;
    T2s = 0;
    ERR2s = 10e-4;
    GAP2s =1;
    k2s = 1; % he so lap
    
    
    
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
                elseif i>=102&&i<=104&&j>=75&&j<=77%i>=80&&i<=82&&j>=44&&j<=48
                    p2(i,j)=Precess2;
                    %end
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
                elseif i>=62&&i<=144&&j>=31&&j<=121
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
    for j=75:77%102:104   73:79%100:106
        %Qd21=Qd21+(1/2)*ld^2*(h2(100,72)-((1/3)*(h2(100,72)^3)*(3*Precess2-4*p2(99,j)+p2(98,j))/(2*deltaphi)))*deltalanda;
        Qd21=Qd21+(1/2)*ld^2*(h2(102,74)-((1/3)*(h2(102,74)^3)*(3*Precess2-4*p2(101,j)+p2(100,j))/(2*deltaphi)))*deltalanda;
        %Qd22=Qd22+(1/2)*ld^2*(-h2(106,72)-((1/3)*(h2(106,72)^3)*(-3*Precess2+4*p2(107,j)-p2(108,j))/(2*deltaphi)))*deltalanda;
        Qd22=Qd22+(1/2)*ld^2*(-h2(104,74)-((1/3)*(h2(104,74)^3)*(-3*Precess2+4*p2(107,j)-p2(106,j))/(2*deltaphi)))*deltalanda;
    end
    for j=74%j=72
        for i=102:104%
            %Qd23=Qd23-(1/6)*(h2(i,44))^3*(3*Precess2-4*p2(i,44)+p2(i,43))*(deltaphi/deltalanda);
            Qd23=Qd23-(1/6)*(h2(i,74))^3*(3*Precess2-4*p2(i,74)+p2(i,73))*(deltaphi/deltalanda);
        end
    end
    Q2d=Qd21+Qd22+Qd23;
    %for i=51
    % for j=16:76
    %Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(51,16))^3*((-3*p2s(51,j)+4*p2s(50,j)-p2s(49,j))/(2*deltaphi))*(deltalanda));
    %Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(111,16))^3*((3*p2s(111,j)-4*p2s(112,j)+p2s(113,j))/(2*deltaphi))*(deltalanda));
    %end
    %end
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62
        for j=31:121
            Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(62,31))^3*((-3*p2s(62,j)+4*p2s(61,j)-p2s(60,j))/(2*deltaphi))*(deltalanda));
            Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(144,31))^3*((3*p2s(144,j)-4*p2s(145,j)+p2s(146,j))/(2*deltaphi))*(deltalanda));
        end
    end
    %for j=16
    %   for i=51:111
    %   Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
    %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
    % end
    %end
    for j=31
        for i=62:144
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
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62:144
        for j=31:121
            Msr2=Msr2-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
        end
    end
    for i=1:61
        for j=1:n+1
            Mso12=Mso12-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
        end
    end
    for i=145:m+1
        for j=1:n+1
            Mso22=Mso22-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=1:30
            Mso32=Mso32-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=122:n+1
            Mso42=Mso42-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
        end
    end
    Ms2=Msr2+Mso12+Mso22+Mso32+Mso42;
    
    
   
    kk2=kk2+1;
    condition21 = abs((Q2d-Q2s)/Q2d);
    condition22 = abs((Q2d-Q2s)/Q2s);
    if(condition21 > err && condition22 > err && kk2<10)
        %Kiểm tra nếu
        condition23 = Q2d-Q2s;
        if(condition23 > 0)
            Precess21 = Precess2;
        else
            Precess22 = Precess2;
        end
        Precess2 = (Precess21+Precess22)/2;
        continue
    end
    
    conditionM21 = abs((Md2-Ms2)/Md2);
    conditionM22 = abs((Md2-Ms2)/Ms2);
    if(conditionM21 > err && conditionM22 > err)
        conditionM23 = Md2-Ms2;
        if(conditionM23 < 0)
            anpha21 = anpha2;
        else
            anpha22 = anpha2;
        end
        anpha2 = (anpha21+anpha22)/2;
        
        Precess21=1;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
        Precess22=2;%0.9;%0.7;%1;%
        Precess2 = (Precess21+Precess22)/2;
        kk2=0;
        
        cc2=cc2+1;
        
        
        Md2;
        Ms2;
        hs_anpha2 = (anpha2*180)/pi;
        hs_anpha21 = (anpha21*180)/pi;
        hs_anpha22 = (anpha22*180)/pi;
        
        disp('---------------');

        if(cc2==10)
            break
        end
        
        continue
    end
end








%% TINH CHO PAD3
%  do day mang dau
while true
    %-----------------------------------
    capa3=Rt*anpha3/cl;
    
    p03 = zeros(m+1,n+1);
    p3 = zeros(m+1,n+1);
    a3 = zeros(m+1,n+1);
    b3 = zeros(m+1,n+1);
    c3 = zeros(m+1,n+1);
    d3 = zeros(m+1,n+1);
    e3 = zeros(m+1,n+1);
    f3 = zeros(m+1,n+1);
    
    p03s = zeros(m+1,n+1);
    p3s  = zeros(m+1,n+1);
    a3s  = zeros(m+1,n+1);
    b3s  = zeros(m+1,n+1);
    c3s  = zeros(m+1,n+1);
    d3s  = zeros(m+1,n+1);
    e3s  = zeros(m+1,n+1);
    f3s  = zeros(m+1,n+1);
    
    S3 = 0;
    T3= 0;
    ERR3 = 10e-4;
    GAP3 =1;
    k3 = 1; % he so lap
    
    S3s = 0;
    T3s = 0;
    ERR3s = 10e-4;
    GAP3s =1;
    k3s = 1; % he so lap
    
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
                elseif i>=50&&i<=52&&j>=75&&j<=77%i>=80&&i<=82&&j>=44&&j<=48
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
    
    % chuong trinh tinh ap suat tinh
    while GAP3s>ERR3s
        k3s = k3s+1;
        for i =1:m+1
            for j =1:n+1
                if j==1||j==n+1||i==1||i==m+1
                    p3s(i,j)=0;
                elseif i>=62&&i<=144&&j>=31&&j<=121
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
    for j=75:77%102:104   73:79%100:106
        %Qd21=Qd21+(1/2)*ld^2*(h2(100,72)-((1/3)*(h2(100,72)^3)*(3*Precess2-4*p2(99,j)+p2(98,j))/(2*deltaphi)))*deltalanda;
        Qd31=Qd31+(1/2)*ld^2*(h3(50,74)-((1/3)*(h3(50,74)^3)*(3*Precess3-4*p3(49,j)+p2(48,j))/(2*deltaphi)))*deltalanda;
        %Qd22=Qd22+(1/2)*ld^2*(-h2(106,72)-((1/3)*(h2(106,72)^3)*(-3*Precess2+4*p2(107,j)-p2(108,j))/(2*deltaphi)))*deltalanda;
        Qd32=Qd32+(1/2)*ld^2*(-h3(52,74)-((1/3)*(h3(52,74)^3)*(-3*Precess3+4*p3(55,j)-p2(54,j))/(2*deltaphi)))*deltalanda;
    end
    for j=74%j=72
        for i=50:52%
            %Qd23=Qd23-(1/6)*(h2(i,44))^3*(3*Precess2-4*p2(i,44)+p2(i,43))*(deltaphi/deltalanda);
            Qd33=Qd33-(1/6)*(h3(i,74))^3*(3*Precess3-4*p3(i,74)+p3(i,73))*(deltaphi/deltalanda);
        end
    end
    Q3d=Qd31+Qd32+Qd33;
    %for i=51
    % for j=16:76
    %Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(51,16))^3*((-3*p2s(51,j)+4*p2s(50,j)-p2s(49,j))/(2*deltaphi))*(deltalanda));
    %Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(111,16))^3*((3*p2s(111,j)-4*p2s(112,j)+p2s(113,j))/(2*deltaphi))*(deltalanda));
    %end
    %end
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62
        for j=31:121
            Qs31=Qs31-(1/2)*ld^2*((1/3)*(h3s(62,31))^3*((-3*p3s(62,j)+4*p3s(61,j)-p3s(60,j))/(2*deltaphi))*(deltalanda));
            Qs32=Qs32-(1/2)*ld^2*(-(1/3)*(h3s(144,31))^3*((3*p3s(144,j)-4*p3s(145,j)+p3s(146,j))/(2*deltaphi))*(deltalanda));
        end
    end
    %for j=16
    %   for i=51:111
    %   Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
    %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
    % end
    %end
    for j=31
        for i=62:144
            Qs33=Qs33+(1/6)*((h3s(i,j))^3)*(3*p3s(i,j)-4*p3s(i,j-1)+p3s(i,j-2))*(deltaphi/deltalanda);
            %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
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
    %i>=62&&i<=144&&j>=31&&j<=121
    for i=62:144
        for j=31:121
            Msr3=Msr3-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
        end
    end
    for i=1:61
        for j=1:n+1
            Mso13=Mso13-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
        end
    end
    for i=145:m+1
        for j=1:n+1
            Mso23=Mso23-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=1:30
            Mso33=Mso33-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
        end
    end
    for i=62:144
        for j=122:n+1
            Mso43=Mso43-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
        end
    end
    Ms3=Msr3+Mso13+Mso23+Mso33+Mso43;
    
    
    
    %     %%
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
    
    
    
    kk3=kk3+1;
    condition31 = abs((Q3d-Q3s)/Q3d);
    condition32 = abs((Q3d-Q3s)/Q3s);
    if(condition31 > err && condition32 > err && kk3<10)
        %Kiểm tra nếu
        condition33 = Q3d-Q3s;
        if(condition33 > 0)
            Precess31 = Precess3;
        else
            Precess32 = Precess3;
        end
        Precess3 = (Precess31+Precess32)/2;
        continue
    end
    
    conditionM31 = abs((Md3-Ms3)/Md3);
    conditionM32 = abs((Md3-Ms3)/Ms3);
    if(conditionM31 > err && conditionM32 > err)
        conditionM33 = Md3-Ms3;
        if(conditionM33 < 0)
            anpha31 = anpha3;
        else
            anpha32 = anpha3;
        end
        anpha3 = (anpha31+anpha32)/2;
        
        Precess31=0;%0.87;%0.86;%0.85;%0.86;%0.84;%0.82;%0.84;%0.9;%
        Precess32=2.5;%0.9;%0.7;%1;%
        Precess3 = (Precess31+Precess32)/2;
        kk3=0;
        
        cc3=cc3+1;
        
        
        Md3;
        Ms3;
        hs_anpha3 = (anpha3*180)/pi;
        hs_anpha31 = (anpha31*180)/pi;
        hs_anpha32 = (anpha32*180)/pi;
        
        disp('---------------');
        
        
        
        if(cc3==10)
            break
        end
        
        continue
    end
    
end

%     if(condition11 < err && condition12 < err && conditionM11 < err && conditionM12 < err)
%         break
%     end
    % if(condition1 < err && condition2 < err && conditionM21 < err && conditionM22 < err)
    %     break
    % end





Md1
Ms1
hs_anpha1 = anpha1*180/pi
Q1d
Q1s
Precess1
kk1
cc1


Md2
Ms2
hs_anpha2 = anpha2*180/pi
Q2d
Q2s
Precess2
kk2
cc2

Md3
Ms3
hs_anpha3 = anpha3*180/pi
Q3d
Q3s
Precess3
kk3
cc3




% p1(103,76);
% p2(103,76);
% p3(103,76);
% p1(72,76);
% p3(50,76);