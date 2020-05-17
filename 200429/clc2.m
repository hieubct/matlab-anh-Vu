clc
clear
%% tinh cho pad dau tien
m=160;
n=90;
phiTC = 0*pi/3; % goc dau cua pad
phi12 = 2*pi/3;% goc sau cua pad
deltaphi1=(phi12-phiTC)/m;
deltalanda=2/n;
L=150e-3;
RZ=115e-3;% ban kinh pad
RP=100.21e-3;% ban kinh lap rap
CP=0.21e-3;% khe ho lap ráp
anpha=0.09*pi/180;% góc xoay cua pad
capa=RP*anpha/CP;
lrz=2*RZ/L;
%% khai bao cac ma tran he so ban dau
    p0 = zeros(m+1,n+1);
    p  = zeros(m+1,n+1);
    a  = zeros(m+1,n+1);
    b  = zeros(m+1,n+1);
    c  = zeros(m+1,n+1);
    d  = zeros(m+1,n+1);
    e  = zeros(m+1,n+1);
    f  = zeros(m+1,n+1);
    %%xac dinh do day mang dau
     for i=1:m+1
        for j=1:n+1
            phi(i,j) = phiTC+(i-1)*deltaphi1;
            phic(i,j)= phiTC+(i-1+1/2)*deltaphi1;
            phit(i,j)= phiTC+(i-1-1/2)*deltaphi1;
            h(i,j) = -capa*sin(phiTC-phi(i,j));
            hc(i,j) = -capa*sin(phiTC-phic(i,j));
            ht(i,j) =- capa*sin(phiTC-phit(i,j));
            %h(i,j)   = 1+exilon*cosd(phi(i,j)-phibd)-M*cosd(phi(i,j)-beta1)-(xicma(q)/posi)*sind(phi(i,j)-beta1);
            %hc(i,j)   = 1+exilon*cosd(phic(i,j)-phibd)-M*cosd(phic(i,j)-beta1)-(xicma(q)/posi)*sind(phic(i,j)-beta1);
            %ht(i,j)   = 1+exilon*cosd(phit(i,j)-phibd)-M*cosd(phit(i,j)-beta1)-(xicma(q)/posi)*sind(phit(i,j)-beta1);
        end
     end
    S = 0; 
    T = 0;
    ERR = 10e-4;
    GAP =1;
    k = 1; % he so lap
    while GAP>ERR
        k = k+1;
        for i =1:m+1
            for j =1:n+1
                if j==1||j==n+1||i==1||i==m+1
                    p(i,j)=0;
                    %p(i,j)=0;%p(i,j+1);
                elseif i>=40&&i<=60&&j>=20&&j<=30
                   p(i,j)=1;
                else
                    a(i,j)= hc(i,j)^3;
                   b(i,j)= ht(i,j)^3;
                    c(i,j)= lrz^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    d(i,j)= lrz^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                   f(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                  p(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)...  % he phuong trinh tinh ap suat
                    +d(i,j)*p0(i,j-1)-f(i,j))/e(i,j);
                       %if p(i,j)<=0
                            %p(i,j)=0;
                       %end
                end
           end
        end
    for i= 1:m+1
        for j= 1:n+1
            S=S+abs(p(i,j)-p0(i,j));
            T=T+abs(p(i,j));
            GAP=S/T;
        end
    end
     p0=p;
    end