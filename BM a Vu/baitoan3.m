clc
clear
L=0.065;% chieu dai o
d=0.0645;% duong kinh truc
R=d/2;
ld=d/L;
m=90;% so phan tu theo chu vi
n=40;% so phan tu theo truc
phi1=0; % gocs ban dau
phi2=2*pi;% goc ket thuc
deltaphi=(phi2-phi1)/m;% buoc goc
deltalanda=2/n;% buoc doc truc
nguy=0.2;
omega=200;
k=1;% he so hoi tu
cm=0.2e-3;
c0=0.5e-3;
xicma=cm/c0;

e=0.16e-3;

P0=nguy*omega*(d/c0)^2/4;
Fx=0;
Fy=0;
ss=1;
v=1;
%while abs(ss)>0.0001
%for phihoa=-0.9

%-------------------------------------------
a1 = 1000;
a2 = 1e-1;
W=108889;
%exilon= 0.197:0.0001:0.199;
exilon= 0.78:0.01:0.81;
phihoa=0.227:0.001:0.230;% goc thai do
EXI=[];
PHIH=[];
%-------------------------------------------





    v=v+1;
% khai baos cac ma tran ban dau ( cac ma tran 41x31)
for y=1:length(exilon)
    for z=1:length(phihoa)
        p0 = zeros(m+1,n+1);
        p  = zeros(m+1,n+1);
        a  = zeros(m+1,n+1);
        b  = zeros(m+1,n+1);
        c  = zeros(m+1,n+1);
        d  = zeros(m+1,n+1);
        e  = zeros(m+1,n+1);
        f  = zeros(m+1,n+1);
        P  = zeros(m+1,n+1);
        % tinh do day mang dau( do day mang dau tinh theo exilon va goc phi)
        %exilon=0.2;%:0.2:0.8 %3 gia tri exilon 
        for i=1:m/2
            for j=1:n+1


            phi(i,j)=(i-1)*deltaphi;%cac goc phi i,j
            phit(i,j)=(i-1-1/2)*deltaphi;% cac goc truoc phi i,j 1/2
            phic(i,j)=(i-1+1/2)*deltaphi;% cac goc sau phi i,j 1/2
            %if 0<=phi(i,j)<pi
               h(i,j)=(((1+xicma)/(2*xicma))+((1-xicma)/(2*xicma))*cos(phi(i,j))-exilon(y)*sin(phihoa(z)-phi(i,j)));
               ht(i,j)=(((1+xicma)/(2*xicma))+((1-xicma)/(2*xicma))*cos(phit(i,j))-exilon(y)*sin(phihoa(z)-phit(i,j)));
               hc(i,j)=(((1+xicma)/(2*xicma))+((1-xicma)/(2*xicma))*cos(phic(i,j))-exilon(y)*sin(phihoa(z)-phic(i,j)));
           % elseif pi<=phi(i,j)<2*pi
                 %h(i,j)=cm*(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phi(i,j))-exilon*sin(phihoa-phi(i,j)));
                  %ht(i,j)=cm*(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phit(i,j))-exilon*sin(phihoa-phit(i,j)));
                   %hc(i,j)=cm*(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phic(i,j))-exilon*sin(phihoa-phic(i,j)));
            %end
                end
        end
        for i=(m/2+1):m+1
            for j=1:n+1
            phi(i,j)=(i-1)*deltaphi;%cac goc phi i,j
            phit(i,j)=(i-1-1/2)*deltaphi;% cac goc truoc phi i,j 1/2
            phic(i,j)=(i-1+1/2)*deltaphi;% cac goc sau phi i,j 1/2
           h(i,j)=(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phi(i,j))-exilon(y)*sin(phihoa(z)-phi(i,j)));
                  ht(i,j)=(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phit(i,j))-exilon(y)*sin(phihoa(z)-phit(i,j)));
                  hc(i,j)=(((1+xicma)/(2*xicma))-((1-xicma)/(2*xicma))*cos(phic(i,j))-exilon(y)*sin(phihoa(z)-phic(i,j)));
            end
        end
        %% chuong trinh tinh ap suat
        % khai bao cac bien ban dau
        S = 0; 
        T = 0;
        %ERR = 10e-3;
        ERR = e-4;
        GAP =1;
        k = 1; % he so lap
        while GAP>ERR
            k = k+1;
            for i =1:m+1
                for j =1:n+1
                    if i==1||i==m+1||j==1||j==n+1
                        p(i,j)=0;
                    else
                        a(i,j)= hc(i,j)^3;
                        b(i,j)= ht(i,j)^3;
                        c(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                        d(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                        e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                        f(i,j)= 3*deltaphi*(hc(i,j)-ht(i,j));
                        p(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)...  % he phuong trinh tinh ap suat
                            +d(i,j)*p0(i,j-1)-f(i,j))/e(i,j);
                        if p(i,j)<=0
                            p(i,j)=0;
                        end
                    end
                 end
            end
            % Kiem tra dieu kien hoi tu
            for i= 1:m+1
                for j= 1:n+1
                    S=S+abs(p(i,j)-p0(i,j));
                    T=T+abs(p(i,j));
                    GAP=S/T;
                end
            end
            p0=p;
        end
        P=p*omega*nguy*(R*R/(c0*c0));
        for i= 1:m+1
           for j=1:n+1
               Fx=Fx+p(i,j)*cos((3*pi/2)-phi(i,j))*deltaphi*deltalanda;% luc theo truc x
               Fy=Fy+p(i,j)*sin((3*pi/2)-phi(i,j))*deltaphi*deltalanda;% luc theo truc y
               FX=Fx*P0;
               FY=Fy*P0;
               ss=FX/Fy;
           end
        end
        %end
        %end
        %x=1:1:361;
        %y=1:1:161;
        %z=p(x,y);
        %figure
        %mesh(z);
        FX;
        FY;
        if (abs(W-FX)<=a1) && (FY/FX< a2)
            EXI = [EXI exilon(y)]
            PHIH = [PHIH phihoa(z)]
        end
    end
end
