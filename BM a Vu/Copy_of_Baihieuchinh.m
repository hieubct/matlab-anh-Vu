%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHUONG TRINH GIAI BAI TOAN THUAT PHONG %
%               GIAI DOAN                %
%-------------------------------------------------
function dao_dong_day_4
clear all;
close all;
clf;
warning on all
set(gcf,'Visible','off')              % turns current figure "off"
%set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
%% Cac thong so dau vao
global n  m md g k mu dl C1 K1 pk Fth dt B k1 ti pt teta0 P mk pa Jz Ltk Lta Vb nu R Tdo Cp vk b K0k fi2
global deltaY Ltkp alphaP ckd Ltk0
L1    = 70.0;                           % m - Chieu dai day
n     = 70;                             % So phan tu duoc chia
teta0 = 70*pi/180;                      % Radian - Goc phong ban dau
dl    = L1/n;                           % m - Do dai phan tu
g     = 9.81;                           % m/s2 - Gia toc trong truong
C1    = 1e4;                           % N - Do cung cua day
K1    = 3;                              % Ns - He so can nhot cua day
mu    = 0.05;                          % kg/m - Mat do tren met dai cua day
v0    = 10;                             % m/s - Van toc ban dau cua dau keo 10.8277
md    = 1.2;                            % kg - khoi luong dau keo
pa    = 101325.0;                       % Ap suat khi quyen
Vb    = 0.00073;                        % The tich cua binh khi
nu    = 29e-3;                          % kg/mol cua khong khi
R     = 8.314;                          % [m3�Pa�K-1�mol-1]-Hang so khi ly tuong  
Tdo   = 298;

dt    = 0.001;                          % s - Buoc tinh
ti    = 0;
k1    = 1.4;                            % He so doan nhiet
fi2   = 0.98;                           % he so ton that loa phut

pk    = 25e6  ;                         % N/m^2  %310 bar = 3.1e7 Pa = 3.1e7 N/m^2
mk    = nu*pk*Vb/(R*Tdo);               % 338.84e-3; %Kg - khoi luong khi trong binh
m     = md+mk;                          % kg - Khoi luong toan bo dau dan va khi nen
vk    = Vb/mk;                          % m^3/kg 
Fth   = pi*4e-6;                        % m^2 - dien tich tiet dien toi han
K0k   = sqrt(2.0*k1/(k1+1))*(2.0/(k1+1))^(1.0/(k1-1));
b     = fi2*K0k*k1*Fth/mk*sqrt(pk/vk);
B     = (k1-1)*b/(2*k1);
pt    = pk;
Cp    = k1*(2/(k1+1))^(k1/(k1-1))*sqrt((k1+1)*(1-(pa/pk)^((k1-1)/k1))/(k1-1));
P     = [];

Jz   = 77129e-6; %kg.m2 ~ 0.8*0.34^2/12
Lta  = 0.;
Ltk  = 0.11;
Ltk0 = Ltk;
deltaY = 0.015;
Ltkp   = sqrt(Ltk^2+deltaY^2);
alphaP = atan(deltaY/Ltk);
ckd    = 0.0001;
h0     = -.0;
%% MO PHONG

figure(1)
axis equal
xlim([-10 80])
ylim([-10 40])
h=line(0,0,'color','b','linew',2);
hold on
title('KS chuyen dong cua day dan hoi','fontangle','italic','fontweight','bold','color','r')
ylabel('y(m)','color','b')
xlabel('x(m)','color','b')
grid on
grid minor
h1=plot(0,0,'r.','markersize',5);
h2=plot(0,0,'k.','markersize',15);
h3=line(0,0,'color',[0.3 0.3 0.3],'linew',4);
set(gcf,'renderer','zbuffer')
%% 

t0 = 0;
% state vector 
y0 = [0 v0*cos(teta0) 0 v0*sin(teta0) 0 teta0];
k  = 0;           % number of strings elements + 1
j  = 0;
L  = 0;
t(1)=0;
TTx0(1)= 0;   
TTy0(1)= 0;
XX(1)=0;        
TTxm(1)= 0;   
TTym(1)= 0;
Pda(1)=0;
Gieta(1)=0;
Mza(1) = 0;
%%

for i=1:n-1
    k=k+1;
    L=0;
    while L<=dl                         % keo tung phan tu cua chuoi len
        j=j+1;
        [y,Tx,Ty,Pd,Mz]=runKut4(@odef,t0,y0,dt);    % @odef is a construct handle for odef function
            X=[0 y(1:6:6*i)];
            Y=[0 y(3:6:6*i)];
            xGtl = y(6*i-5);
            yGtl = y(6*i-3);
            gieta= y(6*i); 
            X(end) = xGtl + Ltk*cos(gieta) + deltaY*sin(gieta);
            Y(end) = yGtl + Ltk*sin(gieta) - deltaY*cos(gieta);
        for ii=3:6:6*i
            if y(ii)<h0
                y(ii)=h0;
                y(ii-1)=h0;
                y(ii+1)=h0;
            end
        end

        set(h,'XData',X,'YData',Y);
        set(h1,'XData',X,'YData',Y);
        set(h2,'XData',X(end),'YData',Y(end));
        %X(end)
        plot(xGtl,yGtl,'.c','markersize',2);
        vtx = Ltk*cos(gieta);
        vty = Ltk*sin(gieta);
        set(h3,'XData',[xGtl-20*vtx xGtl+13*vtx],'YData',[yGtl-20*vty yGtl+13*vty]);
        drawnow
        L = sqrt(y(1)^2+y(3)^2);
        t0= t0+dt;
        y0= y;
        
        t(j)=t0;
        Gieta(j)=gieta;
        Mza(j) = Mz;
        Pda(j) = Pd;
        TTx0(j)=Tx(1);
        TTy0(j)=Ty(1);
        
        TTxm(j) = Tx(end-1);
        TTym(j) = Ty(end-1);
        
        XX(j)=y(end-5);
        YY(j)=y(end-3);
        Vx(j)=y(end-4);
        Vy(j)=y(end-2);
        Vx_goc=y(2);
        Vy_goc=y(4);
        if Vy(end)==0;
            C1=0;
            break
        end
    end
    
    m=[dl*mu m]; % m(end) = m dau dan
    
    y0=[0 y(2) 0 y(4) y(5) y(6) y0];
    if(Vy(end)==0);
        C1=0;
        break
    end
    
end

%%

k = k+1; % -> k = 50, m(k) = m(end)

zz= max(y(3:6:6*(k-1)));

%% Tinh phan nay khi da keo het day len

while zz>0  % Xac dinh khi day roi het xuong dat
     j=j+1;
     [y,Tx,Ty,Pd,Mz]=runKut4(@odef,t0,y0,dt);
            X=[0 y(1:6:6*k)];
            Y=[0 y(3:6:6*k)];
            xGtl = y(6*k-5);
            yGtl = y(6*k-3);
            gieta= y(6*k); 
            X(end) = xGtl + Ltk*cos(gieta) + deltaY*sin(gieta);
            Y(end) = yGtl + Ltk*sin(gieta) - deltaY*cos(gieta);
     for ii=3:6:6*k
        if y(ii)<0  
                y(ii)=h0; 
                y(ii-1)=h0;
                y(ii+1)=h0;
        end
     end
     set(h,'XData',X,'YData',Y);
     set(h1,'XData',X,'YData',Y);
     set(h2,'XData',X(end),'YData',Y(end));

     plot(xGtl,yGtl,'.c','markersize',2);
     vtx = Ltk*cos(gieta);
     vty = Ltk*sin(gieta);
     set(h3,'XData',[xGtl-20*vtx xGtl+13*vtx],'YData',[yGtl-20*vty yGtl+13*vty]);
     
     drawnow
     y0=y;   
     t0=t0+dt;
     zz=max(y(3:6:6*k));
     t(j)=t0;
     Gieta(j)=gieta;
     Mza(j) = Mz;
     Pda(j) = Pd;
     TTx0(j)=Tx(1);
     TTy0(j)=Ty(1);
     
     TTxm(j) = Tx(end-1);
     TTym(j) = Ty(end-1);
        
     XX(j)=y(end-5);
     YY(j)=y(end-3);
     Vx(j)=y(end-4);
     Vy(j)=y(end-2);
end

V     = sqrt(Vx.^2+Vy.^2);                % Van toc
T_goc = sqrt(TTx0.^2+TTy0.^2);   
T_m   = sqrt(TTxm.^2+TTym.^2); 
figure(2)
plot(t,V,t,Vx,t,Vy)
legend('V','Vx','Vy');
grid

figure(3)
plot(t,T_goc)
title('T tai diem O');

figure(4)
plot(t,T_m)
title('T tai diem buoc TL');

% disp(ti);
% figure(5)
% plot(t,Pda)
% title('P cua TL');
% figure(6)
% plot(t,Gieta)
% title('Gieta cua TL');
figure(6)
plot(t,Mza)
title('Mza t/d vao TL');
%%--------------------------------------------------------------
function [f,Tx,Ty,Pd,Mz]=odef(t,y)
global n md g k mu dl C1 K1 dt B k1 pk Fth ti pt P mk m pa Vb nu R Tdo Jz Ltk Lta Cp
global deltaY Ltkp alphaP ckd Ltk0
f       = zeros(1,6*k);
%thay trong tam dau TL bang diem chot o duoi TL
xGtl    = y(6*k-5);
yGtl    = y(6*k-3);
gieta   = y(6*k);
omega   = y(6*k-1);
Vx      = y(6*k-4);
Vy      = y(6*k-2);

xx      = y(1:6:6*k);
yy      = y(3:6:6*k);
xx(k)   = xGtl + Ltk*cos(gieta) + deltaY*sin(gieta);
yy(k)   = yGtl + Ltk*sin(gieta) - deltaY*cos(gieta);

xdif    = diff([-1e-9 xx]); %position
ydif    = diff([-1e-9 yy]);
l       = sqrt((xdif).^2+(ydif).^2);

xx2     = y(2:6:6*k);    %velocity
yy2     = y(4:6:6*k);
% thay boi diem buoc voi duoi TL
xx2(k)  = Vx - Ltk*sin(gieta)*omega+ deltaY*cos(gieta)*omega;
yy2(k)  = Vy + Ltk*cos(gieta)*omega+ deltaY*sin(gieta)*omega;
x2dif   = diff([0 xx2]);
y2dif   = diff([0 yy2]);

w       = -0.002 ; % van toc gio
Tkdx    = -ckd.*((xx2-w).^2+yy2.^2).*(ydif./l);
Tkdy    =  ckd.*((xx2-w).^2+yy2.^2).*(xdif./l);
 
% day khong chiu nen, suy ra
detl    = l-dl;
nn      = detl<0;
detl(nn)= 0;

% Luc dan hoi cua day theo x
Tdhx    = C1.*(xdif./l).*(detl/dl);
% Luc dan hoi cua day theo y
Tdhy    = C1.*(ydif./l).*(detl/dl); 


CNx     = K1*(xdif./l).*(xdif.*x2dif+ydif.*y2dif)./(dl*l.^2);
nn      = detl<0;
CNx(nn) = 0;
CNy     = K1*(ydif./l).*(xdif.*x2dif+ydif.*y2dif)./(dl*l.^2);
nn      = detl<0;
CNy(nn) = 0;
 
if (pt>pa)
    Cp = k1*(2/(k1+1))^(k1/(k1-1))*sqrt((k1+1)*(1-(pa/pk)^((k1-1)/k1))/(k1-1));
    P  = [P Cp*Fth*pt];
    ti = ti+dt;
    pt = pk*(1+B*ti)^(-2*k1/(k1-1));
    mk = nu*pt*Vb/(R*Tdo); 
    Ltk = Ltk0*(1 - pt/pk*0.25);
    m(end)=md+mk; 
    Pd = P(end);
else
    P  = [P 0];
    Pd = P(end);
end

Tx=Tdhx+CNx;
Ty=Tdhy+CNy;

kmu = 0.14; %vi dau cuoi cuon tu do
if ~(k==n & l(1)>=dl)
    Tx(1)=kmu*mu*y(2)*sqrt(y(2)^2+y(4)^2); 
    Ty(1)=kmu*mu*y(4)*sqrt(y(2)^2+y(4)^2);
end

% fim=atan2(ydif(end),xdif(end))
% Ban dau Tm la luc meski, sau do la hop luc cua luc cang va luc nhot?
Tmx=-Tx(end);
Tmy=-Ty(end);
Tm = [Tmx,Tmy,0]';

Tx(end+1)=0;
Ty(end+1)=0;


GV   = [Vx-w,Vy,0]';
GM  = [Ltk*cos(gieta),Ltk*sin(gieta),0]';
GH  = (GV'*GM)/Ltk^2*GM;
HV  = GV - GH;
Vxx = norm(GH);
Vyy = norm(HV);
Rx  = 3*(6.666666667e-8*Vxx^5-0.6916666667e-5*Vxx^4+0.2483333334e-3*Vxx^3+0.441666666e-3*Vxx^2+0.17e-1*Vxx);
Ry  = 2*(7.916666667e-8*Vyy^5-0.9291666667e-5*Vyy^4+0.3779166666e-3*Vyy^3+0.9979166664e-2*Vyy^2+0.299166667e-1*Vyy);
if(Vxx==0)
    ex = [0,0,0]';
else
    ex = GH/Vxx;
end
if(Vyy==0)
    ey = [0,0,0]';
else
    ey = -HV/Vyy;
end
GA  = [(Ltk-Lta)*cos(gieta),(Ltk-Lta)*sin(gieta),0]';
% Luc khi dong hoc
FA  = -Rx*ex-Ry*ey;

GM1  = [Ltk*cos(gieta) + deltaY*sin(gieta),Ltk*sin(gieta)- deltaY*cos(gieta),0]';

Mzv = cross(GA,FA)+cross(GM1,Tm);
Mz  = Mzv(3);


for i=1:k
    f(6*i-5)= y(6*i-4);
    f(6*i-4)=(Tx(i+1)-Tx(i)+Tkdx(i))/m(i); 
    f(6*i-3)= y(6*i-2);
    f(6*i-2)=(Ty(i+1)-Ty(i)+Tkdy(i))/m(i)-g;
%     f(6*i-1)=Mz/Jz;
    f(6*i)  =y(6*i-1);
end

    f(6*k-5)= y(6*k-4);
    f(6*k-4)= (P(end)*cos(gieta)+FA(1) + Tmx)/m(k);
    f(6*k-3)= y(6*k-2);
    f(6*k-2)= (P(end)*sin(gieta)+FA(2) + Tmy)/m(k)-g;
    f(6*k-1)= Mz/Jz;
    f(6*k)  = y(6*k-1);

     
    
%-------------------------------------------------------------
function [y,Tx,Ty,Pd,Mz] = runKut4(dEqs,tn,y,h)
if size(y,1) > 1 ; y = y'; end 
[K1,~,~]   = feval(dEqs,tn,y);
[K2,~,~]   = feval(dEqs,tn + h/2,y + h*K1/2);
[K3,~,~]   = feval(dEqs,tn + h/2,y + h*K2/2);
[K4,Tx,Ty,Pd,Mz] = feval(dEqs,tn + h,y + h*K3);
y = y + h*(K1 + 2*K2 + 2*K3 + K4)/6;
