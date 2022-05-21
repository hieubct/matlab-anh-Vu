           %PROGRAM TPN_9M22Y
clear all;
clc
%--cac tham so cua dan---
d=0.122; %duong kinh dan:m
tk=1.89; %thoi gian lam viec cua dong co: s
Ue=1900; %toc do phut hieu dung: m/s
g=9.81; %gia toc trong truong: m/s2
a=341; %van toc am: m/s
tooc=288.9;% nhiet do tieu chuan tren mat dat
omega=20.45; %khoi luong thuoc phong:
Qo=input('nhap khoi luong ban dau cua dan=');
G=6.328*10^-3; %gradian nhiet do
i=input('he so hinh dang dan=');%0.96
Co=(i*d^2/Qo)*10^3;
h=input('buoc tinh tich phan=');%0.05
     %DIEU KIEN DAU'
t0=0;y0=0;x0=0;v0=41;tetak0=0.785;
%CHUONG TRINH CHINH GIAI HE PHUONG TRINH BANG TICH PHAN SO
t=t0;y=y0;x=x0;v=v0;teta=tetak0;kq=[];
vmax=v;ymax=y;xmax=x;
while y >=0
    Mb=[0.0001 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4];
    Cxb=[0.305 0.305 0.306 0.308 0.308 0.308 0.316 0.333 0.382 0.551 0.616 0.618 0.605 0.578 0.559 0.538 0.521 0.506 0.489 0.477 0.462 0.449 0.439 0.426 0.416 0.406 0.394 0.386 0.376 0.367 0.361 0.352 0.346 0.338 0.331 0.325 0.318 0.312 0.308 0.320];
    Yb=[0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5000 5100 5200 5300 5400 5500 5600 5700 5800 5900 6000 6100 6200 6300 6400 6500 6600 6700 6800 6900 7000 15000];
    pyb=[1 0.9881 0.9764 0.9649 0.9546 0.9423 0.9311 0.92 0.9090 0.8982 0.887 0.8766 0.866 0.8555 0.845 0.835 0.8249 0.8148 0.8048 0.7949 0.7852 0.7756 0.766 0.7565 0.7471 0.7379 0.7288 0.7196 0.7106 0.7017 0.6928 0.6841 0.6754 0.6670 0.6585 0.65 0.6418 0.6336 0.6255 0.6175 0.6095 0.6017 0.5939 0.5861 0.5785 0.571 0.5635 0.5562 0.5488 0.5417 0.5346 0.5276 0.5206 0.5137 0.5069 0.5 0.4934 0.4869 0.4801 0.4737 0.4673 0.4608 0.4546 0.4484 0.4423 0.4363 0.4304 0.4245 0.4187 0.4130 0.4072 0.125];
    if t<=tk
        miu=omega*t/(Qo*tk);
        anfa=1;beta=1;
    else
        miu=omega/Qo;
        anfa=0;beta=1;
    end
    to=tooc-G*y;M=v/a;Cx=interp1(Mb,Cxb,M);Fvto=4.74*(10^(-4))*(v^2)*(tooc/to)*Cx;
    py=interp1(Yb,pyb,y);
    %LAN 1
    y1=y;x1=x;teta1=teta;v1=v;
    dv1=(anfa*omega*Ue/(Qo*tk)-beta*Co*py*Fvto)/(1-miu)-g*sin(teta1);d1v=h*dv1;
    dteta1=-g*cos(teta1)/v1;d1teta=h*dteta1;
    dx1=v1*cos(teta1);d1x=h*dx1;
    dy1=v1*sin(teta1);d1y=h*dy1;
    %LAN 2
    y2=y+d1y/2;x2=x+d1x/2;teta2=teta+d1teta/2;v2=v+d1v/2;
    dv2=(anfa*omega*Ue/(Qo*tk)-beta*Co*py*Fvto)/(1-miu)-g*sin(teta2);d2v=h*dv2;
    dteta2=-g*cos(teta2)/v2;d2teta=h*dteta2;
    dx2=v2*cos(teta2);d2x=h*dx2;
    dy2=v2*sin(teta2);d2y=h*dy2;
    %LAN 3
    y3=y+d2y/2;x3=x+d2x/2;teta3=teta+d2teta/2;v3=v+d2v/2;
    dv3=(anfa*omega*Ue/(Qo*tk)-beta*Co*py*Fvto)/(1-miu)-g*sin(teta3);d3v=h*dv3;
    dteta3=-g*cos(teta3)/v3;d3teta=h*dteta3;
    dx3=v3*cos(teta3);d3x=h*dx3;
    dy3=v3*sin(teta3);d3y=h*dy3;
    %LAN 4
    y4=y+d3y/2;x4=x+d3x/2;teta4=teta+d3teta/2;v4=v+d3v/2;
    dv4=(anfa*omega*Ue/(Qo*tk)-beta*Co*py*Fvto)/(1-miu)-g*sin(teta4);d4v=h*dv4;
    dteta4=-g*cos(teta4)/v4;d4teta=h*dteta4;
    dx4=v4*cos(teta4);d4x=h*dx4;
    dy4=v4*sin(teta4);d4y=h*dy4;
    %TINH CAC GIA TRI DENTA;
    deltav=(d1v+2*d2v+2*d3v+d4v)/6;
    deltateta=(d1teta+2*d2teta+2*d3teta+d4teta)/6;
    deltay=(d1y+2*d2y+2*d3y+d4y)/6;
    deltax=(d1x+2*d2x+2*d3x+d4x)/6;
    %LAP VONG CHAY;
    v=v+deltav;
    teta=teta+deltateta;
    y=y+deltay;
    x=x+deltax;
    t=t+h;
    ketqua=[t  v  x  y  teta*180/pi];
    kq=[kq;ketqua]
    if v>=vmax
        vmax=v;
    end
    if y>=ymax
        ymax=y;
    end
    if x>=xmax
        xmax=x;
    end
end
Xmax=xmax
Ymax=ymax
Vmax=vmax
tetac=abs(teta)
vcham=v
%VE DO THI BIEU DIEN SU PHU THUOC GIUA VAN TOC VA THOI GIAN
V=kq(:,2);T=kq(:,1);X=kq(:,3);Y=kq(:,4);Teta=kq(:,5);
plot(X,Y);
grid('on');
Xlabel('TAM XA');
Ylabel('DO CAO');
figure
plot(T,V)
grid on;
Xlabel('THOI GIAN(s)');
Ylabel('VAN TOC(m/s)');