function chuongtrinhtinhlucmangdau
    clc;    
    global xicma_nghiem FxFy
    xicma_nghiem = [];
    phi = [0 120 120 240 240 360];
    FxFy = [];
    xicma = [-0.0005 0.0005 0 0.001 0 0.005];
    
    tinhtoan(xicma(1),xicma(2),phi(1), phi(2));
    tinhtoan(xicma(3),xicma(4),phi(3), phi(4));
    tinhtoan(xicma(5),xicma(6),phi(5), phi(6));

    xicma_nghiem
    FxFy
    FX=FxFy(1) + FxFy(3) + FxFy(5)
    FY=FxFy(2) + FxFy(4) + FxFy(6)
end

function tinhtoan(xicmaa,xicmab, phid, phic)
    global xicma_nghiem FxFy
    xicmac = (xicmaa + xicmab)/2;
    a=hamgocnghieng(xicmaa, phid, phic);
    b=hamgocnghieng(xicmab, phid, phic);
    c=hamgocnghieng(xicmac, phid, phic);

    if(a*b>0)
        disp('Vo nghiem trong khoang da cho');
        return
    end

    while abs(c)>0.001
        if(c*a < 0)
            xicmab = xicmac;
        elseif(c*b< 0)
            xicmaa = xicmac;
        else
            disp('Vo nghiem trong khoang da cho');
            break;
        end
        xicmac = (xicmaa + xicmab)/2;
        a=hamgocnghieng(xicmaa, phid, phic);
        b=hamgocnghieng(xicmab, phid, phic);
        c=hamgocnghieng(xicmac, phid, phic);
    end
    
    xicma_nghiem = [ xicmac];
    [Fx,Fy] = hamtinhluc(xicmac, phid, phic)
    FxFy = [ Fx Fy];
    
end


%% ham tinh goc nghieng cua pad
function mi1 = hamgocnghieng(xicma, phi1, phi2)
    e = 0.1e-3; % do lech tam ban dau
    c = 0.254e-3; % khe ho cua o truc ( tam truc va tam pad trung nhau)
    cb = 0.1905e-3; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    l = 76.2e-3; % chieu dai o
    R = 116.8095e-3; % ban kinh truc
    d = 2*R; % duong kinh truc
    posi = c/R; % ti le khe ho
    M = 1 -(cb/c);
    ld = l/d; % ti le chieu dai tren duong kinh
    m = 160; % so luoi theo chu vi
    n = 90; % so luoi theo phuong doc truc
    exilon = e/c; % ti le lech tam
    coxi = 0.0005;
    phibd = 1.9*pi/180; % goc lech ban dau cua truc
    % % tinh cac thong so ban dau
    phi21 = phi1*pi/180; % goc dau cua pad
    phi22 = phi2*pi/180;% goc sau cua pad
    beta2 = (phi22+phi21)/2; % toa do cua diem xoay
    deltalanda = 1/n; % so luoi theo phuong doc truc
    deltaphi1 = (phi22-phi21)/m; % chia luoi theo chu vi
    % khai bao cac ma tran he so ban dau
    p0 = zeros(m+1,n+1);
    p = zeros(m+1,n+1);
    a = zeros(m+1,n+1);
    b = zeros(m+1,n+1);
    c = zeros(m+1,n+1);
    d = zeros(m+1,n+1);
    e = zeros(m+1,n+1);
    f = zeros(m+1,n+1);
    %h = zeros(m+1,n+1);
    % tinh do day mang dau
    for i=1:m+1
        for j=1:n+1
            phi(i,j) = phi21+(i-1)*deltaphi1;
            phic(i,j)= phi21+(i-1+1/2)*deltaphi1;
            phit(i,j)= phi21+(i-1-1/2)*deltaphi1;
            h(i,j) = 1+exilon*cos(phi(i,j)-phibd)-M*cos(phi(i,j)-beta2)-(xicma/posi)*sin(phi(i,j)-beta2);
            hc(i,j) = 1+exilon*cos(phic(i,j)-phibd)-M*cos(phic(i,j)-beta2)-(xicma/posi)*sin(phic(i,j)-beta2);
            ht(i,j) = 1+exilon*cos(phit(i,j)-phibd)-M*cos(phit(i,j)-beta2)-(xicma/posi)*sin(phit(i,j)-beta2);
        end
    end
    % chuong trinh tinh ap suat
    % khai bao cac bien ban dau
    S = 0; 
    T = 0;
    ERR = 10e-4;
    GAP =1;
    k = 1; % he so lap
    %vong lap dk tinh ap suat
    while GAP>ERR
        k = k+1;
        mi1=0;
        for i =1:m+1
            for j =1:n+1
                if i==1||i==m+1||j==1||j==n+1
                    p(i,j)=0;
                else
                    a(i,j)= hc(i,j)^3;
                    b(i,j)= ht(i,j)^3;
                    c(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    d(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                    f(i,j)= 3*deltaphi1*(hc(i,j)-ht(i,j));
                    p(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)... % he phuong trinh tinh ap suat
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
        for i=1:m+1
            for j=1:n+1
                mi1= mi1 +p(i,j)*sin(beta2-phi(i,j))*deltaphi1*deltalanda;
            end
        end
    end 
end

%% % ham tinh luc ap suat mang dau1
function [Fx,Fy] = hamtinhluc(xicma, phi1, phi2)
    e = 0.1e-3; % do lech tam ban dau
    c = 0.254e-3; % khe ho cua o truc ( tam truc va tam pad trung nhau)
    cb = 0.1905e-3; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    l = 76.2e-3; % chieu dai o
    R = 116.8095e-3; % ban kinh truc
    d = 2*R; % duong kinh truc
    posi = c/R; % ti le khe ho
    M = 1 -(cb/c);
    ld = l/d; % ti le chieu dai tren duong kinh
    m = 160; % so luoi theo chu vi
    n = 90; % so luoi theo phuong doc truc
    exilon = e/c; % ti le lech tam
    coxi = 0.0005;
    phibd = 1.9*pi/180; % goc lech ban dau cua truc
    Fx=0;
    Fy=0;
    % % tinh cac thong so ban dau
    phi21 = phi1*pi/180; % goc dau cua pad
    phi22 = phi2*pi/180;% goc sau cua pad
    beta2 = (phi22+phi21)/2; % toa do cua diem xoay
    deltalanda = 1/n; % so luoi theo phuong doc truc
    deltaphi1 = (phi22-phi21)/m; % chia luoi theo chu vi
    % khai bao cac ma tran he so ban dau
    p0 = zeros(m+1,n+1);
    p = zeros(m+1,n+1);
    a = zeros(m+1,n+1);
    b = zeros(m+1,n+1);
    c = zeros(m+1,n+1);
    d = zeros(m+1,n+1);
    e = zeros(m+1,n+1);
    f = zeros(m+1,n+1);
    %h = zeros(m+1,n+1);
    % tinh do day mang dau
    for i=1:m+1
        for j=1:n+1
            phi(i,j) = phi21+(i-1)*deltaphi1;
            phic(i,j)= phi21+(i-1+1/2)*deltaphi1;
            phit(i,j)= phi21+(i-1-1/2)*deltaphi1;
            h(i,j) = 1+exilon*cos(phi(i,j)-phibd)-M*cos(phi(i,j)-beta2)-(xicma/posi)*sin(phi(i,j)-beta2);
            hc(i,j) = 1+exilon*cos(phic(i,j)-phibd)-M*cos(phic(i,j)-beta2)-(xicma/posi)*sin(phic(i,j)-beta2);
            ht(i,j) = 1+exilon*cos(phit(i,j)-phibd)-M*cos(phit(i,j)-beta2)-(xicma/posi)*sin(phit(i,j)-beta2);
        end
    end
    % chuong trinh tinh ap suat
    % khai bao cac bien ban dau
    S = 0; 
    T = 0;
    ERR = 10e-4;
    GAP =1;
    k = 1; % he so lap
    %vong lap dk tinh ap suat
    while GAP>ERR
        k = k+1;
        for i =1:m+1
            for j =1:n+1
                if i==1||i==m+1||j==1||j==n+1
                    p(i,j)=0;
                else
                    a(i,j)= hc(i,j)^3;
                    b(i,j)= ht(i,j)^3;
                    c(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    d(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                    f(i,j)= 3*deltaphi1*(hc(i,j)-ht(i,j));
                    p(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)... % he phuong trinh tinh ap suat
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
    for i=1:m+1
        for j=1:n+1
            Fx=Fx+p(i,j)*cos(phi21+(i-1)*deltaphi1)*deltaphi1*deltalanda;
            Fy=Fy+p(i,j)*sin(phi22+(i-1)*deltaphi1)*deltaphi1*deltalanda;
        end
    end
end 


