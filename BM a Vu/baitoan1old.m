

function baitoan
    clc
    clear
    global  e c cb nguy l R d posi M ld m n exilon coxi 
    e = 0.012e-3; % do lech tam ban dau
    c = 0.2e-3; % khe ho cua o truc ( tam truc va tam pad trung nhau)
    cb = 0.06e-3; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    l = 0.1;      % chieu dai o
    R = 0.1;      % ban kinh truc
    d = 2*R;      % duong kinh truc
    posi = c/R;   % ti le khe ho
    M = 1 -(cb/c);
    ld = l/d;     % ti le chieu dai tren duong kinh
    m = 40;       % so luoi theo chu vi
    n = 30;       % so luoi theo phuong doc truc

    exilon = e/c; % ti le lech tam
    mi1=0;
    mi2=0;
    mi3=0;
    
    coxi = 0.05;
    xic1= 0;
    xic2 = 0;
    xic3 = 0;
    %phibdx=0;
    
   % mi=[];
    
    %phibd = 0.5:0.1:2.9;  % goc lech ban dau cua truc
    phibd = 1.9;  % goc lech ban dau cua truc
    xicma = 0.00001:0.00001:0.0004; % goc lech ban dau cua pad
    %for i=1:length(phibd)
    while xic1==0
        for j=1:length(xicma)
            xic1 = tinhtoanpad(phibd,xicma(j),0,120);
        end
    end
    
    while xic2==0
        for j=1:length(xicma)
            xic2 = tinhtoanpad(phibd,xicma(j),120,240);
        end
    end
    
    while xic3==0
        for j=1:length(xicma)
            xic3 = tinhtoanpad(phibd,xicma(j),240,360);
        end
    end
    %end
    xic1
    xic2
    xic3
    %plot(xicma,mi)
end

function xic = tinhtoanpad(phibd,xicma,alfa,beta)
        global  e c d posi M ld m n exilon coxi
        %mi1=0;

        % tinh cho pad dau tien
        phi11 = alfa; % goc dau cua pad
        phi12 = beta;% goc sau cua pad
        beta1 = (phi12+phi11)/2; % toa do cua diem xoay

        fi1=0;
        fi2=0;
        % tinh cac thong so ban dau
        deltalanda = 1/n; % so luoi theo phuong doc truc
        deltaphi1  = (phi12-phi11)/m; % chia luoi theo chu vi
        % khai bao cac ma tran he so ban dau
        p0 = zeros(m+1,n+1);
        p  = zeros(m+1,n+1);
        a  = zeros(m+1,n+1);
        b  = zeros(m+1,n+1);
        c  = zeros(m+1,n+1);
        d  = zeros(m+1,n+1);
        e  = zeros(m+1,n+1);
        f  = zeros(m+1,n+1);
        %h  = zeros(m+1,n+1);
        % tinh do day mang dau

        for i=1:m+1
            for j=1:n+1
                phi(i,j) = phi11+(i-1)*deltaphi1;
                phic(i,j)= phi11+(i-1+1/2)*deltaphi1;
                phit(i,j)= phi11+(i-1-1/2)*deltaphi1;
                h(i,j)   = 1+exilon*cosd(phi(i,j)-phibd)-M*cosd(phi(i,j)-beta1)-(xicma/posi)*sind(phi(i,j)-beta1);
                hc(i,j)   = 1+exilon*cosd(phic(i,j)-phibd)-M*cosd(phic(i,j)-beta1)-(xicma/posi)*sind(phic(i,j)-beta1);
                ht(i,j)   = 1+exilon*cosd(phit(i,j)-phibd)-M*cosd(phit(i,j)-beta1)-(xicma/posi)*sind(phit(i,j)-beta1);
            end
        end
        % chuong trinh tinh ap suat
        % khai bao cac bien ban dau
        S = 0; 
        T = 0;
        ERR = 10e-3;
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
                        c(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        d(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                        f(i,j)= 3*deltaphi1*(hc(i,j)-ht(i,j));
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
        
        for i=1:m+1
            for j=1:n+1
                mi1= mi1 +p(i,j)*sind(beta1-phi(i,j))*deltaphi1*deltalanda;
            end
        end

        if(mi1<=coxi)
            xic = xicma;
        
        else
            xic = 0;
        end
            
    end





        
              %end
% VE DO THI
% x=1:1:41;
% y=1:1:31;
% z=p(x,y);
% figure
% mesh(z);