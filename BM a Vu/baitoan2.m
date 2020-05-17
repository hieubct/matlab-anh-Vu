function baitoan1
    clc
    clear
    
    SoR = zeros(4,10);
    count=1;
    Fx = 0;
    Fy = 0;
    xic1= 0;
    xic2 = 0;
    xic3 = 0;
    %phibdx=0;
    
   % mi=[];
    
    %phibd = 0.5:0.1:2.9;  % goc lech ban dau cua truc
    phibd = 2.95:0.01:3.05;  % goc lech ban dau cua truc
    xicma = 0.0007:0.000001:0.0011; % goc lech  ban dau cua pad
    
    for i=1:length(phibd)
        if(xic1==0 || xic2 ==0 || xic3 == 0)
            xic1= 0;
            xic2 = 0;
            xic3 = 0;
        end
               
        for j=1:length(xicma)
            if xic1==0
                xic1 = tinhtoanpad(phibd(i),xicma(j),0,120);
            end
        
            if (xic2==0)
                xic2 = tinhtoanpad(phibd(i),xicma(j),120,240);
            end
            
        
            if (xic3==0)
                xic3 = tinhtoanpad(phibd(i),xicma(j),240,360);
            end
        end
           
        if(xic1~=0 && xic2~=0 && xic2 ~=0)
            [Fx1 Fy1]=fxPfy(phibd(i),xic1,0,120);
            [Fx2 Fy2]=fxPfy(phibd(i),xic2,120,240);
            [Fx3 Fy3]=fxPfy(phibd(i),xic3,240,360);
            Fx = Fx1+Fx2+Fx3;
            Fy = Fy1+Fy2+Fy3;
            if(Fx/Fy <0.0001)
                %Thu thap ket qua vào bo du lieu. 
                SoR(1,count)=xic1;
                SoR(2,count)=xic2;
                SoR(3,count)=xic3;
                SoR(4,count)=phibd(i);
                count=count+1;
                xic1= 0;
                xic2 = 0;
                xic3 = 0;
            end
        end
        
        
    end
    
%     xic1
%     xic2
%     xic3
%     phibd(i)
    %plot(xicma,mi)
%     [Fx1 Fy1]=fxPfy(phibd,xic1,0,120);
%     [Fx2 Fy2]=fxPfy(phibd,xic2,120,240);
%     [Fx3 Fy3]=fxPfy(phibd,xic3,240,360);
%     Fx = Fx1+Fx2+Fx3
%     Fy = Fy1+Fy2+Fy3
%     Fx1
%     Fx2
%     Fx3
%     Fy1
%     Fy2
%     Fy3
%     Fx/Fy
        SoR
end

function xic = tinhtoanpad(phibd,xicma,alfa,beta)
        %global mi1 e c cb nguy l R d posi M ld m n exilon coxi xic1
        mi1=0;
        coxi = 0.001;

        % tinh cho pad dau tien
       phi11=alfa; % goc dau cua pad
       phi12=beta;% goc sau cua pad
        beta1 = (phi12+phi11)/2; % toa do cua diem xoay

        fi1=0;
        fi2=0;
        n=30;
        m=40;
        c = 0.2e-3;
         e = 0.012e-3; % do lech tam ban dau
         exilon = e/c;
    
    
    cb = 0.06e-3; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    l = 0.1;      % chieu dai o
    R = 0.1;      % ban kinh truc
    d = 2*R;      % duong kinh truc
    %posi = c/R;   % ti le khe ho
    M = 1 -(cb/c);
    posi = c/R;
    ld = l/d;   
        % tinh cac thong so ban dau
        deltalanda = 1/n; % so luoi theo phuong doc truc
        deltaphi1  = (phi12-phi11)/m; % chia luoi theo chu vi
        % khai bao cac ma tran he so ban dau
        P0 = zeros(m+1,n+1);
        P  = zeros(m+1,n+1);
        A  = zeros(m+1,n+1);
        B  = zeros(m+1,n+1);
        C  = zeros(m+1,n+1);
        D  = zeros(m+1,n+1);
        E  = zeros(m+1,n+1);
        F  = zeros(m+1,n+1);
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
                        P(i,j)=0;
                    else
                        A(i,j)= hc(i,j)^3;
                        B(i,j)= ht(i,j)^3;
                        C(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        D(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        E(i,j)= A(i,j)+B(i,j)+C(i,j)+D(i,j);
                        F(i,j)= 3*deltaphi1*(hc(i,j)-ht(i,j));
                        P(i,j)= (A(i,j)*P0(i+1,j)+B(i,j)*P0(i-1,j)+C(i,j)*P0(i,j+1)...  % he phuong trinh tinh ap suat
                            +D(i,j)*P0(i,j-1)-F(i,j))/E(i,j);
                        if P(i,j)<=0
                            P(i,j)=0;
                        end
                    end
                 end
            end
            % Kiem tra dieu kien hoi tu
            for i= 1:m+1
                for j= 1:n+1
                    S=S+abs(P(i,j)-P0(i,j));
                    T=T+abs(P(i,j));
                    GAP=S/T;
                end
            end
            P0=P;
        end
        
        for i=1:m+1
            for j=1:n+1
                mi1= mi1 +P(i,j)*sind(beta1-phi(i,j))*deltaphi1*deltalanda;
            end
        end

        if(mi1<=coxi)
            xic = xicma;
        
        else
            xic = 0;
        end
          % x=1:1:41;
% y=1:1:31;
% z=p(x,y);
 %figure
 %mesh(z);  
end


function [Fx Fy]= fxPfy(phibd,xicma,alfa,beta)
        %global mi1 e c cb nguy l R d posi M ld m n exilon coxi xic1
        %global Fx Fy
        mi1=0;
        coxi = 0.001;
        Fx=0;
        Fy=0;

        % tinh cho pad dau tien
       phi11=alfa; % goc dau cua pad
       phi12=beta;% goc sau cua pad
        beta1 = (phi12+phi11)/2; % toa do cua diem xoay

        fi1=0;
        fi2=0;
        n=30;
        m=40;
        c = 0.2e-3;
         e = 0.012e-3; % do lech tam ban dau
         exilon = e/c;
    
    
    cb = 0.06e-3; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    l = 0.1;      % chieu dai o
    R = 0.1;      % ban kinh truc
    d = 2*R;      % duong kinh truc
    %posi = c/R;   % ti le khe ho
    M = 1 -(cb/c);
    posi = c/R;
    ld = l/d;   
        % tinh cac thong so ban dau
        deltalanda = 1/n; % so luoi theo phuong doc truc
        deltaphi1  = (phi12-phi11)/m; % chia luoi theo chu vi
        % khai bao cac ma tran he so ban dau
        P0 = zeros(m+1,n+1);
        P  = zeros(m+1,n+1);
        A  = zeros(m+1,n+1);
        B  = zeros(m+1,n+1);
        C  = zeros(m+1,n+1);
        D  = zeros(m+1,n+1);
        E  = zeros(m+1,n+1);
        F  = zeros(m+1,n+1);
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
                        P(i,j)=0;
                    else
                        A(i,j)= hc(i,j)^3;
                        B(i,j)= ht(i,j)^3;
                        C(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        D(i,j)= ld^2*(deltaphi1/deltalanda)^2*h(i,j)^3;
                        E(i,j)= A(i,j)+B(i,j)+C(i,j)+D(i,j);
                        F(i,j)= 3*deltaphi1*(hc(i,j)-ht(i,j));
                        P(i,j)= (A(i,j)*P0(i+1,j)+B(i,j)*P0(i-1,j)+C(i,j)*P0(i,j+1)...  % he phuong trinh tinh ap suat
                            +D(i,j)*P0(i,j-1)-F(i,j))/E(i,j);
                        if P(i,j)<=0
                            P(i,j)=0;
                        end
                    end
                 end
            end
            % Kiem tra dieu kien hoi tu
            for i= 1:m+1
                for j= 1:n+1
                    S=S+abs(P(i,j)-P0(i,j));
                    T=T+abs(P(i,j));
                    GAP=S/T;
                end
            end
            P0=P;
        end
        
        for i= 1:m+1
           for j=1:n+1
               Fx=Fx+P(i,j)*sind(phi(i,j))*deltaphi1*deltalanda;% luc theo truc x
               Fy=Fy+P(i,j)*cosd(phi(i,j))*deltaphi1*deltalanda;% luc theo truc y
           end
        end

          % x=1:1:41;
% y=1:1:31;
% z=p(x,y);
 %figure
 %mesh(z);  
end