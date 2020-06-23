
function tinhtoanPad()
        clc;

        %! Neu phuong trinh can bang luu luong khong thoa man
        %! thi thay doi ap suat hoc lom cho truoc.
        %? Neu phuong trinh momen khong thoa man, thi thay doi goc nghieng cua Pad.
        %? Sau do lai tinh toan lai tu dau. bao gom ca tinh toan luu luong.

        % Thiet lap cac thong so ban dau

            %? Tinh do day mang dau
        % vong lap ap suat vung lom pad1
            %? Tinh ap suat dong
            %? Tinh ap suat tinh
            %? Tinh luu luong
            %? Tinh momen
        %todo: Bo ket qua: alpha1,precess1 cho pad1.

        global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
        
        global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
        global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
        global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
        global QQ MM
        

        result = [];
        QQ = [];
        MM = [];


        % anpha1 = 0.02*pi/180;%! goc nghieng cua pad1
        % anpha2 = 0.02*pi/180;%! goc nghieng cua pad2
        % anpha3 = 0.02*pi/180;%! goc nghieng cua pad3 
        anpha1 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad1
        anpha2 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad2
        anpha3 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad3 

        %  Precess1= 1;  % !áp suất vùng lõm1
        % Precess2= 1.8;  % !áp suất vùng lõm2
        % Precess3= 0.5;  % !áp suất vùng lõm3
        Precess1=0.8:0.01:1;  % !áp suất vùng lõm1
        Precess2= 1.6:0.01:2.0;  % !áp suất vùng lõm2
        Precess3= 0.4:0.01:0.6;  % !áp suất vùng lõm3

         %* cac thong so ban dau
        l = 150e-3; % chieu dai o
        R = 100e-3; % ban kinh truc
        Rt=100.21e-3;% ban kinh trong cua pad
        Rn=115e-3; % ban kinh ngoai cua pad
        Rlr=100.09e-3; % ban kinh lap rap
        cl = Rt-R; % khe ho cua o truc ( tam truc va tam pad trung nhau)
        cb = Rlr-R; % khe ho cai dat ( tam truc va tam o trung nhau)
        nguy = 0.063; % do nhot cua dau
        dt = 2*R; % duong kinh truc
        N=2000;
        W=4000;
        lrz=2*Rn/l;
        %posi = c1/R; % ti le khe ho
        %M =0.57;%

        m = 160; % so luoi theo chu vi
        ld = l/dt; % ti le chieu dai tren duong kinh
        n = 90; % so luoi theo phuong doc truc
        %e1 =0.75e-04;% 0.1e-3; % do lech tam ban dau
        exilon = 0.05;%e1/c1; % ti le lech tam
        phibd = 0.5; % goc lech ban dau cua truc
        deltalanda = 2/n; % so luoi theo phuong doc truc

        % thong so pad1
        phi11 = 3*pi/180; % goc dau cua pad
        phi12 = 117*pi/180;% goc sau cua pad
        beta1 = (phi11+phi12)/2; % toa do cua diem xoay
        exilonn1=0;%1.5e-5;% do nang cua pad1
        cbtt1=cb-exilonn1;
        M1 = 1-(cbtt1/cl);
        
        deltaphi1 = (phi12-phi11)/m; %@ chia luoi theo chu vi
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



        % thong so pad2
        phi21 = 123*pi/180; % goc dau cua pad
        phi22 = 237*pi/180;% goc sau cua pad
        beta2 = (phi22+phi21)/2; % toa do cua diem xoay
        exilonn2=0;%7e-6;% do nang cua pad2
        cbtt2=cb-exilonn2;
        M2 = 1-(cbtt2/cl);
        capa2=Rn*anpha2/cl;
        deltaphi2 = (phi22-phi21)/m; %@ chia luoi theo chu vi
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



        % thong so pad3
        phi31 = 243*pi/180; % goc dau cua pad
        phi32 = 357*pi/180;% goc sau cua pad
        beta3 = (phi32+phi31)/2; % toa do cua diem xoay
        exilonn3=0;%3e-6;% do nang cua pad3
        cbtt3=cb-exilonn3;
        M3 = 1-(cbtt3/cl);
        capa3=Rn*anpha3/cl;
        deltaphi3 = (phi32-phi31)/m; %@ chia luoi theo chu vi
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

    global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
    global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s 
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




        % TINH PAD1
        for c=1:length(anpha1)
            
            capa1=Rn*anpha1(c)/cl;
            [phi1,phic1,phit1,h1,hc1,ht1,h1s,hc1s,ht1s] = tinh_do_day_mang_dau(phi11,deltaphi1);
            for k=1:length(Precess1)
                anpha1(c)
                Precess1(k)
                 tinh_ap_suat_dong(deltaphi1,Precess1(k),hc1,ht1,h1);
                 tinh_ap_suat_tinh(deltaphi1,Precess1(k),hc1s,ht1s,h1s);
                 Q = tinh_luu_luong(deltaphi1,h1,h1s);
                 if(Q)
                    M = tinh_mo_men(phi11,phi1, deltaphi1);
                    if(M)
                        result = [result anpha1(c) Precess1(k)];
                    end
                 end
            end
        end
        QQ
        MM
        result
                  
    end



function [phi1,phic1,phit1,h1,hc1,ht1,h1s,hc1s,ht1s] = tinh_do_day_mang_dau(phid,deltaphi)
    %  do day mang dau
    
    global  m n cl exilon phibd  beta1 exilonn1  M1 capa1 

    for i=1:m+1
        for j=1:n+1
            % phan dong
            phi1(i,j) = phid+(i-1)*deltaphi;
            phic1(i,j)= phid+(i-1+1/2)*deltaphi;
            phit1(i,j)= phid+(i-1-1/2)*deltaphi;
            h1(i,j) = 1+exilon*cos(phi1(i,j)-phibd)-M1*cos(phi1(i,j)-beta1)+capa1*sin(phid-phi1(i,j));
            hc1(i,j) = 1+exilon*cos(phic1(i,j)-phibd)-M1*cos(phic1(i,j)-beta1)+capa1*sin(phid-phic1(i,j));
            ht1(i,j) = 1+exilon*cos(phit1(i,j)-phibd)-M1*cos(phit1(i,j)-beta1)+capa1*sin(phid-phit1(i,j));
            % phan tinh
            h1s(i,j)  = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phid-phi1(i,j));
            hc1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phid-phic1(i,j));
            ht1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phid-phit1(i,j));
        end
    end
        
    end



function pDong = tinh_ap_suat_dong(deltaphi,Precess,hc1,ht1,h1)
        % ch??ng trinh tinh ap suat dong
        global  m n ld  deltalanda S1 T1 ERR1
        global p01 p1 a1 b1 c1 d1 e1 f1
        
         % ma tran he so ban dau cho phan hydrodynamic
        p01 = zeros(m+1,n+1);
        p1 = zeros(m+1,n+1);
        a1 = zeros(m+1,n+1);
        b1 = zeros(m+1,n+1);
        c1 = zeros(m+1,n+1);
        d1 = zeros(m+1,n+1);
        e1 = zeros(m+1,n+1);
        f1 = zeros(m+1,n+1);

        GAP1=1;
        k1=1;
        S1 = 0; 
        T1 = 0;
        while GAP1>ERR1
            k1 = k1+1;
            fx1=0;
            fy1=0;
            for i =1:m+1
                for j =1:n+1
                    if i==1||i==m+1||j==1||j==n+1
                        p1(i,j)=0;
                    elseif i>=80&&i<=82&&j>=45&&j<=47
                        p1(i,j)=Precess;    
                    else
                        a1(i,j)= hc1(i,j)^3;
                        b1(i,j)= ht1(i,j)^3;
                        c1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                        d1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                        e1(i,j)= a1(i,j)+b1(i,j)+c1(i,j)+d1(i,j);
                        f1(i,j)= 3*deltaphi*(hc1(i,j)-ht1(i,j));
                        p1(i,j)= (a1(i,j)*p01(i+1,j)+b1(i,j)*p01(i-1,j)+c1(i,j)*p01(i,j+1)+d1(i,j)*p01(i,j-1)-f1(i,j))/e1(i,j); % he phuong trinh tinh ap suat
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
    end


function pTinh = tinh_ap_suat_tinh(deltaphi, Precess,hc1s,ht1s,h1s)
        % chuong trinh tinh ap suat tinh

        global  m n lrz deltalanda S1s T1s ERR1s
        global p01s p1s a1s b1s c1s d1s e1s f1s 
        
        % ma tran he so cho phan hydrostatic
        p01s = zeros(m+1,n+1);
        p1s  = zeros(m+1,n+1);
        a1s  = zeros(m+1,n+1);
        b1s  = zeros(m+1,n+1);
        c1s  = zeros(m+1,n+1);
        d1s  = zeros(m+1,n+1);
        e1s  = zeros(m+1,n+1);
        f1s  = zeros(m+1,n+1);
        
        
        GAP1s=1;
        k1s=1;
        S1s = 0; 
        T1s = 0;
        while GAP1s>ERR1s
            k1s = k1s+1;
            for i =1:m+1
                for j =1:n+1
                    if j==1||j==n+1||i==1||i==m+1
                        p1s(i,j)=0;
                        %p(i,j)=0;%p(i,j+1);
                    elseif i>=65&&i<=95&&j>=35&&j<=55
                        p1s(i,j)=Precess;
                    else
                        a1s(i,j)= hc1s(i,j)^3;
                        b1s(i,j)= ht1s(i,j)^3;
                        c1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
                        d1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
                        e1s(i,j)= a1s(i,j)+b1s(i,j)+c1s(i,j)+d1s(i,j);
                        f1s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                        p1s(i,j)= (a1s(i,j)*p01s(i+1,j)+b1s(i,j)*p01s(i-1,j)+c1s(i,j)*p01s(i,j+1)+d1s(i,j)*p01s(i,j-1)-f1s(i,j))/e1s(i,j);  % he phuong trinh tinh ap suat
                        
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
end



function Q = tinh_luu_luong(deltaphi,h1,h1s)
        % Can bang luu luong

        global  n  ld deltalanda p1 p1s
        global QQ


        Qd11=0;
        Qd12=0;
        Qd13=0;
        Qs11=0;
        Qs12=0;
        Qs13=0;
        for j=(round(n/2-1):(n/2+2))

            Qd11=Qd11+(1/2)*ld^2*((h1(79,44)-(1/3)*(h1(79,44))^3*((-3*p1(78,j)+4*p1(79,j)-p1(80,j)))/(2*deltaphi))*(deltalanda));
            %Qd11=Qd11+(1/2)*ld^2*((h1(79,44)+(1/3)*(h1(79,44))^3*((-3*p1(79,j)+4*p1(78,j)-p1(77,j)))/(2*deltaphi))*(deltalanda));
            Qd12=Qd12+(1/2)*ld^2*((-h1(83,44)-(1/3)*(h1(83,44))^3*((3*p1(84,j)-4*p1(83,j)+p1(82,j)))/(2*deltaphi))*(deltalanda));
            %Qd12=Qd12+(1/4)*ld^2*((h1(83,44)-(1/3)*(h1(83,44))^3*((3*p1(83,j)-4*p1(84,j)+p1(85,j)))/(2*deltaphi))*(deltalanda));


        end
        for j=44
            for i=79:83
                Qd13=Qd13-(1/6)*(h1(i,44))^3*(3*p1(i,44)-4*p1(i,43)+p1(i,42))*(deltaphi/deltalanda);
            end
        end
        Q1d=Qd11+Qd12+Qd13;
        
        for i=51
            for j=16:76
            Qs11=Qs11-(1/2)*ld^2*((1/3)*(h1s(51,16))^3*((-3*p1s(51,j)+4*p1s(50,j)-p1s(49,j))/(2*deltaphi))*(deltalanda)); 
            Qs12=Qs12-(1/2)*ld^2*(-(1/3)*(h1s(111,16))^3*((3*p1s(111,j)-4*p1s(112,j)+p1s(113,j))/(2*deltaphi))*(deltalanda));

            end
        end
        for j=16
            for i=51:111
                Qs13=Qs13+(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j-1)+p1s(i,j-2))*(deltaphi/deltalanda);
                %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
            end
        end
        Q1s=Qs11+Qs12+Qs13;
        
        QQ = [QQ Q1d Q1s];
        % KIEM TRA CAN BANG LUU LUONG.
        if(abs(Q1d-Q1s)<=0.001)
            Q=true;
        else
            Q=false;
        end
        
    end



    
function M = tinh_mo_men(phid,phi1, deltaphi)

    global  m n R Rn  deltalanda p1 p1s MM
     % can bang momen
        Md1=0;
        Msr1=0;
        Mso11=0;
        Mso21=0;
        Mso31=0;
        Mso41=0;
        for i=1:m+1
            for j=1:n+1
                Md1=Md1-R*p1(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
        end
        for i=51:111
            for j=16:76
                Msr1=Msr1-Rn*p1s(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
        end
        for i=1:50
            for j=1:n+1
                Mso11=Mso11-Rn*p1s(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
        end
        for i=112:m+1
            for j=1:n+1
                Mso21=Mso21-Rn*p1s(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
        end
        for i=51:111
            for j=1:15
                Mso31=Mso31-Rn*p1s(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
        end
            for i=51:111
            for j=77:n+1
                Mso41=Mso41-Rn*p1s(i,j)*sin(phid-phi1(i,j))*deltaphi*deltalanda;
            end
            end
        Ms1=Msr1+Mso11+Mso11+Mso31+Mso41;
        MM = [MM Md1 Ms1];
        if(abs(Md1-Ms1)<=0.001)
            M=true;
        else
            M=false;
        end
        
    end

