

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


         anpha1 = 0.02*pi/180;%! goc nghieng cua pad1
        % anpha2 = 0.02*pi/180;%! goc nghieng cua pad2
        % anpha3 = 0.02*pi/180;%! goc nghieng cua pad3 
        % anpha1 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad1
        anpha2 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad2
        anpha3 = 0.01*pi/180:(0.001*pi/180):0.03*pi/180;%! goc nghieng cua pad3 

         Precess1= 1;  % !áp suất vùng lõm1
        % Precess2= 1.8;  % !áp suất vùng lõm2
        % Precess3= 0.5;  % !áp suất vùng lõm3
        % Precess1=0.8:0.01:1.2;  % !áp suất vùng lõm1
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
                 tinh_ap_suat_dong(deltaphi1,Precess1(k),hc1,ht1,h1);
                 tinh_ap_suat_tinh(deltaphi1,Precess1(k),hc1s,ht1s,h1s);
                 Q = tinh_luu_luong(deltaphi1,h1,h1s);
                 M = tinh_mo_men(phi11,phi1, deltaphi1);
                 if(Q&&M)
                     result = [result anpha1(c) Precess1(k)];
                 end
            end
        end
        QQ
        MM
        result

        % % TINH PAD2
        % for c=1:size(anpha2)
        %     capa2=Rn*anpha2(c)/cl;
        %     [phi1,phic1,phit1,h1,hc1,ht1,h1s,hc1s,ht1s] = tinh_do_day_mang_dau(phi21,deltaphi2);
        %     % tinh_ap_suat_dong(deltaphi1);
        %     % tinh_ap_suat_tinh(deltaphi1);
        %     % check = tinh_luu_luong(deltaphi1);
        %     % while(!check)
        %     %     check = tinh_luu_luong(deltaphi1);
        %     % end
        %     % tinh_mo_men();
        % end


        % % TINH PAD2
        % for c=1:size(anpha3)
        %     capa3=Rn*anpha3(c)/cl;
        %     [phi1,phic1,phit1,h1,hc1,ht1,h1s,hc1s,ht1s] = tinh_do_day_mang_dau(phi31,deltaphi3);
        %     % tinh_ap_suat_dong(deltaphi1);
        %     % tinh_ap_suat_tinh(deltaphi1);
        %     % check = tinh_luu_luong(deltaphi1);
        %     % while(!check)
        %     %     check = tinh_luu_luong(deltaphi1);
        %     % end
        %     % tinh_mo_men();
        % end

            
    end










function [phi1,phic1,phit1,h1,hc1,ht1,h1s,hc1s,ht1s] = tinh_do_day_mang_dau(phid,deltaphi)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    %  do day mang dau
    
    global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
    global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
    global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
    global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
    global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
    global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s 

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



function tinh_ap_suat_dong(deltaphi,Precess,hc1,ht1,h1)
        % ch??ng trinh tinh ap suat dong
        global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
        global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
        global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
        global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
        global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
        global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s 

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


function tinh_ap_suat_tinh(deltaphi, Precess,hc1s,ht1s,h1s)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
        % chuong trinh tinh ap suat tinh

        global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
        global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
        global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
        global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
        global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
        global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s 


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
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description

        % Can bang luu luong

        global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
        global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
        global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
        global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
        global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
        global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s 
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
        if(abs(Q1s - Q1d)<=0.001)
            Q=true;
        else
            Q=false;
        end
        
    end

function M = tinh_mo_men(phid,phi1, deltaphi)
    %myFun - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    global result m n l R Rt Rn Rlr cl cb nguy dt N W lrz ld exilon phibd  deltalanda
    global phi12 beta1 exilonn1  cbtt1  M1 capa1 deltaphi1  S1 T1 ERR1 GAP1 k1 S1s T1s ERR1s GAP1s k1s
    global phi22 beta2 exilonn2  cbtt2  M2 capa2 deltaphi2  S2 T2 ERR2 GAP2 k2 S2s T2s ERR2s GAP2s k2s
    global phi32 beta3 exilonn3  cbtt3  M3 capa3 deltaphi3  S3 T3 ERR3 GAP3 k3 S3s T3s ERR3s GAP3s k3s
    global p01 p1 a1 b1 c1 d1 e1 f1  p02 p2 a2 b2 c2 d2 e2 f2 p03 p3 a3 b3 c3 d3 e3 f3 
    global p01s p1s a1s b1s c1s d1s e1s f1s  p02s p2s a2s b2s c2s d2s e2s f2s p03s p3s a3s b3s c3s d3s e3s f3s MM
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
        if(abs(Md1-Ms1)<=0.01)
            M=true;
        else
            M=false;
        end
        
    end


%?------------------------------------------------------------------------------------------------------------------------------------------------------


%     %% TINH CHO PAD1
%     %  do day mang dau
%     for i=1:m+1
%         for j=1:n+1
%             % phan dong
%             phi1(i,j) = phi11+(i-1)*deltaphi;
%             phic1(i,j)= phi11+(i-1+1/2)*deltaphi;
%             phit1(i,j)= phi11+(i-1-1/2)*deltaphi;
%             h1(i,j) = 1+exilon*cos(phi1(i,j)-phibd)-M1*cos(phi1(i,j)-beta1)+capa1*sin(phi11-phi1(i,j));
%             hc1(i,j) = 1+exilon*cos(phic1(i,j)-phibd)-M1*cos(phic1(i,j)-beta1)+capa1*sin(phi11-phic1(i,j));
%             ht1(i,j) = 1+exilon*cos(phit1(i,j)-phibd)-M1*cos(phit1(i,j)-beta1)+capa1*sin(phi11-phit1(i,j));
%             % phan tinh
%             h1s(i,j)  = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phi1(i,j));
%             hc1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phic1(i,j));
%             ht1s(i,j) = (exilonn1/cl)*cos(phi1(i,j)-beta1)-capa1*sin(phi11-phit1(i,j));
%         end
%     end


%     % ch??ng trinh tinh ap suat dong
%     while GAP1>ERR1
%         k1 = k1+1;
%         fx1=0;
%         fy1=0;
%         for i =1:m+1
%             for j =1:n+1
%                 if i==1||i==m+1||j==1||j==n+1
%                     p1(i,j)=0;
%                 elseif i>=80&&i<=82&&j>=45&&j<=47
%                     p1(i,j)=Precess1;    
%                 else
%                     a1(i,j)= hc1(i,j)^3;
%                     b1(i,j)= ht1(i,j)^3;
%                     c1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
%                     d1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
%                     e1(i,j)= a1(i,j)+b1(i,j)+c1(i,j)+d1(i,j);
%                     f1(i,j)= 3*deltaphi*(hc1(i,j)-ht1(i,j));
%                     p1(i,j)= (a1(i,j)*p01(i+1,j)+b1(i,j)*p01(i-1,j)+c1(i,j)*p01(i,j+1)+d1(i,j)*p01(i,j-1)-f1(i,j))/e1(i,j); % he phuong trinh tinh ap suat
%                     if p1(i,j)<=0
%                         p1(i,j)=0;
%                     end
%                 end
%             end
%         end
%         % Kiem tra dieu kien hoi tu
%         for i= 1:m+1
%             for j= 1:n+1
%                 S1=S1+abs(p1(i,j)-p01(i,j));
%                 T1=T1+abs(p1(i,j));
%                 GAP1=S1/T1;
%             end
%         end
%         p01=p1;
%     end




% % chuong trinh tinh ap suat tinh
% while GAP1s>ERR1s
%             k1s = k1s+1;
%             for i =1:m+1
%                 for j =1:n+1
%                     if j==1||j==n+1||i==1||i==m+1
%                         p1s(i,j)=0;
%                         %p(i,j)=0;%p(i,j+1);
%                     elseif i>=65&&i<=95&&j>=35&&j<=55
%                     p1s(i,j)=Precess1;
%                     else
%                         a1s(i,j)= hc1s(i,j)^3;
%                     b1s(i,j)= ht1s(i,j)^3;
%                         c1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
%                         d1s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h1s(i,j)^3;
%                         e1s(i,j)= a1s(i,j)+b1s(i,j)+c1s(i,j)+d1s(i,j);
%                     f1s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
%                     p1s(i,j)= (a1s(i,j)*p01s(i+1,j)+b1s(i,j)*p01s(i-1,j)+c1s(i,j)*p01s(i,j+1)+d1s(i,j)*p01s(i,j-1)-f1s(i,j))/e1s(i,j);  % he phuong trinh tinh ap suat
                        
%                         %if p(i,j)<=0
%                                 %p(i,j)=0;
%                         %end
%                     end
%             end
%             end
%         for i= 1:m+1
%             for j= 1:n+1
%                 S1s=S1s+abs(p1s(i,j)-p01s(i,j));
%                 T1s=T1s+abs(p1s(i,j));
%                 GAP1s=S1s/T1s;
%             end
%         end
%         p01s=p1s;
%  end
 



%                 % Can bang luu luong
%                 Qd11=0;
%                 Qd12=0;
%                 Qd13=0;
%                 Qs11=0;
%                 Qs12=0;
%                 Qs13=0;
%                 for j=(n/2-1):(n/2+2)
%                     Qd11=Qd11+(1/2)*ld^2*((h1(79,44)+(1/3)*(h1(79,44))^3*((-3*p1(79,j)+4*p1(78,j)-p1(77,j)))/(2*deltaphi))*(deltalanda));
%                     Qd12=Qd12+(1/4)*ld^2*(h1(83,44)-(1/3)*(h1(83,44))^3*((3*p1(83,j)-4*p1(84,j)+p1(85,j))/(2*deltaphi))*(deltalanda));
                    
%                 end
%                 for j=44
%                     for i=79:83
%                     Qd13=Qd13-(1/6)*(h1(i,44))^3*(3*p1(i,44)-4*p1(i,43)+p1(i,42))*(deltaphi/deltalanda);
%                     end
%                 end
%                 Q1d=Qd11+Qd12+Qd13;
                
%                 for i=51
%                     for j=16:76
%                     Qs11=Qs11-(1/2)*ld^2*((1/3)*(h1s(51,16))^3*((-3*p1s(51,j)+4*p1s(50,j)-p1s(49,j))/(2*deltaphi))*(deltalanda)); 
%                     Qs12=Qs12-(1/2)*ld^2*(-(1/3)*(h1s(111,16))^3*((3*p1s(111,j)-4*p1s(112,j)+p1s(113,j))/(2*deltaphi))*(deltalanda));
%                     end
%                 end
%                 for j=16
%                     for i=51:111
%                     Qs13=Qs13+(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j-1)+p1s(i,j-2))*(deltaphi/deltalanda);
%                     %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
%                     end
%                 end
%                 Q1s=Qs11+Qs12+Qs13;










 
%  % can bang momen
%  Md1=0;
%  Msr1=0;
%  Mso11=0;
%  Mso21=0;
%  Mso31=0;
%  Mso41=0;
%  for i=1:m+1
%      for j=1:n+1
%          Md1=Md1-R*p1(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%  end
%  for i=51:111
%      for j=16:76
%          Msr1=Msr1-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%  end
%   for i=1:50
%      for j=1:n+1
%          Mso11=Mso11-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%   end
%    for i=112:m+1
%      for j=1:n+1
%          Mso21=Mso21-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%    end
%    for i=51:111
%      for j=1:15
%          Mso31=Mso31-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%    end
%     for i=51:111
%      for j=77:n+1
%          Mso41=Mso41-Rn*p1s(i,j)*sin(phi11-phi1(i,j))*deltaphi*deltalanda;
%      end
%     end
%  Ms1=Msr1+Mso11+Mso11+Mso31+Mso41;




































%  %% TINH CHO PAD2
% %  do day mang dau
% for i=1:m+1
%     for j=1:n+1
% % phan dong
% phi2(i,j) = phi21+(i-1)*deltaphi;
% phic2(i,j)= phi21+(i-1+1/2)*deltaphi;
% phit2(i,j)= phi21+(i-1-1/2)*deltaphi;
% h2(i,j) = 1+exilon*cos(phi2(i,j)-phibd)-M2*cos(phi2(i,j)-beta2)+capa2*sin(phi21-phi2(i,j));
% hc2(i,j) = 1+exilon*cos(phic2(i,j)-phibd)-M2*cos(phic2(i,j)-beta2)+capa2*sin(phi21-phic2(i,j));
% ht2(i,j) = 1+exilon*cos(phit2(i,j)-phibd)-M2*cos(phit2(i,j)-beta2)+capa2*sin(phi21-phit2(i,j));
% % phan tinh
% h2s(i,j)  = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phi2(i,j));
% hc2s(i,j) = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phic2(i,j));
% ht2s(i,j) = (exilonn2/cl)*cos(phi2(i,j)-beta2)-capa2*sin(phi21-phit2(i,j));
%     end
% end
% % ch??ng trinh tinh ap suat dong
% while GAP2>ERR2
% k2 = k2+1;
% fx1=0;
% fy1=0;
% for i =1:m+1
% for j =1:n+1
% if i==1||i==m+1||j==1||j==n+1
% p2(i,j)=0;
% elseif i>=80&&i<=82&&j>=45&&j<=47
% p2(i,j)=Precess2;    
% else
% a2(i,j)= hc2(i,j)^3;
% b2(i,j)= ht2(i,j)^3;
% c2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
% d2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
% e2(i,j)= a2(i,j)+b2(i,j)+c2(i,j)+d2(i,j);
% f2(i,j)= 3*deltaphi*(hc2(i,j)-ht2(i,j));
% p2(i,j)= (a2(i,j)*p02(i+1,j)+b2(i,j)*p02(i-1,j)+c2(i,j)*p02(i,j+1)... % he phuong trinh tinh ap suat
% +d2(i,j)*p02(i,j-1)-f2(i,j))/e2(i,j);
% if p2(i,j)<=0
% p2(i,j)=0;
% end
% end
% end
% end
% % Kiem tra dieu kien hoi tu
% for i= 1:m+1
% for j= 1:n+1
% S2=S2+abs(p2(i,j)-p02(i,j));
% T2=T2+abs(p2(i,j));
% GAP2=S2/T2;
% end
% end
% p02=p2;
% end
% % chuong trinh tinh ap suat tinh
%  while GAP2s>ERR2s
%         k2s = k2s+1;
%         for i =1:m+1
%             for j =1:n+1
%                 if j==1||j==n+1||i==1||i==m+1
%                     p2s(i,j)=0;
%                     %p(i,j)=0;%p(i,j+1);
%                 elseif i>=65&&i<=95&&j>=35&&j<=55
%                    p2s(i,j)=Precess2;
%                 else
%                     a2s(i,j)= hc2s(i,j)^3;
%                     b2s(i,j)= ht2s(i,j)^3;
%                     c2s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h2s(i,j)^3;
%                     d2s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h2s(i,j)^3;
%                     e2s(i,j)= a2s(i,j)+b2s(i,j)+c2s(i,j)+d2s(i,j);
%                    f2s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
%                   p2s(i,j)= (a2s(i,j)*p02s(i+1,j)+b2s(i,j)*p02s(i-1,j)+c2s(i,j)*p02s(i,j+1)...  % he phuong trinh tinh ap suat
%                     +d2s(i,j)*p02s(i,j-1)-f2s(i,j))/e2s(i,j);
%                 end
%            end
%         end
%     for i= 1:m+1
%         for j= 1:n+1
%             S2s=S2s+abs(p2s(i,j)-p02s(i,j));
%             T2s=T2s+abs(p2s(i,j));
%             GAP2s=S2s/T2s;
%         end
%     end
%      p02s=p2s;
%  end
 
%  Qd21=0;
%  Qd22=0;
%  Qd23=0;
%  Qs21=0;
%  Qs22=0;
%  Qs23=0;
 
%  for j=(n/2-1):(n/2+2)
%      Qd21=Qd21+(1/2)*ld^2*((h2(79,44)+(1/3)*(h2(79,44))^3*((-3*p2(79,j)+4*p2(78,j)-p2(77,j)))/(2*deltaphi))*(deltalanda));
%       Qd22=Qd22+(1/4)*ld^2*(h2(83,44)-(1/3)*(h2(83,44))^3*((3*p2(83,j)-4*p2(84,j)+p2(85,j))/(2*deltaphi))*(deltalanda));
%  end
%  for j=44
%      for i=79:83
%       Qd23=Qd23-(1/6)*(h2(i,44))^3*(3*p2(i,44)-4*p2(i,43)+p2(i,42))*(deltaphi/deltalanda);
%      end
%  end
%  Q2d=Qd21+Qd22+Qd23;
 
%  for i=51
%      for j=16:76
%      Qs21=Qs21-(1/2)*ld^2*((1/3)*(h2s(51,16))^3*((-3*p2s(51,j)+4*p2s(50,j)-p2s(49,j))/(2*deltaphi))*(deltalanda)); 
%      Qs22=Qs22-(1/2)*ld^2*(-(1/3)*(h2s(111,16))^3*((3*p2s(111,j)-4*p2s(112,j)+p2s(113,j))/(2*deltaphi))*(deltalanda));
%      end
%  end
%  for j=16
%      for i=51:111
%       Qs23=Qs23+(1/6)*((h2s(i,j))^3)*(3*p2s(i,j)-4*p2s(i,j-1)+p2s(i,j-2))*(deltaphi/deltalanda);
%       %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
%      end
%  end
%  Q2s=Qs21+Qs22+Qs23;
%  % can bang momen
%  Md2=0;
%  Msr2=0;
%  Mso12=0;
%  Mso22=0;
%  Mso32=0;
%  Mso42=0;
%  for i=1:m+1
%      for j=1:n+1
%          Md2=Md2-R*p2(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%  end
%  for i=51:111
%      for j=16:76
%          Msr2=Msr2-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%  end
%   for i=1:50
%      for j=1:n+1
%          Mso12=Mso12-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%   end
%    for i=112:m+1
%      for j=1:n+1
%          Mso22=Mso22-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%    end
%    for i=51:111
%      for j=1:15
%          Mso32=Mso32-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%    end
%     for i=51:111
%      for j=77:n+1
%          Mso42=Mso42-Rn*p2s(i,j)*sin(phi21-phi2(i,j))*deltaphi*deltalanda;
%      end
%     end
%  Ms2=Msr2+Mso12+Mso22+Mso32+Mso42;































%  %% TINH CHO PAD3
% %  do day mang dau
% for i=1:m+1
%     for j=1:n+1
% % phan dong
% phi3(i,j) = phi31+(i-1)*deltaphi;
% phic3(i,j)= phi31+(i-1+1/2)*deltaphi;
% phit3(i,j)= phi31+(i-1-1/2)*deltaphi;
% h3(i,j) = 1+exilon*cos(phi3(i,j)-phibd)-M3*cos(phi3(i,j)-beta3)+capa3*sin(phi31-phi3(i,j));
% hc3(i,j) = 1+exilon*cos(phic3(i,j)-phibd)-M3*cos(phic3(i,j)-beta3)+capa3*sin(phi31-phic3(i,j));
% ht3(i,j) = 1+exilon*cos(phit3(i,j)-phibd)-M3*cos(phit3(i,j)-beta3)+capa3*sin(phi31-phit3(i,j));
% % phan tinh
% h3s(i,j)  = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phi3(i,j));
% hc3s(i,j) = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phic3(i,j));
% ht3s(i,j) = (exilonn3/cl)*cos(phi3(i,j)-beta3)-capa3*sin(phi31-phit3(i,j));
%     end
% end
% % ch??ng trinh tinh ap suat dong
% while GAP3>ERR3
% k3 = k3+1;
% fx3=0;
% fy3=0;
% for i =1:m+1
% for j =1:n+1
% if i==1||i==m+1||j==1||j==n+1
% p3(i,j)=0;
% elseif i>=80&&i<=82&&j>=45&&j<=47
% p3(i,j)=Precess3;    
% else
% a3(i,j)= hc3(i,j)^3;
% b3(i,j)= ht3(i,j)^3;
% c3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
% d3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
% e3(i,j)= a3(i,j)+b3(i,j)+c3(i,j)+d3(i,j);
% f3(i,j)= 3*deltaphi*(hc3(i,j)-ht3(i,j));
% p3(i,j)= (a3(i,j)*p03(i+1,j)+b3(i,j)*p03(i-1,j)+c3(i,j)*p03(i,j+1)... % he phuong trinh tinh ap suat
% +d3(i,j)*p03(i,j-1)-f3(i,j))/e3(i,j);
% if p3(i,j)<=0
% p3(i,j)=0;
% end
% end
% end
% end
% % Kiem tra dieu kien hoi tu
% for i= 1:m+1
% for j= 1:n+1
% S3=S3+abs(p3(i,j)-p03(i,j));
% T3=T3+abs(p3(i,j));
% GAP3=S3/T3;
% end
% end
% p03=p3;
% end
% % chuong trinh tinh ap suat tinh
%  while GAP3s>ERR3s
%         k3s = k3s+1;
%         for i =1:m+1
%             for j =1:n+1
%                 if j==1||j==n+1||i==1||i==m+1
%                     p3s(i,j)=0;
%                     %p(i,j)=0;%p(i,j+1);
%                 elseif i>=65&&i<=95&&j>=35&&j<=55
%                    p3s(i,j)=Precess3;
%                 else
%                     a3s(i,j)= hc3s(i,j)^3;
%                     b3s(i,j)= ht3s(i,j)^3;
%                     c3s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h3s(i,j)^3;
%                     d3s(i,j)= lrz^2*(deltaphi/deltalanda)^2*h3s(i,j)^3;
%                     e3s(i,j)= a3s(i,j)+b3s(i,j)+c3s(i,j)+d3s(i,j);
%                    f3s(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
%                   p3s(i,j)= (a3s(i,j)*p03s(i+1,j)+b3s(i,j)*p03s(i-1,j)+c3s(i,j)*p03s(i,j+1)...  % he phuong trinh tinh ap suat
%                     +d3s(i,j)*p03s(i,j-1)-f3s(i,j))/e3s(i,j);
%                        %if p(i,j)<=0
%                             %p(i,j)=0;
%                        %end
%                 end
%            end
%         end
%     for i= 1:m+1
%         for j= 1:n+1
%             S3s=S3s+abs(p3s(i,j)-p03s(i,j));
%             T3s=T3s+abs(p3s(i,j));
%             GAP3s=S3s/T3s;
%         end
%     end
%      p03s=p3s;
%  end
%  Qd31=0;
%  Qd32=0;
%  Qd33=0;
%  Qs31=0;
%  Qs32=0;
%  Qs33=0;
 
%  for j=(n/2-1):(n/2+2)
%      Qd31=Qd31+(1/2)*ld^2*((h3(79,44)+(1/3)*(h3(79,44))^3*((-3*p3(79,j)+4*p3(78,j)-p3(77,j)))/(2*deltaphi))*(deltalanda));
%       Qd32=Qd32+(1/4)*ld^2*(h3(83,44)-(1/3)*(h3(83,44))^3*((3*p3(83,j)-4*p3(84,j)+p3(85,j))/(2*deltaphi))*(deltalanda));
%  end
%  for j=44
%      for i=79:83
%       Qd33=Qd33-(1/6)*(h3(i,44))^3*(3*p3(i,44)-4*p3(i,43)+p3(i,42))*(deltaphi/deltalanda);
%      end
%  end
%  Q3d=Qd31+Qd32+Qd33;
 
%  for i=51
%      for j=16:76
%      Qs31=Qs31-(1/2)*ld^2*(-(1/3)*(h3s(51,16))^3*((-3*p3s(51,j)+4*p3(50,j)-p3(49,j))/(2*deltaphi))*(deltalanda)); 
%      Qs32=Qs32-(1/2)*ld^2*(-(1/3)*(h3s(111,16))^3*((3*p3s(111,j)-4*p3s(112,j)+p3s(113,j))/(2*deltaphi))*(deltalanda));
%      end
%  end
%  for j=16
%      for i=51:111
%       Qs33=Qs33+(1/6)*((h3s(i,j))^3)*(3*p3s(i,j)-4*p3s(i,j-1)+p3s(i,j-2))*(deltaphi/deltalanda);
%       %Qs13=(1/6)*((h1s(i,j))^3)*(3*p1s(i,j)-4*p1s(i,j+1)+p1s(i,j+2))*(deltaphi/deltalanda);
%      end
%  end
%  Q3s=Qs31+Qs32+Qs33;
%  % tinh can bang momen
%   Md3=0;
%  Msr3=0;
%  Mso13=0;
%  Mso23=0;
%  Mso33=0;
%  Mso43=0;
%  for i=1:m+1
%      for j=1:n+1
%          Md3=Md3-R*p3(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%  end
%  for i=51:111
%      for j=16:76
%          Msr3=Msr3-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%  end
%   for i=1:50
%      for j=1:n+1
%          Mso13=Mso13-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%   end
%    for i=112:m+1
%      for j=1:n+1
%          Mso23=Mso23-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%    end
%    for i=51:111
%      for j=1:15
%          Mso33=Mso33-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%    end
%     for i=51:111
%      for j=77:n+1
%          Mso43=Mso43-Rn*p3s(i,j)*sin(phi31-phi3(i,j))*deltaphi*deltalanda;
%      end
%     end
%  Ms3=Msr3+Mso13+Mso23+Mso33+Mso43;
%  %% tinh luc trên các pad
%  % trên pad 1
% fx1=0;
% fy1=0;
% for i=1:m+1
% for j=1:n+1
% fx1=fx1-(p1(i,j)-p1s(i,j))*sin(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
% fy1=fy1-(p1(i,j)-p1s(i,j))*cos(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
% end
% end
% % tren pad2
% fx2=0;
% fy2=0;
% for i=1:m+1
% for j=1:n+1
% fx2=fx2-(p2(i,j)-p2s(i,j))*sin(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
% fy2=fy2-(p2(i,j)-p2s(i,j))*cos(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
% end
% end
% % trên pad 3
% fx3=0;
% fy3=0;
% for i=1:m+1
% for j=1:n+1
% fx3=fx3-(p3(i,j)-p3s(i,j))*sin(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
% fy3=fy3-(p3(i,j)-p3s(i,j))*cos(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
% end
% end
% fx=fx1+fx2+fx3;
% fy=fy1+fy2+fy3;