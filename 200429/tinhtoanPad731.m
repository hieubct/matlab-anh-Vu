% *HEADER
% ? Chương trình gửi kèm bài báo.



% 

function tinhtoanPad731()
    clc;
    format long g
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
    global m n l R Rt Rn Rlr cb nguy dt N W lrz ld exilon phibd  deltalanda
    global  cbtt1 cl capa S T ERR GAP Ss Ts ERRs GAPs
    % global p0 p a b c d e f
    % global p0s ps as bs cs ds es fs
    global QQ1 MM1 QQ2 MM2 QQ3 MM3
 
    
    % anpha1 = 2.6180e-4;%! goc nghieng cua pad1
    % anpha2 = 2.6180e-4;%! goc nghieng cua pad2
    % anpha3 = 0.015*pi/180:(0.001*pi/180):0.025*pi/180;%! goc nghieng cua pad3
 
    % Precess1=0.86:0.02:0.9;  % !�p su?t v�ng l?m1
    % Precess2=0.86:0.02:0.9;  % !�p su?t v�ng l?m1
    % Precess3= 0.4:0.02:0.6;  % !�p su?t v�ng l?m3
 
    % anpha1 = 0.005*pi/180:(0.001*pi/180):0.007*pi/180;%! goc nghieng cua pad1
    % anpha2 = 0.005*pi/180:(0.001*pi/180):0.007*pi/180;%! goc nghieng cua pad2
    % anpha3 = 0.005*pi/180:(0.001*pi/180):0.007*pi/180;%! goc nghieng cua pad3
    % Precess1=0.33:0.01:0.4;  % !�p su?t v�ng l?m1
    % Precess2=0.74:0.01:0.78;  % !�p su?t v�ng l?m1
    % Precess3= 0.34:0.01:0.6;  % !�p su?t v�ng l?m3
    
    
    
    
    anpha1 = 0.002*pi/180:(0.0005*pi/180):0.035*pi/180;%! goc nghieng cua pad1
    anpha2 = 0.002*pi/180:(0.0005*pi/180):0.035*pi/180;%! goc nghieng cua pad2
    anpha3 = 0.002*pi/180:(0.0005*pi/180):0.035*pi/180;%! goc nghieng cua pad3
    Precess1=0.01:0.005:1;  % !�p su?t v�ng l?m1
    Precess2=0.01:0.005:1;  % !�p su?t v�ng l?m1
    Precess3= 0.01:0.005:1;  % !�p su?t v�ng l?m3
 
 
    result1 = [];
    QQ1 = [];
    MM1 = [];
    result2 = [];
    QQ2 = [];
    MM2 = [];
    result3 = [];
    QQ3 = [];
    MM3 = [];
 
 
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
    exilon = 0.25;%e1/c1; % ti le lech tam
    phibd = 0.8; % goc lech ban dau cua truc
    deltalanda = 2/n; % so luoi theo phuong doc truc
    % thong so pad1
    phi11 = 3*pi/180; % goc dau cua pad
    phi12 = 117*pi/180;% goc sau cua pad
    % phi11 = 123*pi/180; % goc dau cua pad
    % phi12 = 237*pi/180;% goc sau cua pad
    beta1 = (phi11+phi12)/2; % toa do cua diem xoay
    exilonn1=0;%1.5e-5;% do nang cua pad1
    cbtt1=cb-exilonn1;
 
    deltaphi1 = (phi12-phi11)/m; %@ chia luoi theo chu vi
    S = 0;
    T = 0;
 
 
    ERR = 10e-4;
    GAP =1;
    k = 1; % he so lap
    Ss = 0;
    Ts = 0;
    ERRs = 10e-4;
    GAPs =1;
    ks = 1; % he so lap
 
 
    % thong so pad2
    phi21 = 123*pi/180; % goc dau cua pad
    phi22 = 237*pi/180;% goc sau cua pad
    beta2 = (phi22+phi21)/2; % toa do cua diem xoay
    exilonn2=0;%7e-6;% do nang cua pad2
    cbtt2=cb-exilonn2;
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
 
%     % ma tran he so ban dau cho phan hydrodynamic
%     p0 = zeros(m+1,n+1);
%     p = zeros(m+1,n+1);
%     a = zeros(m+1,n+1);
%     b = zeros(m+1,n+1);
%     c = zeros(m+1,n+1);
%     d = zeros(m+1,n+1);
%     e = zeros(m+1,n+1);
%     f = zeros(m+1,n+1);
 
 
    % % ma tran he so cho phan hydrostatic
    % p0s = zeros(m+1,n+1);
    % ps  = zeros(m+1,n+1);
    % as  = zeros(m+1,n+1);
    % bs  = zeros(m+1,n+1);
    % cs  = zeros(m+1,n+1);
    % ds  = zeros(m+1,n+1);
    % es  = zeros(m+1,n+1);
    % fs  = zeros(m+1,n+1);
 
 
 
 
 
 
 
    ss=1;
    while (ss>0.001 && phibd>0)
        fx=0;
        fy=0;
        fx1=0;
        fy1=0;
        fx2=0;
        fy2=0;
        fx3=0;
        fy3=0;
        pDong=0;
        pDong1 = zeros(m+1,n+1);
        pDong2 = zeros(m+1,n+1);
        pDong3 = zeros(m+1,n+1);
        if(fy~=0)
            phibd=phibd-atan(fx/fy);
        end
        % TINH PAD1-----------------------------------------------------------------
        M = 1-(cbtt1/cl);
        cc=1;
        Q=false;
        Momen=false;
        Pad1=false;
        while(cc<=length(anpha1) && Momen==false)
            capa=Rn*anpha1(cc)/cl;
            [phi,phic,phit,h,hc,ht,hs,hcs,hts] = tinh_do_day_mang_dau(phi11,deltaphi1,M,beta1, exilonn1);
            kk=1;
            while(kk<=length(Precess1) && Q==false )
                pDong = tinh_ap_suat_dong_one(deltaphi1,Precess1(kk),hc,ht,h);
                pTinh = tinh_ap_suat_tinh_one(deltaphi1,Precess1(kk),hcs,hts,hs);
                Q = tinh_luu_luong_one(deltaphi1,h,hs,pDong,pTinh);
                if(Q)
                    Momen = tinh_mo_men(phi11,phi,deltaphi1,pDong,pTinh,1);
                    if(Momen)
                        Pad1=true;
                        pDong1=pDong;
                        result1 = [[anpha1(cc) Precess1(kk)].'];
                    end
                end
                kk=kk+1;
            end
            cc=cc+1;
        end
        assignin('base','QQ1',QQ1);
        assignin('base','MM1',MM1);
        assignin('base','result1',result1);
        assignin('base','pDong',pDong);
        assignin('base','pTinh',pTinh);
        assignin('base','deltaphi1',deltaphi1);
        assignin('base','h',h);
        assignin('base','hs',hs);
 
 
        
        % TINH PAD2-----------------------------------------------------------------
        M = 1-(cbtt2/cl);
        cc=1;
        Q=false;
        Momen=false;
        Pad2=false;
        while(cc<=length(anpha2) && Momen==false && Pad1==true)
            capa=Rn*anpha2(cc)/cl;
            [phi,phic,phit,h,hc,ht,hs,hcs,hts] = tinh_do_day_mang_dau(phi21,deltaphi2,M,beta2, exilonn2);
            kk=1;
            while(kk<=length(Precess2) && Q==false)
                pDong = tinh_ap_suat_dong(deltaphi2,Precess2(kk),hc,ht,h);
                pTinh = tinh_ap_suat_tinh(deltaphi2,Precess2(kk),hcs,hts,hs);
                Q = tinh_luu_luong(deltaphi2,h,hs,pDong,pTinh);
                if(Q)
                    % anpha2(cc)
                    Momen = tinh_mo_men(phi21,phi,deltaphi2,pDong,pTinh,2);
                    if(Momen)
                        Pad2=true;
                        pDong2=pDong;
                        result2 = [[anpha2(cc) Precess2(kk)].'];
                    end
                end
                kk=kk+1;
            end
            cc=cc+1;
        end
 
        QQ2
        MM2
        result2
 
 
        % TINH PAD 3-----------------------------------------------------------------
        M = 1-(cbtt3/cl);
        cc=1;
        Q=false;
        Momen=false;
        Pad3 = false;
        while(cc<=length(anpha3) && Momen==false && Pad1==true && Pad2==true)
            capa=Rn*anpha3(cc)/cl;
            [phi,phic,phit,h,hc,ht,hs,hcs,hts] = tinh_do_day_mang_dau(phi31,deltaphi3,M,beta3, exilonn3);
            kk=1;
            while(kk<=length(Precess3) && Q==false)
                pDong = tinh_ap_suat_dong_three(deltaphi3,Precess3(kk),hc,ht,h);
                pTinh = tinh_ap_suat_tinh_three(deltaphi3,Precess3(kk),hcs,hts,hs);
                Q = tinh_luu_luong_three(deltaphi3,h,hs,pDong,pTinh);
                if(Q)
                    %  anpha3(cc)
                    Momen = tinh_mo_men(phi31,phi,deltaphi3,pDong,pTinh,3);
                    if(Momen)
                        Pad3=true;
                        pDong3=pDong;
                        result3 = [[anpha3(cc) Precess3(kk)].'];
                    end
                end
                kk=kk+1;
            end
            cc=cc+1;
        end
        QQ3
        MM3
        result3
        if(Pad1==true && Pad2==true && Pad3==true)
            % %t�nh fx, fy pad1
            % [fxp,fyp]=tinh_luc(m,n,pDong1,phi11,deltaphi1,deltalanda);
            % fx=fxp;
            % fy=fyp;
            % %t�nh fx, fy pad2
            % [fxp,fyp]=tinh_luc(m,n,pDong2,phi21,deltaphi2,deltalanda);
            % fx=fx+fxp;
            % fy=fy+fyp;
            % %t�nh fx, fy pad3
            % [fxp,fyp]=tinh_luc(m,n,pDong3,phi31,deltaphi3,deltalanda);
            % fx=fx+fxp;
            % fy=fy+fyp;
            %[fx,fy]=tinh_luc(m,n,pDong1,phi11,deltaphi1,deltalanda)+tinh_luc(m,n,pDong2,phi21,deltaphi2,deltalanda)+tinh_luc(m,n,pDong3,phi31,deltaphi3,deltalanda);
            [fx1,fy1]=tinh_luc(m,n,pDong1,phi11,deltaphi1,deltalanda);
            [fx2,fy2]=tinh_luc(m,n,pDong2,phi21,deltaphi2,deltalanda);
            [fx3,fy3]=tinh_luc(m,n,pDong3,phi31,deltaphi3,deltalanda);
            fx=fx1+fx2+fx3;
            fy=fy1+fy2+fy3;
            ss=abs(fx/fy);
        else
            break;
        end
    end
    phibd
end
 
 
function [phi,phic,phit,h,hc,ht,hs,hcs,hts] = tinh_do_day_mang_dau(phid,deltaphi,M,beta,exilonn)
    %  do day mang dau
    global  m n cl exilon phibd capa
    for i=1:m+1
        for j=1:n+1
            % phan dong
            phi(i,j) = phid+(i-1)*deltaphi;
            phic(i,j)= phid+(i-1+1/2)*deltaphi;
            phit(i,j)= phid+(i-1-1/2)*deltaphi;
            h(i,j) = 1+exilon*cos(phi(i,j)-phibd)-M*cos(phi(i,j)-beta)+capa*sin(phid-phi(i,j));
            hc(i,j) = 1+exilon*cos(phic(i,j)-phibd)-M*cos(phic(i,j)-beta)+capa*sin(phid-phic(i,j));
            ht(i,j) = 1+exilon*cos(phit(i,j)-phibd)-M*cos(phit(i,j)-beta)+capa*sin(phid-phit(i,j));
            % phan tinh
            hs(i,j)  = (exilonn/cl)*cos(phi(i,j)-beta)-capa*sin(phid-phi(i,j));
            hcs(i,j) = (exilonn/cl)*cos(phi(i,j)-beta)-capa*sin(phid-phic(i,j));
            hts(i,j) = (exilonn/cl)*cos(phi(i,j)-beta)-capa*sin(phid-phit(i,j));
        end
    end
end
 
% ap suat dong cho  pad1
function pp1 = tinh_ap_suat_dong_one(deltaphi,Precess,hc,ht,h)
    % ch??ng trinh tinh ap suat dong
    global  m n ld  deltalanda ERR
    % global p0 p a b c d e f
    % ma tran he so ban dau cho phan hydrodynamic
    p0 = zeros(m+1,n+1);
    pp1 = zeros(m+1,n+1);
    a = zeros(m+1,n+1);
    b = zeros(m+1,n+1);
    c = zeros(m+1,n+1);
    d = zeros(m+1,n+1);
    e = zeros(m+1,n+1);
    f = zeros(m+1,n+1);
    GAP=1;
    k=1;
    S = 0;
    T = 0;
    while GAP>ERR
        k = k+1;
        fx=0;
        fy=0;
        for i =1:m+1
            for j =1:n+1
                if i==1||i==m+1||j==1||j==n+1
                    pp1(i,j)=0;
                elseif i>=55&&i<=59&&j>=45&&j<=47
                    pp1(i,j)=Precess;
                else
                    a(i,j)= hc(i,j)^3;
                    b(i,j)= ht(i,j)^3;
                    c(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    d(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                    f(i,j)= 3*deltaphi*(hc(i,j)-ht(i,j));
                    pp1(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)+d(i,j)*p0(i,j-1)-f(i,j))/e(i,j); % he phuong trinh tinh ap suat
                    if pp1(i,j)<=0
                        pp1(i,j)=0;
                    end
                end
            end
        end
        % Kiem tra dieu kien hoi tu
        for i= 1:m+1
            for j= 1:n+1
                S=S+abs(pp1(i,j)-p0(i,j));
                T=T+abs(pp1(i,j));
                GAP=S/T;
            end
        end
        p0=pp1;
    end
end
 
%ap suat dong cho pad2
 
function p = tinh_ap_suat_dong(deltaphi,Precess,hc,ht,h)
    % ch??ng trinh tinh ap suat dong
    global  m n ld  deltalanda ERR
    % global p0 p a b c d e f
    % ma tran he so ban dau cho phan hydrodynamic
    p0 = zeros(m+1,n+1);
    p = zeros(m+1,n+1);
    a = zeros(m+1,n+1);
    b = zeros(m+1,n+1);
    c = zeros(m+1,n+1);
    d = zeros(m+1,n+1);
    e = zeros(m+1,n+1);
    f = zeros(m+1,n+1);
    GAP=1;
    k=1;
    S = 0;
    T = 0;
    while GAP>ERR
        k = k+1;
        fx=0;
        fy=0;
        for i =1:m+1
            for j =1:n+1
                if i==1||i==m+1||j==1||j==n+1
                    p(i,j)=0;
                elseif i>=79&&i<=83&&j>=45&&j<=47
                    p(i,j)=Precess;
                else
                    a(i,j)= hc(i,j)^3;
                    b(i,j)= ht(i,j)^3;
                    c(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    d(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                    f(i,j)= 3*deltaphi*(hc(i,j)-ht(i,j));
                    p(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)+d(i,j)*p0(i,j-1)-f(i,j))/e(i,j); % he phuong trinh tinh ap suat
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
end
 
 
%ap suat dong cho pad3
function pp3 = tinh_ap_suat_dong_three(deltaphi,Precess,hc,ht,h)
    % ch??ng trinh tinh ap suat dong
    global  m n ld  deltalanda ERR
    % global p0 p a b c d e f
    % ma tran he so ban dau cho phan hydrodynamic
    p0 = zeros(m+1,n+1);
    pp3 = zeros(m+1,n+1);
    a = zeros(m+1,n+1);
    b = zeros(m+1,n+1);
    c = zeros(m+1,n+1);
    d = zeros(m+1,n+1);
    e = zeros(m+1,n+1);
    f = zeros(m+1,n+1);
    GAP=1;
    k=1;
    S = 0;
    T = 0;
    while GAP>ERR
        k = k+1;
        fx=0;
        fy=0;
        for i =1:m+1
            for j =1:n+1
                if i==1||i==m+1||j==1||j==n+1
                    pp3(i,j)=0;
                elseif i>=37&&i<=41&&j>=45&&j<=47
                    pp3(i,j)=Precess;
                else
                    a(i,j)= hc(i,j)^3;
                    b(i,j)= ht(i,j)^3;
                    c(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    d(i,j)= ld^2*(deltaphi/deltalanda)^2*h(i,j)^3;
                    e(i,j)= a(i,j)+b(i,j)+c(i,j)+d(i,j);
                    f(i,j)= 3*deltaphi*(hc(i,j)-ht(i,j));
                    pp3(i,j)= (a(i,j)*p0(i+1,j)+b(i,j)*p0(i-1,j)+c(i,j)*p0(i,j+1)+d(i,j)*p0(i,j-1)-f(i,j))/e(i,j); % he phuong trinh tinh ap suat
                    if pp3(i,j)<=0
                        pp3(i,j)=0;
                    end
                end
            end
        end
        % Kiem tra dieu kien hoi tu
        for i= 1:m+1
            for j= 1:n+1
                S=S+abs(pp3(i,j)-p0(i,j));
                T=T+abs(pp3(i,j));
                GAP=S/T;
            end
        end
        p0=pp3;
    end
end
 
% ap suat tinh cho pad1
 
function ps1 = tinh_ap_suat_tinh_one(deltaphi, Precess,hcs,hts,hs)
    % chuong trinh tinh ap suat tinh
    global  m n lrz deltalanda Ss Ts ERRs
    % global p0s ps as bs cs ds es fs
    % ma tran he so cho phan hydrostatic
    p0s = zeros(m+1,n+1);
    ps1  = zeros(m+1,n+1);
    as  = zeros(m+1,n+1);
    bs  = zeros(m+1,n+1);
    cs  = zeros(m+1,n+1);
    ds  = zeros(m+1,n+1);
    es  = zeros(m+1,n+1);
    fs  = zeros(m+1,n+1);
    GAPs=1;
    ks=1;
    Ss = 0;
    Ts = 0;
    while GAPs>ERRs
        ks = ks+1;
        for i =1:m+1
            for j =1:n+1
                if j==1||j==n+1||i==1||i==m+1
                    ps1(i,j)=0;
                    %p(i,j)=0;%p(i,j+1);
                elseif i>=53&&i<=109&&j>=23&&j<=69
                    ps1(i,j)=Precess;
                else
                    as(i,j)= hcs(i,j)^3;
                    bs(i,j)= hts(i,j)^3;
                    cs(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    ds(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    es(i,j)= as(i,j)+bs(i,j)+cs(i,j)+ds(i,j);
                    fs(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                    ps1(i,j)= (as(i,j)*p0s(i+1,j)+bs(i,j)*p0s(i-1,j)+cs(i,j)*p0s(i,j+1)+ds(i,j)*p0s(i,j-1)-fs(i,j))/es(i,j);  % he phuong trinh tinh ap suat
                    %if p(i,j)<=0
                    %p(i,j)=0;
                    %end
                end
            end
        end
        for i= 1:m+1
            for j= 1:n+1
                Ss=Ss+abs(ps1(i,j)-p0s(i,j));
                Ts=Ts+abs(ps1(i,j));
                GAPs=Ss/Ts;
            end
        end
        p0s=ps1;
    end
end
 
 
 
 
 
function ps = tinh_ap_suat_tinh(deltaphi, Precess,hcs,hts,hs)
    % chuong trinh tinh ap suat tinh
    global  m n lrz deltalanda Ss Ts ERRs
    % global p0s ps as bs cs ds es fs
    % ma tran he so cho phan hydrostatic
    p0s = zeros(m+1,n+1);
    ps  = zeros(m+1,n+1);
    as  = zeros(m+1,n+1);
    bs  = zeros(m+1,n+1);
    cs  = zeros(m+1,n+1);
    ds  = zeros(m+1,n+1);
    es  = zeros(m+1,n+1);
    fs  = zeros(m+1,n+1);
    GAPs=1;
    ks=1;
    Ss = 0;
    Ts = 0;
    while GAPs>ERRs
        ks = ks+1;
        for i =1:m+1
            for j =1:n+1
                if j==1||j==n+1||i==1||i==m+1
                    ps(i,j)=0;
                    %p(i,j)=0;%p(i,j+1);
                elseif i>=53&&i<=109&&j>=23&&j<=69
                    ps(i,j)=Precess;
                else
                    as(i,j)= hcs(i,j)^3;
                    bs(i,j)= hts(i,j)^3;
                    cs(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    ds(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    es(i,j)= as(i,j)+bs(i,j)+cs(i,j)+ds(i,j);
                    fs(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                    ps(i,j)= (as(i,j)*p0s(i+1,j)+bs(i,j)*p0s(i-1,j)+cs(i,j)*p0s(i,j+1)+ds(i,j)*p0s(i,j-1)-fs(i,j))/es(i,j);  % he phuong trinh tinh ap suat
                    %if p(i,j)<=0
                    %p(i,j)=0;
                    %end
                end
            end
        end
        for i= 1:m+1
            for j= 1:n+1
                Ss=Ss+abs(ps(i,j)-p0s(i,j));
                Ts=Ts+abs(ps(i,j));
                GAPs=Ss/Ts;
            end
        end
        p0s=ps;
    end
end
 
 
function ps3 = tinh_ap_suat_tinh_three(deltaphi, Precess,hcs,hts,hs)
    % chuong trinh tinh ap suat tinh
    global  m n lrz deltalanda Ss Ts ERRs
    % global p0s ps as bs cs ds es fs
    % ma tran he so cho phan hydrostatic
    p0s = zeros(m+1,n+1);
    ps3  = zeros(m+1,n+1);
    as  = zeros(m+1,n+1);
    bs  = zeros(m+1,n+1);
    cs  = zeros(m+1,n+1);
    ds  = zeros(m+1,n+1);
    es  = zeros(m+1,n+1);
    fs  = zeros(m+1,n+1);
    GAPs=1;
    ks=1;
    Ss = 0;
    Ts = 0;
    while GAPs>ERRs
        ks = ks+1;
        for i =1:m+1
            for j =1:n+1
                if j==1||j==n+1||i==1||i==m+1
                    ps3(i,j)=0;
                    %p(i,j)=0;%p(i,j+1);
                elseif i>=53&&i<=109&&j>=23&&j<=69
                    ps3(i,j)=Precess;
                else
                    as(i,j)= hcs(i,j)^3;
                    bs(i,j)= hts(i,j)^3;
                    cs(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    ds(i,j)= lrz^2*(deltaphi/deltalanda)^2*hs(i,j)^3;
                    es(i,j)= as(i,j)+bs(i,j)+cs(i,j)+ds(i,j);
                    fs(i,j)=0;%3*deltaphi1*(hc(i,j)-ht(i,j));
                    ps3(i,j)= (as(i,j)*p0s(i+1,j)+bs(i,j)*p0s(i-1,j)+cs(i,j)*p0s(i,j+1)+ds(i,j)*p0s(i,j-1)-fs(i,j))/es(i,j);  % he phuong trinh tinh ap suat
                    %if p(i,j)<=0
                    %p(i,j)=0;
                    %end
                end
            end
        end
        for i= 1:m+1
            for j= 1:n+1
                Ss=Ss+abs(ps3(i,j)-p0s(i,j));
                Ts=Ts+abs(ps3(i,j));
                GAPs=Ss/Ts;
            end
        end
        p0s=ps3;
    end
end
 
% can bawng luu luowng pad2
 
function Q = tinh_luu_luong(deltaphi,h,hs,p,ps)
    % Can bang luu luong
    global  n  ld deltalanda
    global QQ2 
 
    Qd1=0;
    Qd2=0;
    Qd3=0;
    Qs1=0;
    Qs2=0;
    Qs3=0;
    for j=(round(n/2-1):(n/2+2))
        Qd1=Qd1+(1/2)*ld^2*((h(79,44)-(1/3)*(h(79,44))^3*((-3*p(77,j)+4*p(78,j)-p(79,j)))/(2*deltaphi))*(deltalanda));
        Qd2=Qd2+(1/2)*ld^2*((-h(83,44)-(1/3)*(h(83,44))^3*((3*p(85,j)-4*p(84,j)+p(83,j)))/(2*deltaphi))*(deltalanda));
 
    end
    for j=44
        for i=79:83
            Qd3=Qd3-(1/6)*(h(i,44))^3*(3*p(i,44)-4*p(i,43)+p(i,42))*(deltaphi/deltalanda);
        end
    end
    Qd=Qd1+Qd2+Qd3;
    for i=53
        for j=23:69
            Qs1=Qs1-(1/2)*ld^2*((1/3)*(hs(53,23))^3*((-3*ps(53,j)+4*ps(52,j)-ps(51,j))/(2*deltaphi))*(deltalanda));
            Qs2=Qs2-(1/2)*ld^2*(-(1/3)*(hs(109,23))^3*((3*ps(109,j)-4*ps(110,j)+ps(111,j))/(2*deltaphi))*(deltalanda));
        end
    end
    for j=23
        for i=53:109
            Qs3=Qs3+(1/6)*((hs(i,j))^3)*(3*ps(i,j)-4*ps(i,j-1)+ps(i,j-2))*(deltaphi/deltalanda);
        end
    end
    Qs=Qs1+Qs2+Qs3;
 
    
   
        %QQ1 = [QQ1 [Qd Qs].'];
    %elseif(index==2)
       QQ2 = [QQ2 [Qd Qs].'];
    %elseif(index==3)
        %QQ3 = [QQ3 [Qd Qs].'];
    %end
 
    
    % KIEM TRA CAN BANG LUU LUONG.
    if(abs(Qd-Qs)<=0.0001)
        Q=true;
    else
        Q=false;
    end
end
 
% can bawng luu luong pad1
 
function Q1 = tinh_luu_luong_one(deltaphi,h,hs,p,ps)
    % Can bang luu luong
    global  n  ld deltalanda
    global QQ1
 
    Qd11=0;
    Qd12=0;
    Qd13=0;
    Qs11=0;
    Qs12=0;
    Qs13=0;
    for j=(round(n/2-1):(n/2+2))
        Qd11=Qd11+(1/2)*ld^2*((h(55,44)-(1/3)*(h(55,44))^3*((-3*p(54,j)+4*p(55,j)-p(56,j)))/(2*deltaphi))*(deltalanda));
        Qd12=Qd12+(1/2)*ld^2*((-h(59,44)-(1/3)*(h(59,44))^3*((3*p(60,j)-4*p(59,j)+p(58,j)))/(2*deltaphi))*(deltalanda));
         
    end
    for j=44
        for i=55:59
            Qd13=Qd13-(1/6)*(h(i,44))^3*(3*p(i,44)-4*p(i,43)+p(i,42))*(deltaphi/deltalanda);
        end
    end
    Qd1=Qd11+Qd12+Qd13;
    for i=53
        for j=23:69
            Qs11=Qs11-(1/2)*ld^2*((1/3)*(hs(53,23))^3*((-3*ps(53,j)+4*ps(52,j)-ps(51,j))/(2*deltaphi))*(deltalanda));
            Qs12=Qs12-(1/2)*ld^2*(-(1/3)*(hs(109,23))^3*((3*ps(109,j)-4*ps(110,j)+ps(111,j))/(2*deltaphi))*(deltalanda));
        end
    end
    for j=23
        for i=53:109
            Qs13=Qs13+(1/6)*((hs(i,j))^3)*(3*ps(i,j)-4*ps(i,j-1)+ps(i,j-2))*(deltaphi/deltalanda);
        end
    end
    Qs1=Qs11+Qs12+Qs13;
 
    QQ1 = [QQ1 [Qd1 Qs1].'];
        
    % KIEM TRA CAN BANG LUU LUONG.
    if(abs(Qd1-Qs1)<=0.0001)
        Q1=true;
    else
        Q1=false;
    end
end
 
% can bang luu luong pad3
 
function Q3 = tinh_luu_luong_three(deltaphi,h,hs,p,ps)
    % Can bang luu luong
    global  n  ld deltalanda
    global QQ3
 
    Qd31=0;
    Qd32=0;
    Qd33=0;
    Qs31=0;
    Qs32=0;
    Qs33=0;
    for j=(round(n/2-1):(n/2+2))
        Qd31=Qd31+(1/2)*ld^2*((h(37,44)-(1/3)*(h(37,44))^3*((-3*p(36,j)+4*p(37,j)-p(38,j)))/(2*deltaphi))*(deltalanda));
        Qd32=Qd32+(1/2)*ld^2*((-h(41,44)-(1/3)*(h(41,44))^3*((3*p(42,j)-4*p(41,j)+p(40,j)))/(2*deltaphi))*(deltalanda));
        
    end
    for j=44
        for i=37:41
            Qd33=Qd33-(1/6)*(h(i,44))^3*(3*p(i,44)-4*p(i,43)+p(i,42))*(deltaphi/deltalanda);
        end
    end
    Qd3=Qd31+Qd32+Qd33;
    for i=53
        for j=23:69
            Qs31=Qs31-(1/2)*ld^2*((1/3)*(hs(53,23))^3*((-3*ps(53,j)+4*ps(52,j)-ps(51,j))/(2*deltaphi))*(deltalanda));
            Qs32=Qs32-(1/2)*ld^2*(-(1/3)*(hs(109,23))^3*((3*ps(109,j)-4*ps(110,j)+ps1(111,j))/(2*deltaphi))*(deltalanda));
        end
    end
    for j=23
        for i=53:109
            Qs33=Qs33+(1/6)*((hs(i,j))^3)*(3*ps(i,j)-4*ps(i,j-1)+ps(i,j-2))*(deltaphi/deltalanda);
        end
    end
    Qs3=Qs31+Qs32+Qs33;
 
    QQ3 = [QQ3 [Qd3 Qs3].'];
        
    % KIEM TRA CAN BANG LUU LUONG.
    if(abs(Qd3-Qs3)<=0.0001)
        Q3=true;
    else
        Q3=false;
    end
end
 
 
function Momen = tinh_mo_men(phid,phi,deltaphi,p,ps,index)
    global  m n R Rn  deltalanda MM1 MM2 MM3
    % can bang momen
    Md=0;
    Msr=0;
    Mso1=0;
    Mso2=0;
    Mso3=0;
    Mso4=0;
 
    for i=1:m+1
        for j=1:n+1
            Md=Md-R*p(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    for i=53:109
        for j=23:69
            Msr=Msr-Rn*ps(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    for i=1:52
        for j=1:n+1
            Mso1=Mso1-Rn*ps(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    for i=110:m+1
        for j=1:n+1
            Mso2=Mso2-Rn*ps(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    for i=53:109
        for j=1:22
            Mso3=Mso3-Rn*ps(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    for i=53:109
        for j=69:n+1
            Mso4=Mso4-Rn*ps(i,j)*sin(phid-phi(i,j))*deltaphi*deltalanda;
        end
    end
    Ms=Msr+Mso1+Mso2+Mso3+Mso4;
 
    if(index==1)
        MM1 = [MM1 [Md Ms].'];
    elseif(index==2)
        MM2 = [MM2 [Md Ms].'];
    elseif(index==3)
        MM3 = [MM3 [Md Ms].'];
    end
 
    
    if(abs(Md-Ms)<=0.05)
        Momen=true;
    else
        Momen=false;
    end
end
 
 
 
function [fxp,fyp] = tinh_luc(m,n,p, phid, deltaphi, deltalanda)
    fxp=0;
    fyp=0;
    for i=1:m+1
        for j=1:n+1
            fxp=fxp-p(i,j)*sin(phid+(i-1)*deltaphi)*deltaphi*deltalanda;
            fyp=fyp-p(i,j)*cos(phid+(i-1)*deltaphi)*deltaphi*deltalanda;
        end
    end
end


