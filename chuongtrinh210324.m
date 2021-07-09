clc
clear
clearvars;
exilon=0.05:0.1:0.45;
phi_all_1=zeros(161,91,15);
phi_all_2=zeros(161,91,15);
phi_all_3=zeros(161,91,15);
p_all_1=zeros(161,91,15);
p_all_2=zeros(161,91,15);
p_all_3=zeros(161,91,15);
h_all_1=zeros(161,91,15);
h_all_2=zeros(161,91,15);
h_all_3=zeros(161,91,15);



for count=1:size(exilon,2)
    strindex = strcat(int2str(count),' exilon = ',num2str(exilon(count)));
    R = 100e-3; % ban kinh truc
    l = 150e-3; % chieu dai o
    Rt=100.2e-3;% ban kinh trong cua pad
    Rn=115e-3; % ban kinh ngoai cua pad
    Rlr=100.15e-3; % ban kinh lap rap
    c1 = Rt-R; % khe ho cua o truc ( tam truc va tam pad trung nhau)
    cb = Rlr-R; % khe ho cai dat ( tam truc va tam o trung nhau)
    nguy = 0.063; % do nhot cua dau
    d1 = 2*R; % duong kinh truc
    %posi = c1/R; % ti le khe ho
    M =0.5;% 1 -(cb/c1);
    ld = l/d1; % ti le chieu dai tren duong kinh
    m = 160; % so luoi theo chu vi
    n = 90; % so luoi theo phuong doc truc
    %e1 =0.75e-04;% 0.1e-3; % do lech tam ban dau
    %exilon = 0.35;%0.05;%e1/c1; % ti le lech tam
    %for v=1:length(exilon)
    phibd = 0.8887; % goc lech ban dau cua truc
    
    % tinh cho pad dau tien
    phi11 = 3*pi/180; % goc dau cua pad
    phi12 = 117*pi/180;% goc sau cua pad
    beta1 = (phi11+phi12)/2; % toa do cua diem xoay
    % tinh cho pad dau thu 2
    phi21 = 123*pi/180; % goc dau cua pad
    phi22 = 237*pi/180;% goc sau cua pad
    beta2 = (phi22+phi21)/2; % toa do cua diem xoay
    % tinh cho pad dau thu 3
    phi31 = 243*pi/180; % goc dau cua pad
    phi32 = 357*pi/180;% goc sau cua pad
    beta3 = (phi32+phi31)/2; % toa do cua diem xoay
    
    % tinh cac thong so ban dau
    deltalanda = 1/n; % so luoi theo phuong doc truc
    deltaphi = (phi22-phi21)/m; % chia luoi theo chu vi
    % khai bao cac ma tran he so ban dau
    fx=0.01;
    fy=1;
    ss=abs(fx/fy);
    %%tinh goc trang thai
    while ss>0.001
        phibd=phibd-atan(fx/fy);
        % tinh cho pad dau tien
        fx=0;
        fy=0;
        phi11 = 3*pi/180; % goc dau cua pad
        phi12 = 117*pi/180;% goc sau cua pad
        beta1 = (phi11+phi12)/2; % toa do cua diem xoay
        % tinh cho pad dau thu 2
        phi21 = 123*pi/180; % goc dau cua pad
        phi22 = 237*pi/180;% goc sau cua pad
        beta2 = (phi22+phi21)/2; % toa do cua diem xoay
        % tinh cho pad dau thu 3
        phi31 = 243*pi/180; % goc dau cua pad
        phi32 = 357*pi/180;% goc sau cua pad
        beta3 = (phi32+phi31)/2; % toa do cua diem xoay
        
        % tinh cac thong so ban dau
        deltalanda = 1/n; % so luoi theo phuong doc truc
        deltaphi = (phi22-phi21)/m; % chia luoi theo chu vi
        % khai bao cac ma tran he so ban dau
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
        %h = zeros(m+1,n+1);
        % tinh do day mang dau
        
        for i=1:m+1
            for j=1:n+1
                % ma pad 1
                phi1(i,j) = phi11+(i-1)*deltaphi;
                phic1(i,j)= phi11+(i-1+1/2)*deltaphi;
                phit1(i,j)= phi11+(i-1-1/2)*deltaphi;
                h1(i,j) = 1+exilon(count)*cos(phi1(i,j)-phibd)-M*cos(phi1(i,j)-beta1);
                hc1(i,j) = 1+exilon(count)*cos(phic1(i,j)-phibd)-M*cos(phic1(i,j)-beta1);
                ht1(i,j) = 1+exilon(count)*cos(phit1(i,j)-phibd)-M*cos(phit1(i,j)-beta1);
                % ma pad 2
                phi2(i,j) = phi21+(i-1)*deltaphi;
                phic2(i,j)= phi21+(i-1+1/2)*deltaphi;
                phit2(i,j)= phi21+(i-1-1/2)*deltaphi;
                h2(i,j) = 1+exilon(count)*cos(phi2(i,j)-phibd)-M*cos(phi2(i,j)-beta2);
                hc2(i,j) = 1+exilon(count)*cos(phic2(i,j)-phibd)-M*cos(phic2(i,j)-beta2);
                ht2(i,j) = 1+exilon(count)*cos(phit2(i,j)-phibd)-M*cos(phit2(i,j)-beta2);
                % ma pad 3
                phi3(i,j) = phi31+(i-1)*deltaphi;
                phic3(i,j)= phi31+(i-1+1/2)*deltaphi;
                phit3(i,j)= phi31+(i-1-1/2)*deltaphi;
                h3(i,j) = 1+exilon(count)*cos(phi3(i,j)-phibd)-M*cos(phi3(i,j)-beta3);
                hc3(i,j) = 1+exilon(count)*cos(phic3(i,j)-phibd)-M*cos(phic3(i,j)-beta3);
                ht3(i,j) = 1+exilon(count)*cos(phit3(i,j)-phibd)-M*cos(phit3(i,j)-beta3);
            end
        end
        % chuong trinh tinh ap suat
        % khai bao cac bien ban dau
        S1 = 0;
        T1 = 0;
        ERR1 = 10e-4;
        GAP1 =1;
        k1 = 1; % he so lap
        
        S2 = 0;
        T2 = 0;
        ERR2 = 10e-4;
        GAP2 =1;
        k2 = 1; % he so lap
        
        S3 = 0;
        T3= 0;
        ERR3 = 10e-4;
        GAP3 =1;
        k3 = 1; % he so lap
        while GAP1>ERR1
            k1 = k1+1;
            fx1=0;
            fy1=0;
            for i =1:m+1
                for j =1:n+1
                    if i==1||i==m+1||j==1||j==n+1
                        p1(i,j)=0;
                    else
                        a1(i,j)= hc1(i,j)^3;
                        b1(i,j)= ht1(i,j)^3;
                        c1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                        d1(i,j)= ld^2*(deltaphi/deltalanda)^2*h1(i,j)^3;
                        e1(i,j)= a1(i,j)+b1(i,j)+c1(i,j)+d1(i,j);
                        f1(i,j)= 3*deltaphi*(hc1(i,j)-ht1(i,j));
                        p1(i,j)= (a1(i,j)*p01(i+1,j)+b1(i,j)*p01(i-1,j)+c1(i,j)*p01(i,j+1)... % he phuong trinh tinh ap suat
                            +d1(i,j)*p01(i,j-1)-f1(i,j))/e1(i,j);
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
        for i=1:m+1
            for j=1:n+1
                fx1=fx1-p1(i,j)*sin(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
                fy1=fy1-p1(i,j)*cos(phi11+(i-1)*deltaphi)*deltaphi*deltalanda;
            end
        end
        % tinh cho pad 2
        while GAP2>ERR2
            k2 = k2+1;
            fx2=0;
            fy2=0;
            for i =1:m+1
                for j =1:n+1
                    if i==1||i==m+1||j==1||j==n+1
                        p2(i,j)=0;
                    else
                        a2(i,j)= hc2(i,j)^3;
                        b2(i,j)= ht2(i,j)^3;
                        c2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
                        d2(i,j)= ld^2*(deltaphi/deltalanda)^2*h2(i,j)^3;
                        e2(i,j)= a2(i,j)+b2(i,j)+c2(i,j)+d2(i,j);
                        f2(i,j)= 3*deltaphi*(hc2(i,j)-ht2(i,j));
                        p2(i,j)= (a2(i,j)*p02(i+1,j)+b2(i,j)*p02(i-1,j)+c2(i,j)*p02(i,j+1)... % he phuong trinh tinh ap suat
                            +d2(i,j)*p02(i,j-1)-f2(i,j))/e2(i,j);
                        if p2(i,j)<=0
                            p2(i,j)=0;
                        end
                    end
                end
            end
            % Kiem tra dieu kien hoi tu
            for i= 1:m+1
                for j= 1:n+1
                    S2=S2+abs(p2(i,j)-p02(i,j));
                    T2=T2+abs(p2(i,j));
                    GAP2=S2/T2;
                end
            end
            p02=p2;
        end
        for i=1:m+1
            for j=1:n+1
                fx2=fx2-p2(i,j)*sin(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
                fy2=fy2-p2(i,j)*cos(phi21+(i-1)*deltaphi)*deltaphi*deltalanda;
            end
        end
        % tinh cho pad 2
        while GAP3>ERR3
            k3 = k3+1;
            fx3=0;
            fy3=0;
            for i =1:m+1
                for j =1:n+1
                    if i==1||i==m+1||j==1||j==n+1
                        p3(i,j)=0;
                    else
                        a3(i,j)= hc3(i,j)^3;
                        b3(i,j)= ht3(i,j)^3;
                        c3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
                        d3(i,j)= ld^2*(deltaphi/deltalanda)^2*h3(i,j)^3;
                        e3(i,j)= a3(i,j)+b3(i,j)+c3(i,j)+d3(i,j);
                        f3(i,j)= 3*deltaphi*(hc3(i,j)-ht3(i,j));
                        p3(i,j)= (a3(i,j)*p03(i+1,j)+b3(i,j)*p03(i-1,j)+c3(i,j)*p03(i,j+1)... % he phuong trinh tinh ap suat
                            +d3(i,j)*p03(i,j-1)-f3(i,j))/e3(i,j);
                        if p3(i,j)<=0
                            p3(i,j)=0;
                        end
                    end
                end
            end
            % Kiem tra dieu kien hoi tu
            for i= 1:m+1
                for j= 1:n+1
                    S3=S3+abs(p3(i,j)-p03(i,j));
                    T3=T3+abs(p3(i,j));
                    GAP3=S3/T3;
                end
            end
            p03=p3;
        end
        for i=1:m+1
            for j=1:n+1
                fx3=fx3-p3(i,j)*sin(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
                fy3=fy3-p3(i,j)*cos(phi31+(i-1)*deltaphi)*deltaphi*deltalanda;
            end
        end
        fx=fx1+fx2+fx3;
        fy=fy1+fy2+fy3;
        ss=abs(fx/fy);
    end
    p11=transpose(p1);
    p21=transpose(p2);
    p31=transpose(p3);
    [x1,y1]=max(p11);
    [x2,y2]=max(p21);
    [x3,y3]=max(p31);
    [cx11,cot1]=max(x1);
    [cx12,cot2]=max(x2);
    [cx13,cot3]=max(x3);
    
    % Ve do thi
    for i=1:m+1
        phi11(i,1) = phi1(i,1)*180/pi;
        phi22(i,1) = phi2(i,1)*180/pi;
        phi33(i,1) = phi3(i,1)*180/pi;
    end
    
    
    figure('Name',strcat('Do thi 3D p voi cout = ',strindex),'NumberTitle','off')
    [X,Y]=meshgrid(1:1:91,phi11(:,1));
    Z=p1;
    mesh(X,Y,Z)
    hold on
    [X2,Y2]=meshgrid(1:1:91,phi22(:,1));
    Z2=p2;
    mesh(X2,Y2,Z2)
    mesh(X,Y,Z)
    [X3,Y3]=meshgrid(1:1:91,phi33(:,1));
    Z3=p3;
    mesh(X3,Y3,Z3)
    colorbar
    hold off
    
    
    figure('Name',strcat('Do thi 2D phi voi count = ',strindex),'NumberTitle','off')
    plot(transpose(phi11),max(transpose(p1)))
    hold on
    plot(transpose(phi22),max(transpose(p2)))
    plot(transpose(phi33),max(transpose(p3)))
    xlabel('circumferential coordinate');
    ylabel('dimensionless pressure')
    grid on
    hold off
    
    figure('Name',strcat('Do thi 2D h voi count = ',strindex),'NumberTitle','off')
    plot(h1,'DisplayName','h1');
    hold on
    plot(h2,'DisplayName','h2');
    plot(h3,'DisplayName','h3');
    hold off
    
    
    figure('Name',strcat('Do thi 3D h voi count = ',strindex),'NumberTitle','off')
    [X,Y]=meshgrid(1:1:91,phi11(:,1));
    ZH=h1;
    mesh(X,Y,ZH)
    hold on
    [X2,Y2]=meshgrid(1:1:91,phi22(:,1));
    ZH2=h2;
    mesh(X2,Y2,ZH2)
    %mesh(X,Y,ZH2)
    [X3,Y3]=meshgrid(1:1:91,phi33(:,1));
    ZH3=h3;
    mesh(X3,Y3,ZH3)
    colorbar
    %end
    hold off
    
    % Luu lai toan bo cac gia tri
    phi_all_1(:,1,count)=phi11;
    phi_all_2(:,1,count)=phi22;
    phi_all_3(:,1,count)=phi33;
    p_all_1(:,:,count) = p1;
    p_all_2(:,:,count) = p2;
    p_all_3(:,:,count) = p3;
    h_all_1(:,:,count) = h1;
    h_all_2(:,:,count) = h2;
    h_all_3(:,:,count) = h3;
    
    clearvars -except exilon phi_all_1 phi_all_2 phi_all_3 p_all_1 p_all_2 p_all_3 h_all_1 h_all_2 h_all_3 count;
end


figure('Name','Do thi 3D p TONG','NumberTitle','off')
hold on
for count=1:size(exilon,2)
    phi1_all = phi_all_1(:,1,count);
    phi2_all = phi_all_2(:,1,count);
    phi3_all = phi_all_3(:,1,count);
    [X,Y]=meshgrid(1:1:91,phi1_all(:,1));
    Z=p_all_1(:,:,count);
    mesh(X,Y,Z)
    
    [X2,Y2]=meshgrid(1:1:91,phi2_all(:,1));
    Z2=p_all_2(:,:,count);
    mesh(X2,Y2,Z2)
    mesh(X,Y,Z)
    [X3,Y3]=meshgrid(1:1:91,phi3_all(:,1));
    Z3=p_all_3(:,:,count);
    mesh(X3,Y3,Z3)
    colorbar
    
end
hold off




figure('Name','Do thi 2D phi TONG','NumberTitle','off')
hold on
grid on
for count=1:size(exilon,2)
    phi1_all = phi_all_1(:,1,count);
    phi2_all = phi_all_2(:,1,count);
    phi3_all = phi_all_3(:,1,count);
    plot(transpose(phi1_all),max(transpose(p_all_1(:,:,count))))
    plot(transpose(phi2_all),max(transpose(p_all_2(:,:,count))))
    plot(transpose(phi3_all),max(transpose(p_all_3(:,:,count))))
    xlabel('circumferential coordinate');
    ylabel('dimensionless pressure')
    
end
hold off




figure('Name','Do thi 2D h TONG','NumberTitle','off')
hold on
for count=1:size(exilon,2)
    plot(h_all_1(:,:,count),'g','DisplayName','h1');
    plot(h_all_2(:,:,count),'r','DisplayName','h2');
    plot(h_all_3(:,:,count),'c','DisplayName','h3');
end
hold off


figure('Name','Do thi 3D h TONG','NumberTitle','off')
hold on
for count=1:size(exilon,2)
    phi1_all = phi_all_1(:,1,count);
    phi2_all = phi_all_2(:,1,count);
    phi3_all = phi_all_3(:,1,count);
    [X,Y]=meshgrid(1:1:91,phi1_all(:,1));
    ZH=h_all_1(:,:,count);
    mesh(X,Y,ZH)
    [X2,Y2]=meshgrid(1:1:91,phi2_all(:,1));
    ZH2=h_all_2(:,:,count);
    mesh(X2,Y2,ZH2)
    %mesh(X,Y,ZH2)
    [X3,Y3]=meshgrid(1:1:91,phi3_all(:,1));
    ZH3=h_all_3(:,:,count);
    mesh(X3,Y3,ZH3)
    colorbar
end
hold off
