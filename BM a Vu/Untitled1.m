clc
clear
m=360;
n=57;
e=0.5;
deltaphi = pi/360;
deltalanda=1/57;
beta = 1.2; % h? s? l?p
H = zeros(1,720);
% tinh phan bo do day mang dau
for i = 1:720
    a = (i-1)*deltaphi;
    H(i) = 1 + e*cos(a);
  
end
% tinh toan cac he so cho phuong trinh
for i=2:359
    A(1,i)=H(i+1)^3;%A
    A(2,i)=H(i-1)^3;%B
    A(3,i)=(deltaphi/deltalanda)^2*H(i)^3;%C
    A(4,i)=A(3,i);%D
    A(5,i)=A(1,i)+A(2,i)+A(3,i)+A(4,i);%E
    A(6,i)=3*deltaphi*(H(i+1)-H(i-1));%F  
end
A(1,360)=H(1)^3;
A(2,360)=H(359)^3;
A(3,360)=(deltaphi/deltalanda)^2*H(360)^3;
A(4,360)=A(3,360);
A(5,360)=A(1,360)+A(2,360)+A(3,360)+A(4,360);
A(6,360)=3*deltaphi*(H(1)-H(359));
% dinh nghia gia tri ap suat bien va ap suat ban dau
P=ones(57,359);
P1=zeros(57,1);
P2=zeros(1,361);
P=[P1 P P1];
P=[P2;P];
PN=zeros(58,361);%matran trung gian dung de tinh toan
sum=0;
sum1=0;
interactionnum=0;
delta=1;
while delta>=e^-3
    interactionnum=interactionnum+1;
    Pold=P;  
    for j=2:57
        for i=2:360
            if P(i,j)<=0
                P(i,j)=0;
                P(i,j+1)=0;
            else
                P(j,i)=(A(1,i)*P(j,i+1)+A(2,i)*PN(j,i-1)+A(3,i)*P(j+1,i)+A(4,i)*PN(j-1,i)-A(6,i))/A(5,i)-P(j,i)+P(j,i);
              
                P(j,i)=PN(j,i);
            end           
        end
    end

%tinh toan phan bo ap suat
for i=2:360
    if P(58,i)<=0
        P(58,i)=0;
        P(58,i+1)=0;
    else
        PN(58,i)=beta*((A(1,i)*P(58,i+1)+A(2,i)*PN(58,i-1)+2*A(3,i)*PN(57,i)-A(6,i))/A(5,i)-P(58,i))+P(58,i);
                P(58,i)=PN(58,i); 
    end
end
for j=2:58
    for i=2:360
        sum=sum+abs(P(j,i)-Pold(j,i));
        sum1=sum1+abs(P(j,1));
    end
end
delta=sum/sum1;
end
% tinh toan phan doi xung con lai
Pdau=P(1:57,:);
Pduoi=flipud(Pdau);
P=[Pdau;Pduoi];
% tim dieu kien bien reynol
for j=2:114
    for i=2:360
        if P(j,i)==0;
            %Rey(j-1)-i;
            break
        end
    end
end
% ve do thi
x=1:1:361;
y=1:1:114;
z=P(y,x);
figure
mesh(x,y,z);