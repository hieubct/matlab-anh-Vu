clear ;
m = 3;
S = 0.0109;
g = 9.81;
%t0 = 0;
%tf = 4.5;
V0 = 0;
h=0.01;
H=0;
% phuong trình vi phân dv/dt 
f = @(t, y) g -(0.5*1.225*y*y * S *(-0.0003*y + 0.4269))/m;
 
% Khoang thoi gian gian
tspan = [0:h:10];
 
% dieu kien ban dau
y0 = 0;
 
% Giai phuong trình bang ode45
[t, y] = ode45(f, tspan, y0);
 
%deltaS = y*h + 0.5*g*h*h;
deltaH = y*h + 0.5*g*h*h;
KQ= [t y deltaH];
H(1) = deltaH(1);
for i=1:length(deltaH)-1
    H = [H H(i)+deltaH(i+1)];
end
% do thi
figure(1);
plot(t, y);
xlabel('t');
ylabel('y');
title('do thi v theo t');
grid on;
figure(2);
plot(y,H);
xlabel('y');
ylabel('H');
title('do thi h theo y');
grid on;
