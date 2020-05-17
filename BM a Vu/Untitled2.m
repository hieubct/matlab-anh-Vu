clear
clc
v=[4000 3000 2000 1000];
v0=[v 830];
v2=[v 540];
v4=[v 420];
v6=[v 340];
v8=[v 230];
n0=[0.2 0.25 0.357 0.614 0.7];
n2=[0.2 0.218 0.275 0.44 0.7];
n4=[0.2 0.21 0.24 0.36 0.7];
n6=[0.2 0.208 0.23 0.3 0.7];
n8=[0.2 0.206 0.215 0.24 0.7];

Y=[v0;v2;v4;v6;v8]
X=[n0;n2;n4;n6;n8]
line=['r','b','g','k','m'];
figure(1)
grid on
hold on
for i=1:4
    xi = linspace(min(X(i,:)), max(X(i,:)), 150);                     % Evenly-Spaced Interpolation Vector
    yi = interp1(X(i,:), Y(i,:), xi, 'spline', 'extrap')

    plot(X(i,:), Y(i,:), 'bp')

    plot(xi, yi,line(i))
end

hold off

% xlabel('X')
% ylabel('Y')
% legend('Original Data', 'Interpolation', 'Location', 'NW')
