clear;
close all;
theta = linspace(0,2*pi,360);
rho = linspace(80,120,64);
R = 100;

Eri = @(alpha_0,alpha_i)abs(sin(alpha_i).*(1./sin(alpha_i)-1./sin(alpha_0)));

[U,V] = meshgrid(rho,theta);

[X,Y] = pol2cart(V,U);
cos_a2 = @(x1,y1,x,y)(-(x1-x).*x-(y1-y).*y).^2./(((x1-x).^2+(y1-y).^2).*(x.^2+y.^2));
Z = real(cos_a2(R,0,X,Y));

hold on;
%mesh(X,Y,Z);
xlabel('x');
ylabel('y');
zlabel('$\cos^2\alpha$','Interpreter','latex');
theta0 = linspace(2*pi/9,2*pi-2*pi/9,8);
x0 = R * cos(theta0);
y0 = R * sin(theta0);
z0 = cos_a2(R,0,x0,y0);
zsplit = [0,(z0(1:3)+z0(2:4))/2,1];
for i = 1:4
    n = find(zsplit(i)>Z|Z>zsplit(i+1));
    x = X;
    x(n) = inf;
    y = Y;
    y(n) = inf;
    z = Z;
    z(n) = inf;
    mesh(x,y,z);
end
for i = 1:8
    n = find((X-x0(i)).^2+(Y-y0(i)).^2<100);
    x = X(n);
    y = Y(n);
    z = Z(n);
    plot3(x,y,z,'*');
end
hold off;
tests =... 
[27.2476  107.7756;
-45.9049   86.6060;
-88.5605   36.2793;
-91.7764  -30.9440;
-49.0405  -79.1272;
24.8499  -93.0478;
79.9858  -55.9554];

rho1 = 200;
theta1 = 0;
x1 = rho1 * cos(theta1);
y1 = rho1 * sin(theta1);
r2 = 200;
theta2 = 2*pi/9;
x2 = rho1 * cos(theta2);
y2 = r2 * sin(theta2);

xs = zeros([7,1]);
ys = zeros([7,1]);
 
txs = zeros([7,1]);
tys = zeros([7,1]);

for i = 1:7
    x = tests(i,1);
    y = tests(i,2);
    xs(i) = x;
    ys(i) = y;
    %模拟接受的信号
    alpha1 = acos((x*(x-x1)+y*(y-y1))/sqrt((x^2+y^2)*((x-x1)^2+(y-y1)^2)));
    alpha2 = acos((x*(x-x2)+y*(y-y2))/sqrt((x^2+y^2)*((x-x2)^2+(y-y2)^2)));
    alpha3 = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)^2+(y-y1)^2)*((x-x2)^2+(y-y2)^2)));
    pos = get_position(rho1,theta1,r2,theta2,alpha1,alpha2,alpha3);
    txs(i) = pos(1,1);
    tys(i) = pos(1,2);
end
kk=[xs,ys,txs,tys,abs((xs-txs)),abs((ys-tys))];
fprintf('第一种方法\n');
disp(kk);

res = zeros([7,6]);
R = 100;
x1 = R;
y1 = 0;
x2 = R * cos(4*2*pi/9);
y2 = R * sin(4*2*pi/9);
tests =... 
[27.2476  107.7756;
-45.9049   86.6060;
-88.5605   36.2793;
-91.7764  -30.9440;
-49.0405  -79.1272;
24.8499  -93.0478;
79.9858  -55.9554];
for i = 1:7
    x = tests(i,1);
    y = tests(i,2);
    cos_a21 = ((x1-x)*x+(y1-y)*y)^2/(((x1-x)^2+(y1-y)^2)*(x^2+y^2));
    cos_a22 = ((x2-x)*x+(y2-y)*y)^2/(((x2-x)^2+(y2-y)^2)*(x^2+y^2));
    cos_a23 = ((x1-x)*(x2-x)+(y1-y)*(y2-y))^2/...
        (((x1-x)^2+(y1-y)^2)*((x2-x)^2+(y2-y)^2));
    rs = search_position(4*2*pi/9,acos(cos_a21),acos(cos_a22),acos(cos_a23),false);
    res(i,:)=[x,y,rs(1),rs(2),abs((rs(1)-x)),abs((rs(2)-y))];
end

fprintf('第二种方法\n');
disp(res);

function pos_list = get_position(rho1,theta1,rho2,theta2,alpha1,alpha2,alpha3)
%获取飞行器的位置
eps = 0.01;
%转化为直角坐标
x1 = rho1 * cos(theta1);
y1 = rho1 * sin(theta1);
x2 = rho2 * cos(theta2);
y2 = rho2 * sin(theta2);
[rho_rs, theta_rs] = find_intersection(rho1, theta1, alpha1, rho2, theta2, alpha2);
[~, sz] = size(rho_rs);
in_list = [];
all_list = [];
%筛选解
for j = 1:sz
    x = rho_rs(j) * cos(theta_rs(j));
    y = rho_rs(j) * sin(theta_rs(j));
    alpha3_test = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)^2+(y-y1)^2)*((x-x2)^2+(y-y2)^2)));
    if in_circle([0,0],[x1,y1],[x2,y2],[x,y]) && abs(alpha3_test-alpha3) < eps
        %较优的解优先输出
        in_list = [in_list;[abs(alpha3_test-alpha3),x,y]];
    else
        all_list = [all_list;[abs(alpha3_test-alpha3),x,y]];
    end
end
if size(in_list,1)~=0
    in_list = sortrows(in_list);
    pos_list = in_list(1,2:3);
else
    all_list = sortrows(all_list);
    pos_list = all_list(:,2:3);
end
end

function [rhox,thetax] = find_intersection(rho1,theta1,alpha1,rho2,theta2,alpha2)
%求交点函数，根据两个在圆上的飞行器的极坐标和两个飞行器与自身与原点之间的夹角确定交点，可能有多个结果，返回交点极坐标
rhox = [];
thetax = [];
xs = linspace(0,2*pi,10800);
ys = rho1*sin(alpha1+inferior_angle(xs, theta1))*sin(alpha2)-rho2*sin(alpha2+inferior_angle(xs, theta2))*sin(alpha1);
%搜索解,theta范围:[0,2*pi]
for i = 2:10800
    if sign(ys(i))~=sign(ys(i-1))
        x = xs(i);
        if alpha1~=0
            new_rx = rho1*sin(alpha1+inferior_angle(x, theta1))/sin(alpha1);
        else
            new_rx = rho2*sin(alpha2+inferior_angle(x, theta2))/sin(alpha2);
        end
        rhox = [rhox, new_rx];
        thetax = [thetax, x];

    end
end
end

function b = in_circle(r,r1,r2,r3)
%判断是否在外接圆内
x = [r1(1);r2(1);r3(1)];
y = [r1(2);r2(2);r3(2)];
s = [sum(r1.*r1);sum(r2.*r2);sum(r3.*r3)];
i = [1;1;1];

det_x = det([i,s,y]);
det_y = det([i,x,s]);
det_2d = 2*det([i,x,y]);

rc = [det_x,det_y]/det_2d;

b = norm(r-rc) < norm(r1-rc);
end

function y = inferior_angle(angle1, angle2)
%求两个方位角之间的小于180度的角
y = pi-abs(pi-abs(angle1-angle2));
end

function result = search_position(thetap,alpha0,alphap,alpha3,back)
eps = 0.3;

R = 100;
n = 8;

Eri = @(alpha_0,alpha_i)abs(cos(alpha_0).^2-cos(alpha_i).^2);

theta = linspace(2*pi/(n+1),2*pi-2*pi/(n+1),n);

alpha = pi/2-inferior_angle(theta,0)/2;
eri = [Eri(alpha0,alpha(1:floor((n+1)/2)));2:floor((n+3)/2)];
eri = sortrows(eri',1);
fy0k = eri(1:3,2);
if back
    fy0k = n+3-fy0k;
end
if (alpha0 > alphap && alpha0 > alpha3) || (alpha3<pi && alpha3>alpha0 && alpha3>alphap)
   fy0k = n+3-fy0k;
end
b = false;
result = [nan,nan,nan,nan];
for i = 1:3
    if b
        break;
    end
    gf0 = grid_find(R,0,R,thetap,R,theta(fy0k(i)-1),0,10,0,2*pi,alpha0,alphap);
    gf0 = gf0(1:3,:);
    for j = 1:3
        if b
            break;
        end
        r = gf0(j,2);
        t = gf0(j,3);
        gf1 = grid_find(R,0,R,thetap,R,theta(fy0k(i)-1), ...
            r-1,r+1,t-pi/8,t+pi/8,alpha0,alphap);
        gf1 = gf1(1:3,:);
        for k = 1:3
            if gf1(k,1) < eps
                b = true;
                rr = gf1(k,2:3);
                result = [R*cos(theta(fy0k(i)-1))+rr(1)*cos(rr(2)),...
                    R*sin(theta(fy0k(i)-1))+rr(1)*sin(rr(2))];
                break
            end
        end
    end
end
end

function M = grid_find(rho1,theta1,rhop,thetap,testrho,testtheta,rho_left,rho_right,theta_left,theta_right,alpha0,alphap)
%在圆形或扇形网格中搜索
M = zeros([25,3]);
i = 1;
[x1,y1] = pol2cart(theta1,rho1);
[xp,yp] = pol2cart(thetap,rhop);
[testx,testy] = pol2cart(testtheta,testrho);
for rho = linspace(rho_left,rho_right,9)
    for theta = linspace(theta_left,theta_right,16)
        [x,y] = pol2cart(theta,rho);
        x = x + testx;
        y = y + testy;
        cos_a2 = ((x1-x)*x+(y1-y)*y)^2/(((x1-x)^2+(y1-y)^2)*(x^2+y^2));
        cos_a2p = ((xp-x)*x+(yp-y)*y)^2/(((xp-x)^2+(yp-y)^2)*(x^2+y^2));
        M(i,1) = abs(cos(alpha0)^2-cos_a2)+abs(cos(alphap)^2-cos_a2p);
        M(i,2) = rho;
        M(i,3) = theta;
        i = i + 1;
    end
end
M = sortrows(M);
end