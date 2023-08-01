%FY02已经到达标准位置
clear;
close all;
R = 100;
target = zeros([9,2]);

for i = 1:9
    target(i,:) = [R,2*(i-1)*pi/9];
end

position = [
    100,0;
    100,360/9;
    112,80.21;
    105,119.75;
    98,159.86;
    112,199.96;
    105,240.07;
    98,280.17;
    112,320.28];

position(:,2) = position(:,2)/180*pi;
target_pol = target;
[target(:,1),target(:,2)] = pol2cart(target(:,2),target(:,1));
[position(:,1),position(:,2)] = pol2cart(position(:,2),position(:,1));

for step = 1:20
    selected = [1,2];
    for i = 3:9
        b_warning = false;
        p2 = position(i,:);
        p3 = position(selected(1),:);
        p4 = position(selected(2),:);
        a1 = get_angle(0,0,...
            p2(1),p2(2),p3(1),p3(2));
        a2 = get_angle(0,0,...
            p2(1),p2(2),p4(1),p4(2));
        a3 = get_angle(p4(1),p4(2),...
            p2(1),p2(2),p3(1),p3(2));
        r1 = target_pol(selected(1),:);
        r2 = target_pol(selected(2),:);
        gp = get_position(r1(1),r1(2),r2(1),r2(2),a1,a2,a3);
        if size(gp,1)~=1
            b_warning = true;
        end
        if ~b_warning
            position(i,:)=position(i,:)+(target(i,:)-gp)/6;
        end
    end
    if mod(step,5)==0 || step==1
        subplot(1,5,floor(step/5)+1);
        hold on;
        scatter(position(:,1),position(:,2));
        scatter(target(:,1),target(:,2),'x');
        hold off;
        axis([-120,120,-120,120]);
        title(['Step:',num2str(step)]);
    end
end
disp([position,target,abs(position-target)]);
function a = get_angle(x1,y1,x2,y2,x3,y3)
    dx1 = x1-x2;
    dy1 = y1-y2;
    dx2 = x3-x2;
    dy2 = y3-y2;
    a = acos((dx1*dx2+dy1*dy2)/sqrt((dx1^2+dy1^2)*(dx2^2+dy2^2)));
end

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