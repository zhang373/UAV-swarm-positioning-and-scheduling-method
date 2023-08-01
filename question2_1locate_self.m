%根据已知三点确定自己的位置
clear;
close all;
PT = zeros([15,2]);
count = 0;
for i = 1:5
    for j = 1:i
        count = count+1;
        PT(count,:)=[-25*(i-1)+50*(j-1),25*sqrt(3)*(i-1)];
    end
end

stable = [1,11,15];
FY =[...
         0         0;
   -6.9447   21.8518;
   23.7085   28.7114;
  -57.8545   96.3084;
  -19.8235   73.3296;
   31.8933   68.3068;
  -85.3896  152.7241;
  -26.7339  124.1420;
    5.7539  158.0393;
   69.4473  150.5731;
 -100.0000  173.2051;
  -43.0805  165.8017;
   22.6309  190.2962;
   47.8973  192.0437;
  100.0000  173.2051];

PP = zeros([15,2]);
for i = 1:15
    if ~ismember(i,stable)
        alpha1 = get_angle(FY(1,1),FY(1,2),FY(i,1),FY(i,2),FY(15,1),FY(15,2));
        alpha2 = get_angle(FY(1,1),FY(1,2),FY(i,1),FY(i,2),FY(11,1),FY(11,2));
        alpha3 = get_angle(FY(11,1),FY(11,2),FY(i,1),FY(i,2),FY(15,1),FY(15,2));
        PP(i,:) = get_position2(200,pi/3,200,2*pi/3,alpha1,alpha2,alpha3);
    else
        PP(i,:)=PT(i,:);
    end
end

disp([FY,PT,PP]);

function pos_list = get_position2(rho1,theta1,rho2,theta2,alpha1,alpha2,alpha3)
x1 = rho1 * cos(theta1);
y1 = rho1 * sin(theta1);
x2 = rho2 * cos(theta2);
y2 = rho2 * sin(theta2);
[rho_rs, theta_rs] = find_intersection2(rho1, theta1, alpha1, rho2, theta2, alpha2);
[~, sz] = size(rho_rs);
in_list = [];
for j = 1:sz
    x = rho_rs(j) * cos(theta_rs(j));
    y = rho_rs(j) * sin(theta_rs(j));
    alpha3_test = acos(((x-x1)*(x-x2)+(y-y1)*(y-y2))/sqrt(((x-x1)^2+(y-y1)^2)*((x-x2)^2+(y-y2)^2)));
    in_list = [in_list;[abs(alpha3_test-alpha3),x,y]];
end
in_list = sortrows(in_list);
pos_list = in_list(1,2:3);
end

function [rhox,thetax] = find_intersection2(rho1,theta1,alpha1,rho2,theta2,alpha2)
rhox = [];
thetax = [];
xs = linspace(0,pi,5400);
ys = rho1*sin(alpha1+inferior_angle(xs, theta1))*sin(alpha2)-rho2*sin(alpha2+inferior_angle(xs, theta2))*sin(alpha1);
for i = 2:5400
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

function y = inferior_angle(angle1, angle2)
y = pi-abs(pi-abs(angle1-angle2));
end

function a = get_angle(x1,y1,x2,y2,x3,y3)
    dx1 = x1-x2;
    dy1 = y1-y2;
    dx2 = x3-x2;
    dy2 = y3-y2;
    a = acos((dx1*dx2+dy1*dy2)/sqrt((dx1^2+dy1^2)*(dx2^2+dy2^2)));
end


