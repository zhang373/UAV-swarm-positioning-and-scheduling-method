%通过FY02,FY01,FY11的坐标定位FY03的坐标
%同理除FY01,FY11,FY15外所有的飞行器都可以知道全局信息
clear;
close all;

FY01 = [0,0];
FY02 = [-6.9447,21.8518];
FY03 = [23.7085,28.7114];
FY11 = [-100.0000,173.2051];

alpha1 = get_angle(FY01(1),FY01(2),FY02(1),FY02(2),FY11(1),FY11(2));
alpha2 = get_angle(FY01(1),FY01(2),FY02(1),FY02(2),FY03(1),FY03(2));
alpha3 = get_angle(FY11(1),FY11(2),FY02(1),FY02(2),FY03(1),FY03(2));

d = 10;

FY02_far = FY02+(FY02-FY01)*d/norm(FY02-FY01);
alpha2_1 = get_angle(FY01(1),FY01(2),FY02_far(1),FY02_far(2),FY03(1),FY03(2));

D = d/sin(alpha2-alpha2_1)*sin(alpha2_1);

F1 = FY01-FY02;
[theta,~] = cart2pol(F1(1),F1(2));
Ealpha_best = inf;
pos = [0,0];
for i = linspace(0,2*pi,3600)
FY03_test = FY02+[D*cos(i+theta),D*sin(i+theta)];
alpha2_test = get_angle(FY01(1),FY01(2),FY02(1),FY02(2),FY03_test(1),FY03_test(2));
alpha3_test = get_angle(FY11(1),FY11(2),FY02(1),FY02(2),FY03_test(1),FY03_test(2));
Ealpha = (alpha3_test-alpha3)^2+(alpha2_test-alpha2)^2;
if Ealpha < Ealpha_best
    Ealpha_best = Ealpha;
    pos = FY03_test;
end
end
FY03_in_FY02 = pos;
disp([FY03;FY03_in_FY02;abs(FY03-FY03_in_FY02)]);

function a = get_angle(x1,y1,x2,y2,x3,y3)
    dx1 = x1-x2;
    dy1 = y1-y2;
    dx2 = x3-x2;
    dy2 = y3-y2;
    a = acos((dx1*dx2+dy1*dy2)/sqrt((dx1^2+dy1^2)*(dx2^2+dy2^2)));
end