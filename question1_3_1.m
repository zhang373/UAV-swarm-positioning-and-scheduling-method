%针对FY02的调整方案,发射定位信号的无人机根据初始两个无人机的信号调整位置
%C通过夹角ACB调整自己的位置
clear;
close all;
FY =[...
   -6.9447   21.8518;
  -85.3896  152.7241;
  -26.7339  124.1420;
    5.7539  158.0393;
   69.4473  150.5731;
  -43.0805  165.8017;
   22.6309  190.2962;
   47.8973  192.0437];
res = zeros([8,3]);
for i = 1:8
    A = [-100,0];
    B = [100,0];
    C = FY(i,:);
    
    min_theta = 0;
    max_theta = 0;
    
    d = 2;
    
    alpha0 = get_angle(A(1),A(2),C(1),C(2),B(1),B(2));
    y_old = get_angle(A(1),A(2),C(1)+d,C(2),B(1),B(2));
    alpha_ab = zeros([2,1]);
    ii=1;
    as = zeros([360,1]);
    for theta = linspace(0,2*pi,360)
        alpha = get_angle(A(1),A(2),C(1)+d*cos(theta),C(2)+d*sin(theta),B(1),B(2));
        y=alpha-alpha0;
        as(ii)=y;
        ii=ii+1;
    end
    theta = linspace(0,2*pi,360);
    asd = diff(sign(as));
    t = theta(asd~=0);
    theta = (t(1)+t(2))/2;
    
    [real_theta,~] = cart2pol(C(1),C(2));
    theta = mod(theta/pi*180,180);
    real_theta = mod(real_theta/pi*180,180);
    delta_theta = abs(theta-real_theta);
    res(i,:)=[theta,real_theta,delta_theta];
end
disp(res);
function a = get_angle(x1,y1,x2,y2,x3,y3)
    dx1 = x1-x2;
    dy1 = y1-y2;
    dx2 = x3-x2;
    dy2 = y3-y2;
    a = acos((dx1*dx2+dy1*dy2)/sqrt((dx1^2+dy1^2)*(dx2^2+dy2^2)));
end
