clear;
close all;
target_i = 3;
target_cos_alpha1 = 0.8;
target_cos_alpha2 = 0.2;

R = 100;
theta = linspace(0,2*pi,10);
theta = theta(2:9);

r1 = [R;0];
r = [R*cos(theta);R*sin(theta)];

M = zeros(2,8,8);
%生成张量M
for i = 1:8
    rr = -r(:,i);
    dr1 = r1 - r(:,i);
    for j = 1:8
        if i~=j
            drx = r(:,j) - r(:,i);
            M(1,i,j) = drx'*rr/(norm(drx)*norm(rr));
            M(2,i,j) = drx'*dr1/(norm(drx)*norm(dr1));
        else
            M(1,i,j) = inf;
            M(2,i,j) = inf;
        end
    end
end

hold on;
M=permute(M,[3,1,2]);
for i = 1:8
    subplot(2,4,i);
    plot(M(:,1,i),M(:,2,i),'*');
    for k = 1:8
        if i~=k
            text(M(k,1,i),M(k,2,i),num2str(k));
        end
    end
    xlabel('$\cos^2{\alpha_2}$','Interpreter','latex');
    ylabel('$\cos^2{\alpha_3}$','Interpreter','latex');
    title(['FY0',num2str(i+1)]);
end
hold off;

d = inf;
target_j = 0;
for j = 1:8
    if j == target_i
        continue
    end
    t = (target_cos_alpha1-M(j,1,target_i))^2+...
        (target_cos_alpha2-M(j,2,target_i))^2;
    if t<d
        d = t;
        target_j = j;
    end
end

target_j