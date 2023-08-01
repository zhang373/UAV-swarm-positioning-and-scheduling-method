%无人机调度，防止碰撞
clear;
close all;
PT = [];
count = 0;
for i = 1:5
    for j = 1:i
        count = count+1;
        if ~ismember(count,[1,11,15])
            PT=[PT;[-25*(i-1)+50*(j-1),25*sqrt(3)*(i-1)]];
        end
    end
end

FY = [...
    -51.5217   16.7795;
   27.8559   46.1572;
  -94.8564   41.7461;
   25.6875  112.2901;
   60.1980   96.8005;
  -39.2831  165.6207;
   23.8277  178.7315;
   67.9484  172.8523;
   65.9515  120.8553;
  -99.9659  123.2392;
    4.0878  177.2929;
   20.7731  143.9781];

near = zeros([12,12]);
chozen_index = ones([12,1]);
chozen = zeros([12,1]);
for i = 1:12
    near(i,:) = find_nearest(FY(i,:),PT);
    chozen(i,1) = near(i,1);
end

subplot(2,2,1);
hold on;
scatter(FY(:,1),FY(:,2));
scatter(PT(:,1),PT(:,2),'x');
for i = 1:12
    plot([FY(i,1),PT(near(i,1),1)],[FY(i,2),PT(near(i,1),2)]);
end
hold off;
title('调整之前');
b_ok = false;
disp(chozen');
while ~b_ok
    b_ok = true;    
    for i = 1:12
        if count_times(chozen(i),chozen)>1
            b_ok = false;
            chozen_index(i)=chozen_index(i)+1;
            chozen(i)=near(i,chozen_index(i));
        end
    end
    disp(chozen');
end
subplot(2,2,2);
hold on;
scatter(FY(:,1),FY(:,2));
scatter(PT(:,1),PT(:,2),'x');
for i = 1:12
    plot([FY(i,1),PT(chozen(i),1)],[FY(i,2),PT(chozen(i),2)]);
end
hold off;
title('调整之后');
fly = [];
no = 0;
while size(fly,2)<12
    no = no+1;
    cross_count = zeros([12,1]);
    set_parent = 1:12;
    for i = 1:12
        if ismember(i,fly)
            continue
        end
        a1 = FY(i,:);
        b1 = PT(chozen(i),:);
        for j = (i+1):12
            if ismember(j,fly)
                continue
            end
            a2 = FY(j,:);
            b2 = PT(chozen(j),:);
            if sign(cross2(a2-a1,b1-a1))~=sign(cross2(b2-a1,b1-a1)) &&...
                    sign(cross2(a1-a2,b2-a2))~=sign(cross2(b1-a2,b2-a2))
                set_parent(set_parent(j))=set_parent(i);
                cross_count(i)=cross_count(i)+1;
                cross_count(j)=cross_count(j)+1;
            end
        end
    end
    max_count_id = zeros([12,1]);
    max_count = -ones([12,1]);
    for i = 1:12
        if ismember(i,fly)
            continue;
        end
        if max_count(set_parent(i))<cross_count(i)
            max_count_id(set_parent(i))=i;
            max_count(set_parent(i))=cross_count(i);
        end
    end
    add_fly = [];
    for i = 1:12
        if max_count_id(i)~=0
            add_fly = [add_fly,max_count_id(i)];
        end
    end
    fprintf('第%d批:\n',no);
    disp(add_fly);
    subplot(2,2,2+no);
    hold on;
    scatter(FY(:,1),FY(:,2));
    scatter(PT(:,1),PT(:,2),'x');
    for i = 1:size(add_fly,2)
        plot([FY(add_fly(i),1),PT(chozen(add_fly(i)),1)],[FY(add_fly(i),2),PT(chozen(add_fly(i)),2)]);
    end
    title(['第',num2str(no),'批']);
    fly = [fly,add_fly];
end
function y = find_nearest(x,A)
n = size(A,1);
y = zeros([n,2]);
for i = 1:n
    dr = norm(x-A(i,:));
    y(i,:)=[dr,i];
end
y=sortrows(y,1);
y=y(:,2);
end

function y = count_times(x,A)
n = size(A,1);
count = 0;
for i = 1:n
    if x==A(i)
        count = count + 1;
    end
end
y = count;
end

function y = cross2(p1,p2)
y  = p1(1)*p2(2)-p1(2)*p2(1);
end