clear
clc
load('assign_save_two');
load('missiles_save_two');
load('targets_save_two');

[~,Nm,~] = size(missiles_save.p);
[~,Nt,~] = size(targets_save.p);

t = 7232;
C =linspecer(7,'qualitative');
axes('NextPlot','replacechildren','ColorOrder',C); 
%% 初始时刻
figure(1)
x = missiles_save.p(t,1,1);
y = missiles_save.p(t,1,2);
scatter(x,y,'b'); hold on;
quiver(x,y, cos(missiles_save.angle(t,1))',sin(missiles_save.angle(t,1))',5,'color','b');
text(x+1,y+1,num2str(1));

for i=1:Nm
    x = missiles_save.p(t,i,1);
    y = missiles_save.p(t,i,2);
    l1=scatter(x,y,'b');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(missiles_save.angle(t,i))',sin(missiles_save.angle(t,i))',5,'color','b');
end

for i=1:Nt
    x = targets_save.p(t,i,1);
    y = targets_save.p(t,i,2);
    l2=scatter(x,y,'r');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(targets_save.angle(t,i))',sin(targets_save.angle(t,i))',5,'color','r');
end

for i=1:Nm
    x1 = missiles_save.p(t,:,1);
    y1 = missiles_save.p(t,:,2);
    x2 = targets_save.p(t,:,1);
    y2 = targets_save.p(t,:,2);
    target = assign_save(t,i);
    plot([x1(i),x2(target)],[y1(i),y2(target)],'color',C(3,:));
end

for i=1:Nm
    plot(missiles_save.p(1:t,i,1),missiles_save.p(1:t,i,2),'b');
end
for i=1:Nt
    plot(targets_save.p(1:t,i,1),targets_save.p(1:t,i,2),'r');
end
legend([l1,l2],'导弹','目标')

figure(2);


% 绘制导弹位置
x = missiles_save.p(1,1,1);
y = missiles_save.p(1,1,2);
scatter(x,y,'b'); hold on;
quiver(x,y, cos(missiles_save.angle(1,1))',sin(missiles_save.angle(1,1))',5,'color','b');
text(x+1,y+1,num2str(1));
% plot(x+r1*cos(theta),y+r1*sin(theta),'g--');
% plot(x+r2*cos(theta),y+r2*sin(theta),'b--');

for i=1:Nm
    x = missiles_save.p(1,i,1);
    y = missiles_save.p(1,i,2);
    scatter(x,y,'b');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(missiles_save.angle(1,i))',sin(missiles_save.angle(1,i))',5,'color','b');
%     plot(x+r1*cos(theta),y+r1*sin(theta),'g--');
%     plot(x+r2*cos(theta),y+r2*sin(theta),'b--');
end

% 绘制目标位置
for i=1:Nt
    x = targets_save.p(1,i,1);
    y = targets_save.p(1,i,2);
    scatter(x,y,'r');
    text(x+1,y+1,num2str(i));
    quiver(x,y, cos(targets_save.angle(1,i))',sin(targets_save.angle(1,i))',5,'color','r');
end

for i=1:Nm
    l1=plot(missiles_save.p(:,i,1),missiles_save.p(:,i,2),'b');
end
for i=1:Nt
    l2=plot(targets_save.p(1:t,i,1),targets_save.p(:,i,2),'r');
end
legend([l1,l2],'导弹','目标')

