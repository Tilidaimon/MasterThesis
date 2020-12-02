clear
clc

rmin = 10;
e = 3;

for x=1:100
    y(x) = (rmin*log(x+e-1)/log(e));
    u(x) = y(x)/x;
end

x=1:100;
figure(1)
plot(x,y,'r.');
xlabel('联盟成员数')
ylabel('任务回报')
grid on
figure(2)
plot(x,u,'b.');
xlabel('联盟成员数')
ylabel('智能体所得回报')
grid on