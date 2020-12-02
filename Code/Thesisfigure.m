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
xlabel('���˳�Ա��')
ylabel('����ر�')
grid on
figure(2)
plot(x,u,'b.');
xlabel('���˳�Ա��')
ylabel('���������ûر�')
grid on