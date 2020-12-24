JSFP_zero = zeros(10,1);
for i=1:10
    JSFP_zero(i) = sum(JSFP_max(i,:)==0)/100;
end

figure(1)
bar(JSFP_zero)
hold on
plot(JSFP_zero,'LineWidth',2)
xlabel('b_{max}')
ylabel('JSFP收敛到不可行解比例')