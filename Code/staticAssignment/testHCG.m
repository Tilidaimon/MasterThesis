clear
clc

%% 导弹与目标位置初始化
Nm = 50;
Nt = 10;

model = DynamicMissileAndTarget(Nm,Nt);
model = RandomSetSituation(model);


%% 分配
Epis = 500;

Times = 10;
for e=5:20
    for r = 1:Times
        
        % HCGSAP算法
        [HCGSAP_assign, HCGSAP_global_save, HCGSAP_time] = HCGFunctionSAP(e,model,Epis);
        HCGSAP_max(e-4,r) = max(HCGSAP_global_save);
        
        [HCGSAP2_assign, HCGSAP2_global_save, HCGSAP2_time] = HCGFunction2(e,model,Epis);
        HCGSAP2_max(e-4,r) = max(HCGSAP2_global_save);
    end
end
Maxhcg = max(max(HCGSAP_max));
Maxhcg2 = max(max(HCGSAP2_max));

M = max(Maxhcg,Maxhcg2);
HCGSAP_max = HCGSAP_max/M;
HCGSAP2_max = HCGSAP2_max/M;
for i=1:16
    mean_box(i) = mean(HCGSAP_max(i,:));
    mean_box2(i) = mean(HCGSAP2_max(i,:));
end

figure(1)
boxplot(HCGSAP_max');
hold on
plot(mean_box,'.-','MarkerSize',10)

figure(2)
boxplot(HCGSAP2_max');
hold on
plot(mean_box2,'.-','MarkerSize',10)


figure(2)
e=5:20;
y = log(e+5)./log(e);
PoA = 1./min(HCGSAP_max,[],2);
    
%     % HCG算法
%     [HCGDMEA_assign, HCGDMEA_global_save, HCGDMEA_time] = HCGFunction(model,Epis);
%     HCGDMEA_max(r) = max(HCGDMEA_global_save);