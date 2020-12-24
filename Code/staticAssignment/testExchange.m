clear
clc

%% 导弹与目标位置初始化
Nm = 10;
Nt = 10;

model = DynamicMissileAndTarget(Nm,Nt);
model = RandomSetSituation(model);

%% 分配
Epis = 500;

Times = 100;
for r=1:Times
    [SAPEx_assign, SAPEx_global_save, SAPEx_time(r)] = SAP_lagrange_function2(model,Epis);
    SAPEx_max(r) = SAPEx_global_save(end);
    
    [SAP_assign, SAP_global_save, SAP_time(r)] = SAP_lagrange_function_noEx(model,Epis);
    SAP_max(r) = SAP_global_save(end);
    
    numLimit = max(model.num_missiles, sum(model.target_require_num_list));
    s = zeros(20,numLimit); % 预分配内存
    for i=1:20
        s(i,:) = randperm(numLimit);
    end

% SAGA算法
    [SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time(r)] = SAGAfunction(s,model,20,500,0.4,0.4, "PMX", "EM",100,0.95);
    SAGA_max(r) = SAGA_global_save(end);
%     [JSFP_assign, JSFP_global_save, JSFP_time(r)] = JSFPfunction(model,Epis);
%     JSFP_max(r) = JSFP_global_save(end);
%     
%     [JSFPEx_assign, JSFPEx_global_save, JSFPEx_time(r)] = JSFPfunction_Ex(model,Epis);
%     JSFPEx_max(r) = JSFPEx_global_save(end);
% 
%     [GRM_assign, GRM_global_save, GRM_time(r)] = GRMfunction(model,Epis);
%     GRM_max(r) = GRM_global_save(end);
%     
%     [GRMEx_assign, GRMEx_global_save, GRMEx_time(r)] = GRMfunction(model,Epis);
%     GRMEx_max(r) = GRMEx_global_save(end);
end

maxU = max([SAPEx_max,SAP_max,SAGA_max]);

SAP_max = SAP_max/maxU;
SAPEx_max = SAPEx_max/maxU;
% JSFP_max = JSFP_max/maxU;
% JSFPEx_max = JSFPEx_max/maxU;
% GRM_max = GRM_max/maxU;
% GRMEx_max = GRMEx_max/maxU;
figure(1)
boxplot([SAP_max;SAPEx_max]')
figure(2)
boxplot([SAP_time;SAPEx_time]')