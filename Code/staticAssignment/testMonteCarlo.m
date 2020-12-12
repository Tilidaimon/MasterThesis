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
for r = 1:Times
    
    % JSFP算法
    [JSFP_assign, JSFP_global_save, JSFP_time] = JSFPfunction(model,Epis);
    JSFP_max(r) = max(JSFP_global_save);
    
    % GRMFMI算法
    [GRM_assign, GRM_global_save] = GRMfunction(model,Epis);
    GRM_max(r) = max(GRM_global_save);
    
    % SAP算法
    [SAP_lagrange2_assign, SAP_global_save, SAP_time] = SAP_lagrange_function2(model,Epis);
    SAP_max(r) = max(SAP_global_save);
    
    % PSO算法
    [PSO_extend_assign, PSO_global_save, PSO_time] = PSOfunction(model);
    %PSO_assign = DecodeFromExtendPlan(model,PSO_extend_assign);
    PSO_max(r) = max(PSO_global_save);
    
    
    numLimit = max(model.num_missiles, sum(model.target_require_num_list));
    s = zeros(20,numLimit); % 预分配内存
    for i=1:20
        s(i,:) = randperm(numLimit);
    end
    
    % SAGA算法
    [SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time] = SAGAfunction(s,model,20,200,0.4,0.4, "PMX", "EM",100,0.95);
    % SAGA_assign = DecodeFromExtendPlan(model,SAGA_extend_plan);
    SAGA_max(r) = max(SAGA_global_save);
    
    % ACOVFSAGA
    popsize=20;steps=200;Pcross=0.4;Pmutate=0.4;
    Top = 50; Dop = 0.95;
    [ACOVFSAGA_plan, ACO_global_save,ACO_mean_save, ACO_time] = ACOVFSAGA(model,2*popsize,steps,Pcross,Pmutate,Top,Dop);
    % ACO_plan = DecodeFromExtendPlan(model,ACOVFSAGA_plan);
    ACO_max(r) = max(ACO_global_save);
    
    % GA
    [GA_extend_plan,GA_global_save,GA_mean_save, GA_time] = GA1(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
    % GA_assign = DecodeFromExtendPlan(model,GA_extend_plan);
    GA_max(r) = max(GA_global_save);
    
    [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
    % IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);
    IGA_max(r) = max(IGA_global_save);
    
    % 归一化
    maxU = max([JSFP_max(r),GRM_max(r),SAP_max(r),PSO_max(r),ACO_max(r),SAGA_max(r),GA_max(r),IGA_max(r)]);
    JSFP_max(r) = JSFP_max(r)/maxU;
    GRM_max(r) = GRM_max(r)/maxU;
    SAP_max(r) = SAP_max(r)/maxU;
    PSO_max(r) = PSO_max(r)/maxU;
    ACO_max(r) = ACO_max(r)/maxU;
    SAGA_max(r) = SAGA_max(r)/maxU;
    GA_max(r) = GA_max(r)/maxU;
    IGA_max(r) = IGA_max(r)/maxU;
    
end

algo_cell = {'JSFP';'GRM';'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
[l,~]=size(algo_cell);
filename = [num2str(model.num_missiles),'vs',num2str(model.num_targets),'/MonteCarlo/'];
for i=1:l
    save([filename,algo_cell{i},'_max'],[algo_cell{i},'_max']);
    %     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
end