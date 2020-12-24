clear
clc

%% 导弹与目标位置初始化
Nm = 20;
Nt = 10;
adiN = Nm-Nt;
model = DynamicMissileAndTarget(Nm,Nt);
model = RandomSetSituation(model);
bmax = 9;
model.target_require_num_list = ones(1,model.num_targets);
    assign_flag = zeros(model.num_targets,1);
    leader = randi(model.num_targets);
    model.target_require_num_list(leader) = model.target_require_num_list(leader) + bmax;
    assign_flag(leader) = 1;
    residual = adiN - bmax;
    for i = 1:model.num_targets
        if assign_flag(i) == 0 && residual>0
            adiMissiles = randi(residual);
            while adiMissiles>bmax
                adiMissiles = randi(residual);
            end
            model.target_require_num_list(i) = model.target_require_num_list(i) + adiMissiles;
            residual = residual - adiMissiles;
        end
    end
    % 虚拟目标增广序列
    k=1;
    model.order_targets = zeros(model.num_missiles,1);
    for i=1:model.num_targets
        for j=1:model.target_require_num_list(i)
            model.order_targets(k) = i;
            k = k+1;
        end
    end
    
    % 目标价值
    model.Targets.value = 1000*ones(model.num_targets,1);
    model.Targets.value = model.target_require_num_list'.*model.Targets.value;
    model.extend_value = model.Targets.value(model.order_targets);

%% 分配
Epis = 500;

Times = 20;
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
    
    % HCGSAP算法
    e=3;
    [HCGSAP_assign, HCGSAP_global_save, HCGSAP_time] = HCGFunctionSAP(e,model,Epis);
    HCGSAP_max(r) = max(HCGSAP_global_save);
    
    % HCG算法
    [HCGDMEA_assign, HCGDMEA_global_save, HCGDMEA_time] = HCGFunction(model,Epis);
    HCGDMEA_max(r) = max(HCGDMEA_global_save);
    
    
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
    
%     [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
%     % IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);
%     IGA_max(r) = max(IGA_global_save);
    

    
    
end
    % 归一化
    maxU = max([JSFP_max,GRM_max,SAP_max,HCGSAP_max,HCGDMEA_max,PSO_max,ACO_max,SAGA_max,GA_max]);
    JSFP_max = JSFP_max/maxU;
    GRM_max = GRM_max/maxU;
    SAP_max = SAP_max/maxU;
    PSO_max = PSO_max/maxU;
    ACO_max = ACO_max/maxU;
    SAGA_max = SAGA_max/maxU;
    HCGSAP_max = HCGSAP_max/maxU;
    HCGDMEA_max = HCGDMEA_max/maxU;
    GA_max = GA_max/maxU;
%     IGA_max = IGA_max/maxU;

algo_cell = {'JSFP';'GRM';'SAP';'HCGSAP';'HCGDMEA';'ACO';'PSO';'GA';'SAGA'};
[l,~]=size(algo_cell);
filename = [num2str(model.num_missiles),'vs',num2str(model.num_targets),'/MonteCarlo/'];
% for i=1:l
%     save([filename,algo_cell{i},'_max'],[algo_cell{i},'_max']);
%     %     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end