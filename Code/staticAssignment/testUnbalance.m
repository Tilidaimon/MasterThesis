clear
clc

%% 导弹与目标位置初始化
Nm = 20;
Nt = 10;
adiN = Nm-Nt;

model = DynamicMissileAndTarget(Nm,Nt);
model = RandomSetSituation(model);

algo_cell = {'JSFP';'GRM';'SAP';'HCGSAP';'HCGDMEA'};
[l,~]=size(algo_cell);
filename = 'TestUnbalance/';
for bmax = 1:adiN
    
    
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
    Epis = 200;
    
    Times = 2;
    
    
    %         IGA_time = zeros(1,Times);
    for r = 1:Times
        %JSFP算法
        [JSFP_assign, JSFP_global_save, JSFP_time(r)] = JSFPfunction(model,Epis);
        JSFP_max(r) = JSFP_global_save(end);
        
        % GRMFMI算法
        [GRM_assign, GRM_global_save, GRM_time(r)] = GRMfunction(model,Epis);
        GRM_max(r) = GRM_global_save(end);
        
        % SAP算法
        tic
        [SAP_lagrange2_assign, SAP_global_save,SAP_time(r)] = SAP_lagrange_function2(model,Epis);
        SAP_max(r) = SAP_global_save(end);
        
        % HCGSAP算法
        e = 3;
        [HCGSAP_assign, HCGSAP_global_save, HCGSAP_time] = HCGFunctionSAP(e,model,Epis);
        HCGSAP_max(r) = HCGSAP_global_save(end);
        
        % HCG算法
        [HCGDMEA_assign, HCGDMEA_global_save, HCGDMEA_time] = HCGFunction(model,Epis);
        HCGDMEA_max(r) = HCGDMEA_global_save(end);
        
        % PSO算法
%         tic
%         [PSO_extend_assign, PSO_global_save, PSO_time(r)] = PSOfunction(model);
%         %PSO_assign = DecodeFromExtendPlan(model,PSO_extend_assign);
%         PSO_max(r) = PSO_global_save(end);
        
        
        numLimit = max(model.num_missiles, sum(model.target_require_num_list));
        s = zeros(20,numLimit); % 预分配内存
        for i=1:20
            s(i,:) = randperm(numLimit);
        end
        
        % SAGA算法
        [SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time(r)] = SAGAfunction(s,model,20,200,0.4,0.4, "PMX", "EM",100,0.95);
        % SAGA_assign = DecodeFromExtendPlan(model,SAGA_extend_plan);
        SAGA_max(r) = SAGA_global_save(end);
        
        % ACOVFSAGA
%         popsize=20;steps=200;Pcross=0.4;Pmutate=0.4;
%         Top = 50; Dop = 0.95;
%         [ACOVFSAGA_plan, ACO_global_save,ACO_mean_save, ACO_time(r)] = ACOVFSAGA(model,2*popsize,steps,Pcross,Pmutate,Top,Dop);
%         % ACO_plan = DecodeFromExtendPlan(model,ACOVFSAGA_plan);
%         ACO_max(r) = ACO_global_save(end);
%         
%         % GA
%         [GA_extend_plan,GA_global_save,GA_mean_save, GA_time(r)] = GA1(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
%         % GA_assign = DecodeFromExtendPlan(model,GA_extend_plan);
%         GA_max(r) = GA_global_save(end);
        
        %             [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
        %             % IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);
        %             IGA_max(r) = max(IGA_global_save);
        
        
        
        
    end
    %归一化
    maxU = max([JSFP_max,GRM_max,SAP_max,HCGSAP_max,HCGDMEA_max,SAGA_max]);
    JSFP_max = JSFP_max/maxU;
    GRM_max = GRM_max/maxU;
    SAP_max = SAP_max/maxU;
%     PSO_max = PSO_max/maxU;
%     ACO_max = ACO_max/maxU;
    SAGA_max = SAGA_max/maxU;
%     GA_max = GA_max/maxU;
    HCGSAP_max = HCGSAP_max/maxU;
    HCGDMEA_max = HCGDMEA_max/maxU;
    %             IGA_max(r) = IGA_max(r)/maxU;
    indexM = Nm/10;
    indexN = Nt/10;
    %         for i=1:l
    %             eval([algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max=',algo_cell{i},'_max']);
    %             eval([algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_time=',algo_cell{i},'_time']);
    %
    %             save([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_max'],[algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_max']);
    %             save([filename,algo_cell{i},'_',num2str(Nm),'vs',num2str(Nt),'_time'],[algo_cell{i},num2str(Nm),'vs',num2str(Nt),'_time']);
    %         end
    
end



%
% for i=1:l
%     save([filename,algo_cell{i},'_max'],[algo_cell{i},'_mean_time']);
%     save([filename,algo_cell{i},'_25tile_time'],[algo_cell{i},'_25tile_time']);
%     %     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
% end