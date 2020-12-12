clear
clc

%% 导弹与目标位置初始化
% Nm = 100;
% Nt = 100;

for Nm = 10:10:100
    for Nt=10:10:100
        model = DynamicMissileAndTarget(Nm,Nt);
        model = RandomSetSituation(model);
        
        
        %% 分配
        Epis = 500;
        
        Times = 100;
        JSFP_time = zeros(1,Times);
        GRM_time = zeros(1,Times);
        SAP_time = zeros(1,Times);
        PSO_time = zeros(1,Times);
        ACO_time = zeros(1,Times);
        SAGA_time = zeros(1,Times);
        GA_time = zeros(1,Times);
        IGA_time = zeros(1,Times);
        for r = 1:Times
            %JSFP算法
            [JSFP_assign, JSFP_global_save, JSFP_time] = JSFPfunction(model,Epis);
            
            % GRMFMI算法
            [GRM_assign, GRM_global_save, GRM_time(r)] = GRMfunction(model,Epis);
            
            % SAP算法
            tic
            [SAP_lagrange2_assign, SAP_global_save,SAP_time(r)] = SAP_lagrange_function2(model,Epis);
            % SAP_max(r) = max(SAP_global_save);
            
            % PSO算法
            tic
            [PSO_extend_assign, PSO_global_save, PSO_time(r)] = PSOfunction(model);
            %PSO_assign = DecodeFromExtendPlan(model,PSO_extend_assign);
            % PSO_max(r) = max(PSO_global_save);
            
            
            numLimit = max(model.num_missiles, sum(model.target_require_num_list));
            s = zeros(20,numLimit); % 预分配内存
            for i=1:20
                s(i,:) = randperm(numLimit);
            end
            
            % SAGA算法
            [SAGA_extend_plan,s,SAGA_global_save,SAGA_mean_save, SAGA_time(r)] = SAGAfunction(s,model,20,200,0.4,0.4, "PMX", "EM",100,0.95);
            % SAGA_assign = DecodeFromExtendPlan(model,SAGA_extend_plan);
            %SAGA_max(r) = max(SAGA_global_save);
            
            % ACOVFSAGA
            popsize=20;steps=200;Pcross=0.4;Pmutate=0.4;
            Top = 50; Dop = 0.95;
            [ACOVFSAGA_plan, ACO_global_save,ACO_mean_save, ACO_time(r)] = ACOVFSAGA(model,2*popsize,steps,Pcross,Pmutate,Top,Dop);
            % ACO_plan = DecodeFromExtendPlan(model,ACOVFSAGA_plan);
            % ACO_max(r) = max(ACO_global_save);
            
            % GA
            [GA_extend_plan,GA_global_save,GA_mean_save, GA_time(r)] = GA1(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
            % GA_assign = DecodeFromExtendPlan(model,GA_extend_plan);
            % GA_max(r) = max(GA_global_save);
            
            [IGA_extend_plan,IGA_global_save,IGA_mean_save, IGA_time] = IGA(model,popsize,steps,Pcross,Pmutate, "PMX", "SIM");
            % IGA_assign = DecodeFromExtendPlan(model,IGA_extend_plan);
            IGA_max(r) = max(IGA_global_save);
            
            % 归一化
            % maxU = max([SAP_max(r),PSO_max(r),ACO_max(r),SAGA_max(r),GA_max(r),IGA_max(r)]);
            %
            % SAP_max(r) = SAP_max(r)/maxU;
            % PSO_max(r) = PSO_max(r)/maxU;
            % ACO_max(r) = ACO_max(r)/maxU;
            % SAGA_max(r) = SAGA_max(r)/maxU;
            % GA_max(r) = GA_max(r)/maxU;
            % IGA_max(r) = IGA_max(r)/maxU;
            
            
        end
        index = Nm/10;
        JSFP_mean_time(index) = mean(JSFP_time);
        GRM_mean_time(index) = mean(GRM_time);
        SAP_mean_time(index) = mean(SAP_time);
        PSO_mean_time(index) = mean(PSO_time);
        ACO_mean_time(index) = mean(ACO_time);
        SAGA_mean_time(index) = mean(SAGA_time);
        GA_mean_time(index) = mean(GA_time);
        IGA_mean_time(index) = mean(IGA_time);
    end
end

algo_cell = {'JSFP';'GRM';'SAP';'ACO';'PSO';'GA';'IGA';'SAGA'};
[l,~]=size(algo_cell);
filename = 'TimeTest/';
for i=1:l
    save([filename,algo_cell{i},'_time'],[algo_cell{i},'_mean_time']);
    %     save([filename,algo_cell{i},'_time'],[algo_cell{i},'_time']);
end