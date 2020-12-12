% 穷举法求精确解(Exhaustive Attack method)
tic;
N = max(model.num_missiles, sum(model.target_require_num_list));
allPerms = perms(1:N);
numOfAll = size(allPerms,1);
fEAM = zeros(numOfAll,1);
for i=1:numOfAll
    plan = allPerms(i,:);
    fEAM(i) =  GlobalUtility(model, plan);
end
plot(fEAM,'o');
[maxF,maxIndex] = max(fEAM)
bestEAMPlan = allPerms(maxIndex,:) 
bestEAM_assign = DecodeFromExtendPlan(model,bestEAMPlan);
clear allPerms;
toc

%% 全局效用
function [Ug] = GlobalUtility(model,plan)
Nm = model.num_missiles;
Nt = model.num_targets;
Ut = zeros(Nt,1);

assign = zeros(model.num_missiles,1);
for i=1:Nm
    assign(i) = model.order_targets(plan(i));
end

for j=1:Nt
    part_missiles = find(assign == j);
    num_missiles = length(part_missiles);
    if num_missiles == 0
        Ut(j) = 0;
    else
        time_max = max(model.Time_to_go(part_missiles,j));
        J_sum = sum(model.Energy_opt(part_missiles,j));
        cost = time_max + J_sum;
        Ut(j) = max(0,model.Targets.value(j) - cost);
    end
    
end
Ug = sum(Ut);
end