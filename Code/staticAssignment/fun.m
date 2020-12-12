function [ fitness ] = fun(model, plan)
fitness=0;
[~,n]=size(individual_pop);
for i=1:n
    fitness=fitness+Val(i,individual_pop(1,i));
end
end

