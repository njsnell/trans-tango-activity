function output = ComputeMeanTrace(input,fpt)

nF = length(input);
nStim = nF/(3*fpt);
output = zeros(1,nF/3);
for s = 1:nStim
    t1 = input(((s-1)*3)*fpt+(1:fpt));
    t2 = input(((s-1)*3+1)*fpt+(1:fpt));
    t3 = input(((s-1)*3+2)*fpt+(1:fpt));     
    output((s-1)*fpt+(1:fpt)) = mean([t1 t2 t3],2);
end
    
end