function [out,toRemove] = NormalizeTraces(in,nStim,fpt)

[C,t] = size(in);
out = zeros(C,t);

for c = 1:C
    for tr = 1:(nStim*3)
        baseR = (tr-1)*fpt + (1:10);
        fullR = (tr-1)*fpt + (1:fpt);
        base = mean(in(c,baseR));
        out(c,fullR) = in(c,fullR) ./ base;
    end
end

toRemove = unique(...
    [find(any(out==0,2)); find(any(isnan(out),2))]);
out(toRemove,:) = [];

end