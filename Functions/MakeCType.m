function [tune,resp] = MakeCType(nStim)

combos = 2 ^ nStim;
tune = cell(nStim+1,1);
resp = cell(nStim,1);
for x = 0:combos-1
    try
        b = de2bi(x,nStim,'left-msb');
    catch
        b_siz = nStim; %ceil(log(nStim+1)/log(2));
        b = NaN(1,b_siz);
        t = dec2bin(x,b_siz);
        for i = 1:length(t)
            b(i) = str2num(t(i));
        end
    end
    b_idx = sum(b)+1; %sum(b==1)
    tune{b_idx} = [tune{b_idx} x];
    for y = find(b==1)
        resp{y} = [resp{y} x];
    end
end

end