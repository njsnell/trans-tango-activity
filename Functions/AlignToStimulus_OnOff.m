function out = AlignToStimulus_OnOff(data,SF1,SF2)

lp = size(data,1);
lt = size(data,2);
ly = size(data,4);
lx = size(data,5);
lf = 25;

out = zeros(ly,lx,lp,lt*lf);

for p = 1:lp
    temp_out = NaN(ly,lx,lt*lf);
    for t = 1:lt
        onF = SF1(t); offF = SF2(t);
        onStart = onF - 10; onEnd = onF + 5;
        offStart = offF - 3; offEnd = offF + 5;
        r1 = onStart:onEnd; r2 = offStart:offEnd;
        out_r1 = (t-1)*25 + (1:16); out_r2 = (t-1)*25 + (17:25);
        for r = 1:16
            temp_out(:,:,out_r1(r)) = data(p,t,r1(r),:,:);
        end
        for r = 1:9
            if r2(r) < 41
                temp_out(:,:,out_r2(r)) = data(p,t,r2(r),:,:);
            end
        end
    end

    out(:,:,p,:) = temp_out;
end

end