function out = AlignToStimulus(data,SF)

lp = size(data,1);
lt = size(data,2);
ly = size(data,4);
lx = size(data,5);
lf = 20;

out = zeros(ly,lx,lp,lt*lf);

for p = 1:lp
    temp_out = zeros(ly,lx,lt*lf);
    for t = 1:lt
        stimF = SF(t);
        startF = stimF - 10;
        endF = startF + lf - 1;
        baseRange = startF:(stimF-1);
        range = startF:endF;
        for f = 1:lf
            curF = range(f);
            if curF < lf + 11
                i = (t-1)*lf + f;
                cur_data = squeeze(data(p,t,curF,:,:));
                temp_out(:,:,i) = cur_data;
            end
        end
        
    end

    out(:,:,p,:) = temp_out;
end

end