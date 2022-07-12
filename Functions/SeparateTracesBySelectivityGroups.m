function c_out = SeparateTracesBySelectivityGroups(T_in,sigma,nStim,fpt,delays)

[C,~] = size(T_in);
c_out = zeros(C,1);
stimL = 3;
if isempty(delays)  % ON only
    on_off = false;
    nConditions = nStim;
else        % ON and OFF
    on_off = true;
    nConditions = nStim*2;
end

for c = 1:C
    curT = T_in(c,:);
    n = zeros(1,nConditions);
    for s = 1:nStim
        stimtype = (s-1)*3;
        if on_off
            for o = 1:2
                x = zeros(1,3);
                for r = 1:3
                    tr = stimtype + r;
                    if o == 1 % "on"
                        baseRange = (tr-1)*fpt + (1:10);
                        stimRange = (tr-1)*fpt + (12:(11+stimL));
                    else % "off"
                        baseRange = (tr-1)*fpt + delays(tr) + (6:10);
                        if delays(tr) + 11 + stimL > fpt
                            stimRange = (tr-1)*fpt + ((12+delays(tr)):fpt);
                        else
                            stimRange = (tr-1)*fpt + delays(tr) + (12:(11+stimL));
                        end
                    end
                    F0 = mean(curT(baseRange));
                    std0 = std(curT(baseRange));
                    F = mean(curT(stimRange));
                    Z = (F-F0)/std0;
                    if Z > sigma && ~isinf(Z) 
                        x(r) = F;
                    end
                end
                if sum(x>0) >= 2
                    n((s-1)*2+o) = 1;
                end
            end
        else
            x = zeros(1,3);
            for r = 1:3
                tr = stimtype + r;
                baseRange = (tr-1)*fpt + (1:10);
                stimRange = (tr-1)*fpt + (12:(11+stimL));
                F0 = mean(curT(baseRange));
                std0 = std(curT(baseRange));
                F = mean(curT(stimRange));
                Z = (F-F0)/std0;
                if Z > sigma && ~isinf(Z) 
                    x(r) = F;
                end
            end
            if sum(x>0) >= 2
                n(s) = 1;
            end
        end
    end
    
    m_base = 2*ones(1,nConditions);
    m_exp = 0:(nConditions-1);
    m = m_base .^ m_exp;
    multiplier_raw = m(end:-1:1)';
    multiplier = multiplier_raw((end+1-(nConditions)):end);
    c_out(c) = n*multiplier; 
    
end

end