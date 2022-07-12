function varargout = AutoselectROIs(SC,cList,c,nSect,colors)
% Copy to work on multiple c-types
varargout = cell(nSect,1);
% colors = hsv(nC); colors = colors([end 1:(end-1)],:);
d1 = 128; d2 = 128; d3 = 12;
%nSect = 2;
step = (d3+1)/nSect;
cen = com(SC,d1,d2,d3);
for s = 1:nSect
    cLo = (s-1)*step; cHi = s*step;
    roi_hiBound = find(cen(:,3) <= cHi); 
    roi_loBound = find(cen(:,3) > cLo);
    roi_cur = intersect(roi_hiBound,roi_loBound);
    cl_cur = find(cList == c);
    cur = intersect(cl_cur,roi_cur);

    out_temp = zeros(d1*d2*d3,1);
    for i = 1:length(cur)
        out_raw = full(SC(:,cur(i)));
        out_norm = mat2gray(out_raw(:));
        out_temp = out_temp + out_norm;
    end
    out_plot = out_temp * colors(c,:);
    out_plot = reshape(out_plot,d1,d2,d3,3);
    out_map = squeeze(max(out_plot,[],3));
    varargout{s} = out_map;
end    

end