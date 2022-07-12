function varargout = SelectROIs(SC,traces,cList,ap,max_c)

% make c_in all of them

global ww; global tt; global pushed; global mean_trace;
pushed = false;
mean_trace = NaN(1,size(traces,2));
colors = hsv(max_c);
d1 = 128; d2 = 128; d3 = 12;
thr = 0.95;
alpha = 0.8;
bg = zeros(128,128);
cen = com(SC,d1,d2,d3);

if strcmp(ap,'a')
    idx = find(cen(:,3) < 6.5);  
elseif strcmp(ap,'p')
    idx = find(cen(:,3) >= 6.5);
else
    error('Choose a or p\n');
end

cen2d = [];
clall = [];
A_temp = zeros(d1*d2*d3,max_c);
for c = 1:max_c
    c_idx = find(cList == c);
    cl = c_idx(ismember(c_idx,idx));
    clall = [clall; cl];
    cen2d = [cen2d; cen(cl,1:2)];
    A_build = zeros(d1*d2*d3,1);
    for i = cl'
        A_raw = full(SC(:,i));
        A_norm = ceil(mat2gray(A_raw(:)));
        A_build(A_norm == 1) = 1;
    end
    A_temp(:,c) = A_build;
end
A_plot = A_temp * colors;
A_plot = reshape(A_plot,d1,d2,d3,3);
A_proj = squeeze(max(A_plot,[],3));

fig1 = figure;
set(fig1,'Position',[80,80,1200,500]);

subplot(2,2,[1 3]);
imshow(bg,[]); 
hold on;
alphamap = alpha .* ceil(max(A_proj,[],3));
image(A_proj,'AlphaData',alphamap);  

wwall = cell(length(clall),1); 
roison = zeros(length(clall),1);
for i = 1:length(clall)
    curROI = reshape(full(SC(:,clall(i))),d1,d2,d3);

    curROI = squeeze(max(curROI,[],3));
    curROI = medfilt2(curROI,[12,12]);
    [temp,ind] = sort(curROI(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first'); 
    subplot(2,2,[1 3]);
    [~,ww] = contour(reshape(curROI,d1,d2),[0,0]+curROI(ind(ff)));
    ww.LineColor = [0.5 0.5 0.5]; ww.LineWidth = 2;
    wwall{i} = ww;
end

fig2 = uifigure;
set(fig2,'Position',[30,200,300,200]);
btn = uibutton(fig2,'Text','Done!','ButtonPushedFcn',@(btn,event) returnTTA(btn,event));

while ~pushed
    chooseROI();
end

close(fig1); close(fig2);

    function returnTTA(btn,evt)
        pushed = true;
        cl_out = find(roison);

        out_temp = zeros(d1*d2*d3,length(cl_out));
        out_traces = traces(clall(cl_out),:);
        for j = 1:length(cl_out)
            cur = cl_out(j);
            out_raw = full(SC(:,clall(cur)));
            out_temp(:,j) = mat2gray(out_raw(:));
        end
        out_plot = out_temp * colors(cList(clall(cl_out)),:);
        out_plot = reshape(out_plot,d1,d2,d3,3);
        out_map = squeeze(max(out_plot,[],3));
        varargout{1} = out_map;
        varargout{2} = out_traces;
    end

    function plotROI(i)
        plotTrace(i);
        if roison(i)
            roison(i) = 0;
            wwall{i}.LineColor = [0.5 0.5 0.5];
        else
            roison(i) = 1;
            wwall{i}.LineColor = [1 1 1];
        end
        plotMeanTrace();
    end

    function plotTrace(i)
        delete(tt);
        A_trace = traces(clall(i),:); 
        figure(fig1);
        subplot(2,2,2);
        tt = plot(A_trace,'b');
    end

    function chooseROI()
        [x,y] = ginput(1);
        try
            nearest = dsearchn(cen2d,[y x]);
            plotROI(nearest);
        catch
        end
    end
end