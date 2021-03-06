function plotMe(out)

warning('off', 'all')
pred_val = out.measures.mu.RH;
data_val = out.measures.mu.HD;

xref = out.model.fitParams.xref;
yref = out.model.fitParams.yref;

% reshape MVL (scaling factor)
data_MVL = out.measures.MVL.HD;
data_MVL(isnan(data_MVL))=0;

model_MVL = out.measures.MVL.RH;
model_MVL(isnan(model_MVL))=0;

% todo: normalize the scaling
% (data scale tends to be much higher)

% define vector orientations
u = cos(pred_val * pi/180); 
v = sin(pred_val * pi/180);
u_data = cos(data_val * pi/180);
v_data = sin(data_val * pi/180);

% scaling factors
sf = model_MVL;
sf_data = data_MVL;

% scale the vectors
uprime = u.*sf; 
vprime = v.*sf;
uprime_data = u_data.*sf_data; 
vprime_data = v_data.*sf_data;


% get rid of nan values in ratemap
rm_vec = reshape(out.data.rxy',100,1);
nan_idx = find(isnan(rm_vec));

% grab bin information
binX = out.info.bin.X;
binY = out.info.bin.Y;

%% PLOT

if sum(~isnan(out.data.rxyN), 'all') > 1
    figure; set(gcf,'color','w');
    hold on;
    set(gca, 'visible', 'off')
    imagescwithnan(out.data.rxyN,jet,[1 1 1])% colorbar
%     [.7 .5 .7]
    alpha(0.2) 
    brighten(.6)
    xlim([0 11]); ylim([0 11]);
    
    xline(0, 'LineWidth', 1, 'Color', [.8 .8 .8]);
    xline(11, 'LineWidth', 1, 'Color', [.8 .8 .8]);
    yline(0, 'LineWidth', 1, 'Color', [.8 .8 .8]);
    yline(11, 'LineWidth', 1, 'Color', [.8 .8 .8]);
    
    % plot data vectors
    modelVecs = quiver(binX, binY, uprime_data, vprime_data, 0);
    set(modelVecs, 'Color', 'b', 'AutoScale', 'off', 'LineWidth',1)

    % plot model vectors
    modelVecs = quiver(binX, binY, uprime, vprime, 0);
    set(modelVecs, 'Color', 'r', 'AutoScale', 'off', 'LineWidth',.35)

    % put in the reference point if its not distance
    if xref>0 & xref<150 & yref>0 & yref<150
        scatter(xref,yref,[10],'k','filled')
    end
else
    [r,c] = size(out.data.rxyN);
    imagesc(zeros(r,c));
end
    pbaspect([1 1 1])
    warning('on', 'all')

end
