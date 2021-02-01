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
% figure; set(gcf,'color','w');
hold on;
% set(gca, 'visible', 'off')
imagescwithnan(out.data.rxy,jet,[.7 .5 .7])% colorbar
alpha(0.2) 
brighten(.6)
xlim([0 11]); ylim([0 11]);

% plot data vectors
modelVecs = quiver(binX, binY, uprime_data, vprime_data, 0);
set(modelVecs, 'Color', 'blue', 'AutoScale', 'off', 'LineWidth',.5)

% plot model vectors
modelVecs = quiver(binX, binY, uprime, vprime, 0);
set(modelVecs, 'Color', 'r', 'AutoScale', 'off', 'LineWidth',.5)

% put in the reference point if its not distance
if xref>0 & xref<150 & yref>0 & yref<150
    scatter(xref,yref,[10],'k','filled')
end

% make the figure a square
pbaspect([1 1 1])
warning('on', 'all')

end
