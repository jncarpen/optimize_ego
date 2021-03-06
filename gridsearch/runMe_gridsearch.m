%% RH model reference point grid search 
% you need -->
% P: position [t x(led1) y(led1) x(led2) y(led2)]
% ST: spiketimes (seconds)
% Z: angular variable (degrees)

% unpackage variables
posx = P(:,2); posy = P(:,3);
posx2 = P(:,4); posy2 = P(:,5);
posx_c = (posx + posx2)./2; posy_c = (posy + posy2)./2;
post = P(:,1); % timestamps
ST = ST(ST < post(end) & ST > post(1));
spiketrain = histcounts(ST, linspace(post(1),post(end),numel(post)+1))';

% initialize number of reference points to test
num_xy_bins = 25; % inside of box
[xref_matrix, yref_matrix] = meshgrid(linspace(-1,10,num_xy_bins),...
    linspace(-1,10,num_xy_bins));
% linearize the meshgrid
xref_vector = reshape(xref_matrix,num_xy_bins^2,1);
yref_vector = reshape(yref_matrix,num_xy_bins^2,1);
numRuns = length(xref_vector);

%% Iterate over each reference point

for refPointIter = 1:numRuns
    % reference point now
    fprintf('Testing reference point %d of %d\n', refPointIter, length(xref_vector));
    ref = [xref_vector(refPointIter), yref_vector(refPointIter)];
    [rhOut{refPointIter}] = modelMe_gridsearch(P, ST, Z, ref);
    error(refPointIter) = rhOut{refPointIter}.model.error;
end

error2 = reshape(error',sqrt(numRuns),sqrt(numRuns));

%% plot results
% initialize colormap
figure; set(gcf,'color','w');
my_colormap = hot; my_colormap = my_colormap(1:end-10,:);
surf(error2); 
colormap(my_colormap); cb=colorbar; cb.FontSize = 12;









