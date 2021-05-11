function [out] = modelMe(P, ST, Z)
%MODELME_V2 
%   INPUTS -
%   P:         position vector [t, x1, y1, x2, y2] or [t, x, y];
%   ST:        spike times vector (s)
%   Z:         angular variable, range 0 to 360 deg.
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position info
t = P(:,1); 
x = P(:,2);
y = P(:,3);
Z = mod((Z+180),360)-180; % shift + mirror

% sampling frequency info
tpf = mode(diff(t)); % time per frame (s)

% remove spikes outside of viable range
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);

% raw spike train
t_edges = linspace(startTime,stopTime,numel(t)+1)';
SpkTrn = histcounts(ST,t_edges)';

% bin arena (10x10)
nBins = 10;
[~, xEdges, yEdges, binX, binY] = histcounts2(x,y,nBins);
xCenter = (diff(xEdges)/2) + xEdges(1:end-1);
yCenter = (diff(yEdges)/2) + yEdges(1:end-1);

% bin centers vector
count = 1;
for cc = 1:length(xCenter)
    for rr = 1:length(yCenter)
        binCenters(count,1:2) = [xCenter(rr), yCenter(cc)];
        count = count+1;
    end
end

% angular bin centers
num_Z_bins = 10; % 36 deg/bin
Z_bins = linspace(0,360,num_Z_bins+1);
Z_bin_ctrs = ((diff(Z_bins)/2) + Z_bins(1:end-1))-180;


%% GENERATE RATEMAPS
r_xy = zeros(nBins,num_Z_bins).*NaN;
r_xyh = zeros(nBins,nBins, num_Z_bins).*NaN;
R_xyh = zeros(nBins,nBins, num_Z_bins).*NaN;
time_H = zeros(nBins,nBins, num_Z_bins).*NaN;
count_H = zeros(nBins,nBins, num_Z_bins).*NaN;
Occ = cell(10,10);

count = 1;
for rr = 1:nBins
    for cc = 1:nBins
        % what bin are we in now?
        x_bin_here(rr,cc) = cc; y_bin_here(rr,cc) = rr;
        
        % coordinates (in cm) of the location of this particular bin
        loc3d(rr,cc,:) = [binCenters(count,1), binCenters(count,2)];
        
        % find frames in which animal occupied this spatial bin
        idx_here = find(rr == binY & cc == binX);
        time_in_bin = length(idx_here)*tpf; % occupancy (s)
        
        % spikes and angular variable Z (rad) in this spatial bin
        spikes_here = SpkTrn(idx_here); 
        Z_here = Z(idx_here);
        
        % compute occupany in each angular bin
        % linear histogram because input constrained between -pi and pi
        Z_edges = linspace(-180, 180, nBins+1);
        [Z_count_here, ~, Z_idx_here] = histcounts(Z_here, Z_edges);
        Z_idx_here(Z_idx_here==0) = NaN; 
        
        % calculate angular occupancy (s) for this spatial bin
        Z_occ_here = Z_count_here .* tpf; 
        
        % @criteria: animal must have visited 50% of bins
        bin_threshold = 0.4; % 400 ms/each
        bin_criteria = round(num_Z_bins*.5);
        num_bins_passed = sum(Z_occ_here>bin_threshold);
        
        if num_bins_passed >= bin_criteria
            % loop through each angular (Z) bin
            for H = 1:length(Z_bin_ctrs)
                % amount of time animal spent in this HD bin (s)
                time_H(rr,cc,H) = Z_occ_here(H);
                count_H(rr,cc,H) = Z_count_here(H);

                % find indices when animal occupied this HD bin
                idx_H = find(Z_idx_here == H);

                % @criteria: animal must have occupied each angular bin for >=500 ms
                if time_H(rr,cc,H) > .5
                    % spiketimes in bin(x,y,H)
                    spk_H = spikes_here(idx_H);

                    % conditional rate, r(x,y,H)
                    r_xyh_here = sum(spk_H)./(count_H(rr,cc,H)*tpf);

                    if isfinite(r_xyh_here)
                        r_xyh(rr,cc,H) = r_xyh_here;
                    end
                end
            end
            % make averaged ratemap
            r_xy(rr,cc) = nanmean(r_xyh(rr,cc,:));
            % conditional (normalized) ratemap
            R_xyh(rr,cc,:) = r_xyh(rr,cc,:)./r_xy(rr,cc);
            % save occupancy map
            Occ{rr,cc} = Z_count_here;
        else
            Occ{rr,cc} = nan;
        end
        count = count + 1; 
    end 
end

%% OPTIMIZATION
% randomly choose some initial conditions
p = choose_initial_conditions(nBins);

% or choose initial conditions manually:
% p.g = 0.25;
% p.thetaP = 0;
% p.xref = nanmax(x)/2;
% p.yref = nanmax(y)/2;

% min firing rate for a bin to be considered
rCutOff = .5;

%initial conditions of the 4 parameters to fit
pFitLM = [p.g; p.thetaP; p.xref; p.yref]; 

% options for fminsearch
options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8);

% perform the optimization
[pFit, ~]=fminsearch(@(pFit)aFitLMNew(pFit,r_xyh,r_xy,rCutOff,nBins),...
    pFitLM,options);

% get model-predicted firing rates for best-fit parameters
[R_xyh_model, fF] = get_Rxyh_model(pFit,r_xyh,r_xy,rCutOff,nBins);

% match nan values across model & data
[RDatanan, rDatanan, rxyDatanan] = matchnans(R_xyh_model, R_xyh, r_xyh, r_xy);


%% VARIANCE EXPLAINED BY MODEL (RH-TUNING)
% variance in r(x,y,H); excluding nans
% mean_Data = mean(rDatanan, 'all', 'omitnan');
% count = 0; 
% fF_Data = 0;
% for rr=1:nBins
%         for cc=1:nBins
%             tc_now = squeeze(rDatanan(rr, cc, :));
%             finiteBins = find(isfinite(tc_now));
%             if ~isempty(finiteBins)
%                 fF_Data = fF_Data + nansum((tc_now(finiteBins) - mean_Data).^2);
%                 count = count + length(finiteBins);
%             end
%         end
% end

% var_Data = fF_Data./count;
var_sub = nanvar(rDatanan(:)-R_xyh_model(:));
var_Data = nanvar(rDatanan(:));
% this line would do basically the same thing:
% var(rDatanan(find(isfinite(rDatanan))))
%  var(modelRateMap - rDatanan(find(isfinite(rDatanan))))
% variance explained by model
VEM = 1-(var_sub./var_Data);
% VEM = 1-(fF./var_Data);

%% VARIANCE EXPLAINED BY PLACE TUNING
 VEP_fF = 0;
    VEP_count = 0;
     for rr=1:nBins
        for cc=1:nBins
            % grab the conditional ratemap now
             tc_now = squeeze(rDatanan(rr, cc, :));
             % find indices of finite bins
             finiteBins = find(isfinite(tc_now));
             if ~isempty(finiteBins)
                 % mean squared error at each bin
                  VEP_fF = VEP_fF + nansum((tc_now(finiteBins)...
                      - rxyDatanan(rr,cc)).^2); 
                  VEP_count = VEP_count + length(finiteBins);
             end
        end
     end
     % subtract using repmat
VEP_fF = VEP_fF/VEP_count;
VEP = 1-(VEP_fF./var_Data);

%% MODULATION STRENGTH
warning('off','all')
for rr = 1:nBins
    for cc = 1:nBins
        % grab angular tuning curve in each spatial bin
        tc_now = reshape(rDatanan(rr,cc,:), nBins, 1);
        tc_now_RH = reshape(R_xyh_model(rr,cc,:), nBins, 1);
        
        % which bins are finite (~nan, ~inf)
        finite_bins = find(isfinite(tc_now));
        
        % grab these bins
        tc_now = tc_now(finite_bins);
        tc_now_RH = tc_now_RH(finite_bins);
        bins_now = deg2rad(Z_bin_ctrs(finite_bins))';
        
        % take circular mean (in RADIANS)
        % note: in jercog paper they sum (which would just be mu_rad * 10)
        [mu_raw, ~, ~] = circ_mean(bins_now, tc_now);
        mu(rr,cc) = rad2deg(mu_raw);
        
        [mu_RH_raw, ~, ~] = circ_mean(bins_now, tc_now_RH);
        mu_RH(rr,cc) = rad2deg(mu_RH_raw);
        
        % mean vector length
        MVL(rr,cc) = circ_r(bins_now, tc_now);
        MVL_RH(rr,cc) = circ_r(bins_now, tc_now_RH);

        
    end
end
warning('on','all')

% linear average of tuning strengths
MVL(MVL==Inf) = NaN; 
tuningStrength_HD = mean(reshape(MVL, nBins.^2,1), 'all', 'omitnan');

MVL_RH(MVL_RH==Inf) = NaN;
tuningStrength_RH = mean(reshape(MVL_RH, nBins^2,1), 'all', 'omitnan');


%% PREPARE OUTPUTS
% MODEL CLASS
out.model.Rxyh = R_xyh_model;
out.model.error = fF;
out.model.fitParams.g = pFit(1);
% out.model.fitParams.thetaP = mod(pFit(2),360)-180;
out.model.fitParams.thetaP = pFit(2);
out.model.fitParams.xref = pFit(3);
out.model.fitParams.yref = pFit(4);

% DATA CLASS
out.data.Rxyh = R_xyh;
out.data.rxyh = r_xyh;
out.data.rxy = r_xy;
out.data.RxyhN = RDatanan;
out.data.rxyhN = rDatanan;
out.data.rxyN = rxyDatanan;
out.data.occ.count = count_H;
out.data.occ.countmat = Occ;
out.data.occ.time = time_H;

% MEASURES
out.measures.VE.place = VEP;
out.measures.VE.RH = VEM;
out.measures.MVL.HD = MVL;
out.measures.MVL.RH = MVL_RH;
out.measures.TS.HD = tuningStrength_HD;
out.measures.TS.RH = tuningStrength_RH;
out.measures.mu.HD = mu;
out.measures.mu.RH = mu_RH;

% INFO
out.info.bin.X = x_bin_here;
out.info.bin.Y = y_bin_here;

%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [RXYHM, fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins)
    %   INPUTS--
    %   'rF'            R_data(x,y,H) (Nbins x Nbins x Nbins)
    %   'rP'            r(x,y) (Nbins x Nbins)
    %   'rCutOff'       min firing rate for a bin to be considered
    %   'Nbins'         number of bins of the discretization
    %   'pFit'          [g; thetaP; xref; yref]; initial conditions of params.
    %   OUTPUTS--
    %   'RXYHM'         model fit (Nbins x Nbins x Nbins)
    %   'fF'            mean squared error between model's predicted
    %                   firing rates & data [rF, R(x,y,H)].
    %   by P. Jercog., notes/slight modification by me
    %   I added the RXYHM output (where the model's predicted rates are saved).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parse pFit input
    g = pFit(1);
    thetaP = pFit(2); % (deg)
    xref = pFit(3);
    yref = pFit(4);
    % initialize loop
    RXYHM = zeros(Nbins,Nbins,Nbins).*nan;
    fF = 0;
        angCount = 0;
        for ii=1:Nbins
            for jj=1:Nbins
                % average FR in spatial bin must exceed cutoff
                if (rP(ii, jj)>rCutOff)
                    % angular tc for spatialbin(ii,jj)
                    rT = squeeze(rF(ii, jj, :));
                    % indices of finite values (ignore nan/inf)
                    iF = find(isfinite(rT));
                    if ~isempty(iF)            
                        a = 180*atan2(yref-ii, xref-jj)/pi - ...
                            (-180 -360/(2*Nbins) + iF*360/Nbins);
                        cFac = cos(pi*(a-thetaP)/180);
                        % expected value
                        cBar = nanmean(cFac);
                        % shape the cosine
                        z = 1+g.*(cFac - cBar);
                        z = z.*(z>0);
                        % variance explained/SSE
                        fF = fF + nansum((z - rT(iF)).^2); 
                        angCount = angCount + length(iF);
                        % model firing rates
                        RXYHM(ii,jj,iF) = z;
                    end
                end
            end
        end
        % average SSE over *finite* bins
        fF = fF/angCount; 
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [fF] = aFitLMNew(pFit,rF,rP,rCutOff,Nbins)
        % call the other function (optimization)
        [~, fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function initial = choose_initial_conditions(total_bins)
        % make a vector of all possible positions
        x_bins = 1:5.:total_bins;
        y_bins = 1:.5:total_bins;
        orientation = -180:1:180;
        % randomly sample from them 
        howMany = 1;
        initial.g = rand(howMany, 1).*2.5;
        initial.thetaP = randsample(orientation,howMany)';
        initial.xref = randsample(x_bins,howMany)';
        initial.yref = randsample(y_bins,howMany)';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [RDatanan, rDatanan, rxyDatanan] = matchnans(RModel, RData, rData, rxyData)
        %MATCHNANS
        [r,c,v] = size(RModel);
        RModel = reshape(RModel, r*c*v, 1);
        RData = reshape(RData, r*c*v, 1);
        rData = reshape(rData, r*c*v, 1);
        nanidx = find(isnan(RModel));
        RData(nanidx) = NaN;
        rData(nanidx) = NaN;
        RDatanan = reshape(RData,r,c,v);
        rDatanan = reshape(rData,r,c,v);
        
        for row = 1:10
            for col = 1:10
                now = squeeze(rDatanan(row,col,:));
                fin = find(isfinite(now));
                rxyDatanan(row,col) = nansum(now(fin))./length(fin);
            end
        end
    end

end

