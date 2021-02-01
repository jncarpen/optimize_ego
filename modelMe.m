function [out] = modelMe(P, ST, Z)
%MODELME_V2 
%   INPUTS -
%   P:         position vector [t, x1, y1, x2, y2];
%   ST:        spike times vector (s)
%   Z:         angular variable, range (-180 to 180)
%   J. Carpenter, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position info
t = P(:,1); 
y = P(:,3); x = P(:,2);

% sampling frequency info
tpf = mode(diff(t)); % time per frame (s)
% fps = 1/tpf; % frames/sec (Hz)

% remove spikes outside of viable range
startTime = t(1); stopTime = t(end);
ST = ST(ST < stopTime & ST > startTime);

% raw spike train
t_edges = linspace(startTime,stopTime,numel(t)+1);
SpkTrn = histcounts(ST,t_edges);

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
r_xy = zeros(nBins,nBins).*NaN;
r_xyh = zeros(nBins,nBins, nBins).*NaN;
R_xyh = zeros(nBins,nBins, nBins).*NaN;

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
                      
        % average rate for this spatial bin (Hz)
        r_xy_here = sum(spikes_here)/(length(idx_here)*tpf);
        
        % make spatial ratemap
        % @criteria: animal must have occupied each 2D spatial bin 
        % for a total of >= 1000 ms
        if time_in_bin >= 1
            r_xy(rr,cc) = r_xy_here; 
        end
        
        % compute occupany in each angular bin
        % linear histogram because input constrained between -pi and pi
        Z_edges = linspace(-180, 180, nBins+1);
        [Z_count_here, ~, Z_idx_here] = histcounts(Z_here, Z_edges);
        Z_idx_here(Z_idx_here==0) = NaN; 
        
        % calculate angular occupancy (s) for this spatial bin
        Z_occ_here = Z_count_here .* tpf; 
        
        % loop through each angular (Z) bin
        for H = 1:length(Z_bin_ctrs)
            % amount of time animal spent in this HD bin (s)
            time_H = Z_occ_here(H);
            count_H = Z_count_here(H);
            
            % find indices when animal occupied this HD bin
            idx_H = find(Z_idx_here == H);
                
            % @criteria: animal must have occupied each angular bin for >=100 ms
            if time_H >= .1 
                % spiketimes in bin(x,y,H)
                spk_H = spikes_here(idx_H);

                % conditional rate, r(x,y,H)
                r_xyh_here = sum(spk_H)./(count_H*tpf);
                test(count) = r_xyh_here;
                
                if isfinite(r_xyh_here)
                    r_xyh(rr,cc,H) = r_xyh_here;
                else
                    r_xyh(rr,cc,H) = NaN;
                end

                % conditional (normalized) rate, R(x,y,H)
                if r_xy_here == 0 && r_xyh_here == 0
                    % since occupancy criteria is met by this step,
                    % rate map should be 0 (not NaN- which you would
                    % get if r_xy_here = r_xyh_here = 0 and 0/0 = NaN).
                    R_xyh(rr,cc,H) = 0;
                else
                    % normalize by average rate in spatial bin
                    R_xyh(rr,cc,H) = r_xyh_here./r_xy_here;
                end
            end
        end  
    end
    count = count + 1; 
end

% smooth ratemaps (?)
r_xy_S = r_xy;
r_xyh_S = r_xyh;
R_xyh_S = R_xyh;

% r_xy_S = smooth2a(r_xy,3,3);
% for i=1:10
%     r_xyh_S(:,:,i) = smooth2a(squeeze(r_xyh(:,:,i)),3,3);
%     R_xyh_S(:,:,i) = smooth2a(squeeze(R_xyh(:,:,i)),3,3);
% end
% for i=1:10
%     for j=1:10
%         r_xyh_S(i,j,:)=smoothdata(squeeze(r_xyh_S(i,j,:)),3,'omitnan');
%         R_xyh_S(i,j,:)=smoothdata(squeeze(R_xyh_S(i,j,:)),3,'omitnan');
%     end
% end

% how correlated are the two maps?
% corrcoef(squeeze(nanmean(r_xyh,3)),r_xy,'rows','complete')


%% OPTIMIZATION
% randomly choose some initial conditions
p = choose_initial_conditions(nBins);

% min firing rate for a bin to be considered
rCutOff = .5;

%initial conditions of the 4 parameters to fit
pFitLM = [p.g; p.thetaP; p.xref; p.yref]; 

% options for fminsearch
options = optimset('Display','off','TolX',1e-8,'TolFun',1e-8);

% perform the optimization
[pFit, ~]=fminsearch(@(pFit)aFitLMNew(pFit,R_xyh_S,r_xy_S,rCutOff,nBins),...
    pFitLM,options);

% get model-predicted firing rates for best-fit parameters
[R_xyh_model, VEM_num] = get_Rxyh_model(pFit,R_xyh_S,r_xy_S,rCutOff,nBins);

% variance in r(x,y,H); excluding nans
mean_rxyh = nanmean(r_xyh_S, 'all');
linear_rxyh = reshape(r_xyh_S, 1000,1);
VE_rxyh_count = 0; Var_rxyh = 0;
for binbin = 1:length(linear_rxyh)
    bin_now = linear_rxyh(binbin);
    if ~isnan(bin_now)
        Var_rxyh = Var_rxyh + nansum((bin_now - mean_rxyh)^2);
        VE_rxyh_count = VE_rxyh_count+1;
    end
end
VEM_den = Var_rxyh./VE_rxyh_count;
% variance explained by model
VEM = 1-(VEM_num./VEM_den);


%% VARIANCE EXPLAINED BY PLACE TUNING
 VEP_fF = 0;
    VEP_count = 0;
     for rr=1:nBins
        for cc=1:nBins
            % grab the conditional ratemap now
             crm_now = squeeze(r_xyh_S(rr, cc, :));
             % find indices of finite bins
             crm_if = find(isfinite(crm_now));
             if ~isempty(crm_if)
                 % mean squared error at each bin
                  VEP_fF = VEP_fF + nansum((crm_now(crm_if) - r_xy_S(rr,cc)).^2); 
                  VEP_count = VEP_count + length(crm_if);
             end
        end
     end
VEP_fF = VEP_fF/VEP_count;
VEP = 1-(VEP_fF./VEM_den);


%% MODULATION STRENGTH
warning('off','all')
for rr = 1:nBins
    for cc = 1:nBins
        % grab angular tuning curve in each spatial bin
        tc_now = reshape(R_xyh_S(rr,cc,:), nBins, 1);
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
out.model.error = VEM_num;
out.model.fitParams.g = pFit(1);
out.model.fitParams.thetaP = mod(pFit(2),360)-180;
out.model.fitParams.xref = pFit(3);
out.model.fitParams.yref = pFit(4);

% DATA CLASS
out.data.Rxyh = R_xyh_S;
out.data.rxyh = r_xyh_S;
out.data.rxy = r_xy_S;
% out.data.RxyhS = R_xyh_S;
% out.data.rxyhS = r_xyh_S;
% out.data.rxyS = r_xy_S;


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

% FOR PABLO
out.PJ.VEMnum = VEM_num;
out.PJ.VEMden = VEM_den;


%% NESTED FUNCTIONS
    function [RXYHM, fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins)
    % INPUT
    %   rF:         R_data(x,y,H) (10x10x10)
    %   rP:         r(x,y) (10x10)
    %   rCutOff:    min firing rate for a bin to be considered
    %   Nbins:      number of bins of the discretization
    %   pFit:       [g; thetaP; xref; yref]; initial conditions of params.
    
    % unroll pFit
    g = pFit(1);
    thetaP = pFit(2); % (deg)
    xref = pFit(3);
    yref = pFit(4);
    
    RXYHM = zeros(Nbins,Nbins,Nbins).*nan;
    fF = 0;
        angCount = 0;
        for ii=1:Nbins
            for jj=1:Nbins
                if (rP(ii, jj)>rCutOff)
                    rT = squeeze(rF(ii, jj, :));
                    iF = find(isfinite(rT));
                    if ~isempty(iF)            
                        a = 180*atan2(yref-ii, xref-jj)/pi - ...
                            (-180 -360/(2*Nbins) + iF*360/Nbins);
                        cFac = cos(pi*(a-thetaP)/180);
                        % expected value
                        cBar = nanmean(cFac);
                        % shape the cosine
                        z = 1+g*(cFac - cBar);
                        z = z.*(z>0);
                        % variance explained
                        fF = fF + nansum((z - rT(iF)).^2);                   
                        angCount = angCount + length(iF);
                        % model firing rates
                        RXYHM(ii,jj,iF) = z;
                    end
                end
            end
        end
        fF = fF/angCount;
    end

    function [fF] = aFitLMNew(pFit,rF,rP,rCutOff,Nbins)
        % call the other function (optimization)
        [~, fF] = get_Rxyh_model(pFit,rF,rP,rCutOff,Nbins);
    end

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
end

