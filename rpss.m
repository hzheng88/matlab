function z = rpss(fcst, obs, edges, enssize)
%  RPSS Calculates ranked probability skill score
% 
% The RPSS meausre the closeness of forecast distribution (fcst) and
% corresponding observation (obs).
% For more info see
%  Weigel, A. P., Liniger, M. A., & Appenzeller, C. (2007a). 
%   The Discrete Brier and Ranked Probability Skill Scores. 
%   Monthly Weather Review, 135(1), 118–124.
%   https://doi.org/10.1175/MWR3280.1
%  Weigel, A. P., Liniger, M. A., & Appenzeller, C. (2007b).
%   Generalization of the discrete brier and ranked probability skill scores for weighted multimodel ensemble forecasts.
%   Monthly Weather Review, 135(7), 2778–2785.
%   https://doi.org/10.1175/MWR3428.1
%  Weigel, A. P., & Mason, S. J. (2011).
%   The generalized discrimination score for ensemble forecasts.
%   Monthly Weather Review, 139(9), 3069–3074.
%   https://doi.org/10.1175/MWR-D-10-05069.1
%
% Usage:
%  rpss = rpss(fcst, obs, edges, enssize);
%
% Inputs
%  obs: vector of observations
%  fcst: Matrix of ensemble forecast of size Nens x Nobs, where Nens is the
%        number of ensemble members, and Nobs is the number of
%        observations.
%  edges: specify the edges of bins. Each bin includes the left edge, but
%         does not include the right edge, except for the last bin which includes
%         both edges.
%  enssize: the effective ensemble size.
%           For equal-weighted ensembles, enssize equals Nens.
%           Refer to Weigel et al. (2007b) for weighted ensembles.

if ~isvector(obs)
    error('RPSS:BadObs', 'OBS is not a vector');
end
if size(fcst,2) ~= length(obs)
    error('RPSS:BadFcst', 'size(FCST,2) does not equal length(obs)');
end
if ~isvector(edges)
    error('RPSS:BadEdges', 'EDGES is not a vector');
end

nobs = length(obs);
nens = size(fcst,1);
nbin = length(edges) - 1;

% which bin do the observations and forecasts fall in?
oi = zeros(size(obs));
fi = zeros(size(fcst));
for ii = 1:nbin-1
    oi(obs>=edges(ii) & obs<edges(ii+1)) = ii;
    fi(fcst>=edges(ii) & fcst<edges(ii+1)) = ii;
end
ii = nbin; % the last bin
oi(obs>=edges(ii) & obs<=edges(ii+1)) = ii;
fi(fcst>=edges(ii) & fcst<=edges(ii+1)) = ii;

% climatology of the observations
c = zeros(1,nbin);
for ii = 1:nbin
    c(ii) = sum(oi==ii) ./ nobs;
end
C = cumsum(c);

% forecast-observation pairs
rps = zeros(size(obs));
rpsc = zeros(size(obs));
p = zeros(1,nbin); % prob of fcst
P = zeros(1,nbin); % cum prob of fcst
o = zeros(1,nbin); % prob of obs
O = zeros(1,nbin); % cum prob of obs
for jj = 1:nobs
    % forecasts
    for ii = 1:nbin
        p(ii) = sum(fi(:,jj)==ii) ./ nens;
    end
    P(:) = cumsum(p);
    % observation
    o(:) = 0;
    o(oi(:,jj)) = 1;
    O(:) = cumsum(o);
    % RPS and RPSc
    rps(jj) = sum((P-O).^2);
    rpsc(jj) = sum((C-O).^2);
end

% debiased rpss (Weigel et al. 2007b)
d0 = (nbin .^ 2 - 1) / 6 / nbin;
z = 1 - mean(rps) / (mean(rpsc) + d0 / enssize);
end