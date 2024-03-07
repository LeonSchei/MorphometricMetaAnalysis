function [b,brmse,sk,n,msz] = regr_xzw(X,z,w);
% 
% [b,brmse,sk,n,msz] = regr_xzw(X,z,w);
% 
% general linear regression call
%
% Input
%   X, nxm dependendant variables
%   z, nx1 observations
%   OPTIONAL w, nx1 weights (0 means observation has no influence)
%
% Output
%   b, mx1 estimated parameters: z^ = X*b;
%   brmse, mx1 estimated variances (root mean square error) of parameter(s)
%      (confidence intervals assume gaussian dist., with bmse estimated variance)
%   sk, the model skill
%   n, the effective dof =(sum(w)/max(w))
%   msz, variance of data

[n,m] = size(X);
nz = length(z);
if nargin==2
    w = ones(n,1);
    nw = n;
else
    nw = length(w);
end
if(nz~=n | nw~=n | nw~=nz)
    fprintf('X and z or w are different lengths\n')
    return
end
b = nan*X(1,:)';
brmse = b;
sk=nan;

% find nans by summing
id = find(isfinite([X,z,w]*ones(m+2,1))==1);
if length(id)<2
    fprintf('n<2 -- exiting\n')
    return
end

% number of dof
n = sum(w(id))/max(w(id));

% convert to weighted space
z = (z).*w;
X = X.*(repmat(w,1,m));

% and compute model-data cov.
XX = (X(id,:)'*X(id,:))/n;
XZ = (z(id)'*X(id,:))/n;
msz = (z(id)'*z(id))/n;

% solve the weighted least squares problem
if(rcond(XX)<eps)
    return;
end
XX_inv = inv(XX);
b = XX_inv*XZ'; 

% model residuals
msr = msz - (b')*XX*(b);
sk = 1-msr/msz;

mse = XX_inv(1)*msr/(n-m);

% and perhaps we want all variance estimates
brmse = sqrt(diag(XX_inv)*msr/(n-m));