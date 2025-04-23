clear all;

%% Analytical considerations
syms x y z om nu sig eps;

% A = [sin(x)*sin(y)*sin(z);0;0];
% B = curl(A,[x,y,z]);
% ccA = curl(nu*B,[x,y,z]);
% 
% D = (1i*om*sig-om^2*eps)*A;
% 
% rhsFull = (D + ccA);
% 
% phi = sin(x)*sin(y)*sin(z);
% kgPhi = (sig+1i*om*eps)*gradient(phi,[x,y,z]);
% rhsEQS = -divergence(kgPhi,[x,y,z]);
% 
% Js = kgPhi + rhsFull;

A = [sin(x)*cos(y)*cos(z);...
     -2*cos(x)*sin(y)*cos(z);...
     cos(x)*cos(y)*sin(z)];

B = curl(A,[x,y,z]);
ccA = curl(nu*B,[x,y,z]);

checkDivA = divergence(A,[x,y,z])

D = (1i*om*sig-om^2*eps)*A;
checkDivD = divergence(D,[x,y,z])

rhsFull = (D + ccA);

phi = cos(x) * cos(y) * cos(z); %sin(x)*sin(y)*sin(z);
kgPhi = (sig+1i*om*eps)*gradient(phi,[x,y,z]);
rhsEQS = -divergence(kgPhi,[x,y,z]);

Js = kgPhi + rhsFull;
checkDivJs = divergence(Js,[x,y,z])


%% Geometry is the pi-cube
srf = nrbsquare([0 0],pi,pi);
vol = nrbextrude(srf,[0 0 pi]);

%% Parameter and materials
om = 50; %Hz

nu_const = 1;
sigma_const = 1;
epsi_const = 1;

nuFun  = @(x, y, z) nu_const*ones(size(x));
sigmaFun = @(x, y, z) sigma_const*ones(size(x));
epsiFun = @(x, y, z) epsi_const*ones(size(x));
kappaFun = @(x, y, z) sigmaFun(x,y,z) + 1i*om*epsiFun(x,y,z);

nu_mat = @(x, y, z) cat(1,reshape(nuFun(x,y,z),[1,size(x)]),...
                          reshape(nuFun(x,y,z),[1,size(x)]),...
                          reshape(nuFun(x,y,z),[1,size(x)]));
kappa_mat = @(x, y, z) cat(1,reshape(kappaFun(x,y,z),[1,size(x)]),...
                             reshape(kappaFun(x,y,z),[1,size(x)]),...
                             reshape(kappaFun(x,y,z),[1,size(x)]));

% Analytical solutions
Afun = @(x, y, z) cat(1, ...
    reshape(sin(x).*cos(y).*cos(z),[1,size(x)]),...
    reshape(-2*cos(x).*sin(y).*cos(z),[1,size(x)]),...
    reshape(cos(x).*cos(y).*sin(z),[1,size(x)]));

Bfun = @(x, y, z) cat(1, ...
    reshape(-3*cos(x).*sin(y).*sin(z),[1,size(x)]),...
    zeros([1,size(x)]),...
    reshape(3*sin(x).*sin(y).*cos(z),[1,size(x)]));

Ufun = @(x, y, z) sin(x).*sin(y).*sin(z);

% negative E-field - just grad(u)
nEfun = @(x, y, z) cat(1,reshape(cos(x).*sin(y).*sin(z),[1,size(x)]),...
                         reshape(sin(x).*cos(y).*sin(z),[1,size(x)]),...
                         reshape(sin(x).*sin(y).*cos(z),[1,size(x)]));

%% Construct right hand sides
ccAfun = @(x, y, z) cat(1,reshape(3*sin(x).*cos(y).*cos(z),[1,size(x)]),...
                          reshape(-6*cos(x).*sin(y).*cos(z),[1,size(x)]),...
                          reshape(3*cos(x).*cos(y).*sin(z),[1,size(x)]));

rhsVecFun = @(x, y, z) nu_mat(x,y,z).*ccAfun(x,y,z) + 1i*om*kappa_mat(x,y,z).*Afun(x,y,z);
rhsScaFun = @(x, y, z) 3*kappa_fun(x,y,z).*Ufun(x,y,z);

% rhsVec = Js - kappa*grad(u) <=> Js = rhsVec + kappa*grad(u)
JsFun = @(x, y, z) rhsVecFun(x, y, z) + kappa_mat(x,y,z) .* nEfun(x,y,z);

%% Dirichlet Boundary conditions
gVecFun = Afun(x,y,z); % Use sp_drchlt_l2_proj() -> Enforcing tangential part is not necessary (is done implicitly in this version)
gScaFun = Ufun(x,y,z);