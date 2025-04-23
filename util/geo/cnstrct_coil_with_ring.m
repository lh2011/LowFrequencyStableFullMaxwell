function [nrbarr, condVec, dirHigh, dirLow, epsVec,regionVec] = cnstrct_coil_with_ring(dimCond)

%% define values
a = dimCond;
b = dimCond/100;

x0 = 0; % x offset if wanted
y0 = 0; % y offset if wanted

x = x0 + [0     a   2*a     3*a     4*a     5*a     6*a     7*a     8*a];
y = y0 + [0     a   2*a     3*a     4*a     5*a     5*a+b   6*a+b   7*a+b   8*a+b   9*a+b];
%         1     2   3       4       5       6       7       8       9       10      11 

%% hardcoded pattern

condVec = [0 0 0 0 0 0 0 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 0];
epsVec = ones(size(condVec));
dirHigh = 9;
dirLow = 22 + numel(condVec);

nrbsurfaces(1) = nrb4surf([x(1), y(1)], [x(2) y(2)], [x(1), y(6)], [x(2) y(6)]);
nrbsurfaces(2) = nrb4surf([x(1), y(6)], [x(2) y(6)], [x(1), y(7)], [x(2) y(7)]);
nrbsurfaces(3) = nrb4surf([x(1), y(7)], [x(2) y(7)], [x(1), y(11)], [x(2) y(10)]);
nrbsurfaces(4) = nrb4surf([x(1), y(11)], [x(2) y(10)], [x(8), y(11)], [x(7) y(10)]);
nrbsurfaces(5) = nrb4surf([x(8), y(11)], [x(7) y(10)], [x(8), y(7)], [x(7) y(7)]);
nrbsurfaces(6) = nrb4surf([x(8), y(7)], [x(7) y(7)], [x(8), y(6)], [x(7) y(6)]);

%nrbsurfaces(7) = nrb4surf([x(8), y(6)], [x(7) y(6)], [x(8), y(4)], [x(7) y(4)]);

nrbsurfaces(7) = nrb4surf([x(8), y(6)], [x(7) y(6)], [x(8), y(5)], [x(7) y(5)]);

nrbsurfaces(8) = nrb4surf([x(8), y(4)], [x(7) y(4)], [x(8), y(3)], [x(7) y(3)]);
nrbsurfaces(9) = nrb4surf([x(8), y(3)], [x(7) y(3)], [x(8), y(2)], [x(7) y(2)]);
nrbsurfaces(10) = nrb4surf([x(8), y(2)], [x(7) y(2)], [x(8), y(1)], [x(7) y(1)]);

%nrbsurfaces(11) = nrb4surf([x(1), y(1)], [x(7) y(1)], [x(2), y(2)], [x(7) y(2)]);
nrbsurfaces(11) = nrb4surf([x(1), y(1)], [x(6) y(1)], [x(2), y(2)], [x(6) y(2)]);

nrbsurfaces(12) = nrb4surf([x(2) y(2)], [x(3) y(3)], [x(2) y(6)], [x(3) y(6)]);
nrbsurfaces(13) = nrb4surf([x(2) y(6)], [x(3) y(6)], [x(2) y(7)], [x(3) y(7)]);
nrbsurfaces(14) = nrb4surf([x(2) y(7)], [x(3) y(7)], [x(2) y(10)], [x(3) y(9)]);

nrbsurfaces(15) = nrb4surf([x(2) y(10)], [x(3) y(9)], [x(7) y(10)], [x(6) y(9)]);
nrbsurfaces(16) = nrb4surf([x(7) y(10)], [x(6) y(9)], [x(7) y(7)], [x(6) y(7)]);
nrbsurfaces(17) = nrb4surf([x(7) y(7)], [x(6) y(7)], [x(7) y(6)], [x(6) y(6)]);

% nrbsurfaces(18) = nrb4surf([x(7) y(6)], [x(6) y(6)], [x(7) y(4)], [x(6) y(5)]);
% nrbsurfaces(19) = nrb4surf([x(7) y(4)], [x(6) y(5)], [x(4) y(4)], [x(5) y(5)]);

nrbsurfaces(18) = nrb4surf([x(7) y(6)], [x(6) y(6)], [x(7) y(5)], [x(6) y(5)]);
nrbsurfaces(19) = nrb4surf([x(6) y(4)], [x(6) y(5)], [x(4) y(4)], [x(5) y(5)]);

nrbsurfaces(20) = nrb4surf([x(4) y(4)], [x(5) y(5)], [x(4) y(6)], [x(5) y(6)] );
nrbsurfaces(21) = nrb4surf([x(4) y(6)], [x(5) y(6)], [x(4) y(7)], [x(5) y(7)] );
nrbsurfaces(22) = nrb4surf([x(4) y(7)], [x(5) y(7)], [x(4) y(8)], [x(5) y(8)] );
 
%nrbsurfaces(23) = nrb4surf([x(7), y(3)], [x(7), y(4)], [x(3), y(3)], [x(4), y(4)]);
nrbsurfaces(23) = nrb4surf([x(6), y(3)], [x(6), y(4)], [x(3), y(3)], [x(4), y(4)]);

nrbsurfaces(24) = nrb4surf([x(3), y(3)], [x(4), y(4)], [x(3), y(6)], [x(4), y(6)]);
nrbsurfaces(25) = nrb4surf([x(3), y(6)], [x(4), y(6)], [x(3), y(7)], [x(4), y(7)]);
nrbsurfaces(26) = nrb4surf([x(3), y(7)], [x(4), y(7)], [x(3), y(9)], [x(4), y(8)]);
nrbsurfaces(27) = nrb4surf([x(3), y(9)], [x(4), y(8)], [x(6), y(9)], [x(5), y(8)]);

nrbsurfaces(28) = nrb4surf([x(6), y(9)], [x(5), y(8)], [x(6), y(7)], [x(5), y(7)]);
nrbsurfaces(29) = nrb4surf([x(6), y(7)], [x(5), y(7)], [x(6), y(6)], [x(5), y(6)]);
nrbsurfaces(30) = nrb4surf([x(6), y(6)], [x(5), y(6)], [x(6), y(5)], [x(5), y(5)]);

%nrbsurfaces(31) = nrb4surf([x(2), y(2)], [x(7), y(2)], [x(3), y(3)], [x(7), y(3)]);
nrbsurfaces(31) = nrb4surf([x(2), y(2)], [x(6), y(2)], [x(3), y(3)], [x(6), y(3)]);

nrbsurfaces(32) = nrb4surf([x(6) y(4)], [x(7) y(4)],[x(6) y(5)], [x(7) y(5)]);

nrbsurfaces(33) = nrb4surf([x(6), y(3)], [x(7), y(3)], [x(6), y(4)], [x(7), y(4)]);

nrbsurfaces(34) = nrb4surf([x(8), y(5)], [x(7) y(5)], [x(8), y(4)], [x(7) y(4)]);

nrbsurfaces(35) = nrb4surf([x(6), y(2)], [x(7), y(2)], [x(6), y(3)], [x(7), y(3)]);
nrbsurfaces(36) = nrb4surf([x(6), y(1)], [x(7), y(1)], [x(6), y(2)], [x(7), y(2)]);

%% extrude surfaces
nrbarr = [];
for i = 1:numel(nrbsurfaces)
    nrbarr = [nrbarr, nrbextrude(nrbsurfaces(i), [0, 0, a])];
end

%% add air layer on bottom
condVec_tmp = zeros(size(nrbsurfaces));
condVec_tmp(22) = 1;
condVec = [condVec, condVec_tmp];
epsVec = [epsVec, ones(size(condVec_tmp))];


nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 -a]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add air layer on top
condVec = [condVec, zeros(size(nrbsurfaces))];
epsVec = [epsVec, ones(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
     nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 a]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add shorted coil

fac = 10;

gap = [13];
coil = [12 13 14 15 16 17 18 31 32 33 35];

condVecShortedCoil = zeros(size(nrbsurfaces));
condVecShortedCoil(coil) = 1;
condVecShortedCoil(gap) = 0;
regionVecShortedCoil = 2 * condVecShortedCoil;
regionVecShortedCoil(gap) = 3;
regionVec = [condVec, regionVecShortedCoil];
condVec = [condVec, condVecShortedCoil];

epsVecShortedCoil = ones(size(nrbsurfaces));
epsVecShortedCoil(gap) = 6e7/(2*pi*150*8.854187e-12);
epsVec = [epsVec, epsVecShortedCoil];

nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbtform(nrbarr(i), vecscale([1 1 1/fac])), vectrans([0 0 2*a]))];
end
nrbarr = [nrbarr, nrbarr_tmp];

%% add air layer on top
condVec = [condVec, zeros(size(nrbsurfaces))];
regionVec = [regionVec, zeros(size(nrbsurfaces))];
epsVec = [epsVec, ones(size(nrbsurfaces))];
nrbarr_tmp = [];
for i=1:numel(nrbsurfaces)
    nrbarr_tmp = [nrbarr_tmp, nrbtform(nrbarr(i), vectrans([0 0 (2+1/fac)*a]))];
end
nrbarr = [nrbarr, nrbarr_tmp];



end
