% routine to construct geometry for academic test example taken 
% from ostrowski paper
function nrbarr = ostrowski_cube(l,d)

a = (l - d) / 2;
  
plane1 = nrb4surf([0 0], [a 0], [0 a], [a a]);
plane2 = nrb4surf([a 0], [a+d 0], [a a], [a+d a]);
plane3 = nrb4surf([a+d 0], [l 0], [a+d a], [l a]);
plane4 = nrb4surf([0 a], [a a], [0 a+d], [a a+d]);

plane5 = nrb4surf([a a], [a+d a], [a a+d], [a+d a+d]);

plane6 = nrb4surf([a+d a], [l a], [a+d a+d], [l a+d]);
plane7 = nrb4surf([0 a+d], [a a+d], [0 l], [a l]);
plane8 = nrb4surf([a a+d], [a+d a+d], [a l], [a+d l]);
plane9 = nrb4surf([a+d a+d], [l a+d], [a+d l], [l l]);

planes_1 = [plane1, plane2, plane3, plane4, plane5, plane6, plane7, plane8, plane9];

cubes_1 = [];
for i=1:numel(planes_1)
    cubes_1 = [cubes_1, nrbextrude(planes_1(i), [0, 0, a])];
end

cubes_2 = [];
for i=1:numel(planes_1)
    cubes_2 = [cubes_2, nrbtform(nrbtform(cubes_1(i), vecscale([1 1 d/a])), vectrans([0 0 a]))];
end
cubes_3 = [];
for i=1:numel(planes_1)
    cubes_3 = [cubes_3, nrbtform(cubes_1(i), vectrans([0 0 a+d]))];
end

nrbarr = [cubes_1, cubes_2, cubes_3];

end