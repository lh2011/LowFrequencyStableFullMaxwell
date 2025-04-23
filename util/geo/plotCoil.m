[nrbarr, condVec, dirHigh, dirLow, epsVec] = cnstrct_paper_example_config(10, 2, 1);

copperColor = '#CC4C03';
dir1Color = '#C9308E';
dir2Color = '#009CDA';
boundingColor = '#008877';

figure(1)
clf();
hold on;
for i=setdiff(find(condVec),[dirLow,dirHigh])
    mynrbplot(nrbarr(i),[0,0,0],copperColor);
end
bnd = nrbextract(nrbarr(dirHigh));
mynrbplot(bnd(2),[1,1],dir1Color);
mynrbplot(bnd(3),[1,1],copperColor);
mynrbplot(bnd(4),[1,1],copperColor);
mynrbplot(bnd(5),[1,1],copperColor);
mynrbplot(bnd(6),[1,1],copperColor);
bnd = nrbextract(nrbarr(dirLow));
mynrbplot(bnd(1),[1,1],copperColor);
mynrbplot(bnd(2),[1,1],copperColor);
mynrbplot(bnd(3),[1,1],copperColor);
mynrbplot(bnd(4),[1,1],copperColor);
mynrbplot(bnd(5),[1,1],dir2Color);

boundingBox = inf*ones(1,6);
for i=1:numel(nrbarr)
    p = nrbeval(nrbarr(i),{[0,1],[0,1],[0,1]});
    xmin = min(p(1,:,:),[],'all');
    xmax = max(p(1,:,:),[],'all');
    ymin = min(p(2,:,:),[],'all');
    ymax = max(p(2,:,:),[],'all');
    zmin = min(p(3,:,:),[],'all');
    zmax = max(p(3,:,:),[],'all');
    
    if isinf(boundingBox(1))
        boundingBox(1) = xmin;
    else
        boundingBox(1) = min([boundingBox(1),xmin]);
    end

    if isinf(boundingBox(2))
        boundingBox(2) = xmax;
    else
        boundingBox(2) = max([boundingBox(2),xmax]);
    end

    if isinf(boundingBox(3))
        boundingBox(3) = ymin;
    else
        boundingBox(3) = min([boundingBox(3),ymin]);
    end

    if isinf(boundingBox(4))
        boundingBox(4) = ymax;
    else
        boundingBox(4) = max([boundingBox(4),ymax]);
    end

    if isinf(boundingBox(5))
        boundingBox(5) = zmin;
    else
        boundingBox(5) = min([boundingBox(5),zmin]);
    end

    if isinf(boundingBox(6))
        boundingBox(6) = zmax;
    else
        boundingBox(6) = max([boundingBox(6),zmax]);
    end
end
plot3([boundingBox(1),boundingBox(2)],[boundingBox(3),boundingBox(3)],[boundingBox(5),boundingBox(5)],Color=boundingColor);
plot3([boundingBox(1),boundingBox(2)],[boundingBox(4),boundingBox(4)],[boundingBox(5),boundingBox(5)],Color=boundingColor);
plot3([boundingBox(1),boundingBox(2)],[boundingBox(3),boundingBox(3)],[boundingBox(6),boundingBox(6)],Color=boundingColor);
plot3([boundingBox(1),boundingBox(2)],[boundingBox(4),boundingBox(4)],[boundingBox(6),boundingBox(6)],Color=boundingColor);

plot3([boundingBox(1),boundingBox(1)],[boundingBox(3),boundingBox(4)],[boundingBox(5),boundingBox(5)],Color=boundingColor);
plot3([boundingBox(2),boundingBox(2)],[boundingBox(3),boundingBox(4)],[boundingBox(5),boundingBox(5)],Color=boundingColor);
plot3([boundingBox(1),boundingBox(1)],[boundingBox(3),boundingBox(4)],[boundingBox(6),boundingBox(6)],Color=boundingColor);
plot3([boundingBox(2),boundingBox(2)],[boundingBox(3),boundingBox(4)],[boundingBox(6),boundingBox(6)],Color=boundingColor);

plot3([boundingBox(1),boundingBox(1)],[boundingBox(3),boundingBox(3)],[boundingBox(5),boundingBox(6)],Color=boundingColor);
plot3([boundingBox(2),boundingBox(2)],[boundingBox(3),boundingBox(3)],[boundingBox(5),boundingBox(6)],Color=boundingColor);
plot3([boundingBox(1),boundingBox(1)],[boundingBox(4),boundingBox(4)],[boundingBox(5),boundingBox(6)],Color=boundingColor);
plot3([boundingBox(2),boundingBox(2)],[boundingBox(4),boundingBox(4)],[boundingBox(5),boundingBox(6)],Color=boundingColor);
axis('vis3d');
axis off;

%Play  around with rotation in figure until you like the perspective
% [a,b] = view() in console to get parameters 
view([63.1810,-45.2400]);

exportgraphics(gcf(),'test.pdf','ContentType','vector');