% // Mesh generated on 08-Mar-2021
% // Local directory: /Users/lukewilson/OneDrive/sp22/coe347/openfoam2/
% // Username: lukewilson 
% 
% 
% Re 20 C
Lf = 4;
Lw = 6;
R  = 1.5; % A:1 A2:1 B:1.5
H  = 6;
v = vertecies(Lf, Lw, R,H)
nyCylinder = 14; % A:10 B:14
nxCylinder = 28; % A:20 B:28 C:40
nz =1;
nHorizontalWake = 50; %A:30 B:42 C:50
nHorizontalFore = 15; %A:15 B:21 
nVertical = 42;       %A:30 B:42 C:50

blocks = block(nxCylinder,nyCylinder, nz, nHorizontalWake,nVertical,nHorizontalFore)
[arcs] = arc(blocks,R)

inlets = inlet(blocks,v)
outlets = outlet(blocks,v)
tops = top (blocks,v)
bottoms = bottom(blocks,v)
cylinders = cylinder(blocks,v)

publish(blocks,v,arcs,inlets,outlets,tops,bottoms,cylinders)

% Now to pring out the file
function  [v] = vertecies(Lf, Lw, R,H)
% setting the coordinate system to be (0,0,0) at the center of the D = 1m
% cylinder
v = zeros(64,4);

% vertecies at surface of cylinder
for i = 1:8
    theta = 2*pi*(1/8*(i-1));
    v(i,1) = 0.5*cos(theta);
    v(i,2) = 0.5*sin(theta);
end

% vertecies at distance R of cylinder surface
for i = 9:16
    theta = 2*pi*(1/8*(i-9));
    v(i,1) = (R)*cos(theta);
    v(i,2) = (R)*sin(theta);
end

% vertecies at boundary of the domain 
vertical_mesh = [0, v(10,2),H,v(10,2),0];
horizontal_mesh_RHS= [v(10,1),0,-v(10,1),-Lf];

for i = 17:19
    v(i,1) = Lw;
    v(i,2) = vertical_mesh(i-16);
end
for i = 20:23
    v(i,2) = H;
    v(i,1) = horizontal_mesh_RHS(i-19);
end
for i = 24:25
    v(i,1) = -Lf;
    v(i,2) = vertical_mesh(i-20);
end

for i = 26:32
    v(i,1) = v(i-(i-25)*2,1);
    v(i,2) = -v(i-(i-25)*2,2);
end

% Set vertecies for z axis 
v(33:64,1:2) = v(1:32,1:2);
for  i = 1:32
    v(i,3) = -0.05;
    v(i+32,3) = 0.05;
    
end

% openFoam Indexing
for i = 1:64
    v(i,4) = i-1;
end

% vx = v(:,1);
% vy = v(:,2);
% plot(vx,vy,'o');

end

function blocks = block(nxCylinder,nyCylinder, nz, nHorizontalWake,nVertical,nHorizontalFore)
% nxCylinder: number of grids in direction of surface of cylinder 
% nyCylinder: number of grids perpendicular to surface of cylinder 
% nz:  1
% nHorizontalWake: number of grids in global x direction in front of
% cylinder
% nVertical: number of grids in global y direction 
% nHorizontalFore:  number of grids in global x direction in back of cylinder


blocks = zeros(20,15);
% blocks around boundary layer top of cylinder
% starting at local bottom left point
for i = 1:7
    % vertecies
    blocks(i,1) = i-1;
    blocks(i,2) = i+7;
    blocks(i,3) = i+8;
    blocks(i,4) = i;

    % grid 
    blocks(i,9) = nyCylinder;
    blocks(i,10) = nxCylinder;
    blocks(i,11) = nz;
    
    blocks(i,12) = 2; % A:2 B:4
    blocks(i,13) = 1;
    blocks(i,14) = 1;
end

for i = 8
    blocks(i,1) = 7;
    blocks(i,2) = 15;
    blocks(i,3) = i;
    blocks(i,4) = 0;

    % grid 
    blocks(i,9) = nyCylinder;
    blocks(i,10) = nxCylinder;
    blocks(i,11) = nz;

    blocks(i,12) = 2; % A:2 B:4
    blocks(i,13) = 1;
    blocks(i,14) = 1;
    
end

% blocks outside boundary layer with arc 
for i = 9
    blocks(i,1) = 8;
    blocks(i,2) = 16;
    blocks(i,3) = 17;
    blocks(i,4) = 9;
    
    % grid 
    blocks(i,9) = nHorizontalWake;
    blocks(i,10) = nxCylinder;
    blocks(i,11) = nz;
    
    blocks(i,12) = 4; % A:4 B:4
    blocks(i,13) = 1; 
    blocks(i,14) = 1;
    
end

for i = 11:12
    blocks(i,1) = i-1;
    blocks(i,2) = i-2;
    blocks(i,3) = i+8;
    blocks(i,4) = i+9;

    % grid 
    blocks(i,9) = nxCylinder;
    blocks(i,10) = nHorizontalWake;
    blocks(i,11) = nz;

    blocks(i,12) = 1; % A:1 B:4
    blocks(i,13) = 4; % A:4
    blocks(i,14) = 1;
    
end

for i = 14:15
    blocks(i,1) = i+10;
    blocks(i,2) = i-2;
    blocks(i,3) = i-3;
    blocks(i,4) = i+9;
    
    % grid 
    blocks(i,9) = nHorizontalFore;
    blocks(i,10) = nxCylinder;
    blocks(i,11) = nz;
end
blocks(14:20,12:14) = [.25 1 1; .24 1 1; .25 .25 1; 1 .25 1; 1 .25 1; 4 .25 1;4 1 1];

for i = 17:18
    blocks(i,1) = i+10;
    blocks(i,2) = i+11;
    blocks(i,3) = i-3;
    blocks(i,4) = i-4;

    % grid 
    blocks(i,9) = nxCylinder;
    blocks(i,10) = nVertical;
    blocks(i,11) = nz;
end

for i = 20
    blocks(i,1) = 15;
    blocks(i,2) = 31;
    blocks(i,3) = 16;
    blocks(i,4) = 8;

    % grid 
    blocks(i,9) = nHorizontalWake;
    blocks(i,10) = nxCylinder;
    blocks(i,11) = nz;
    
end

% blocks in corners of boundary with out arc
blocks(10,1:4) = [9,17,18,19];
blocks(13,1:4) = [23,11,21,22];
blocks(16,1:4) = [26,27,13,25];
blocks(19,1:4) = [29,30,31,15];

blocks(13,12) = .25; % A:4 B:4
blocks(13,13) = 4; % A:4
blocks(13,14) = 1;
% grids 
for i = [10,19]
    blocks(i,9:11) = [nHorizontalWake,nVertical,nz];

    blocks(i,12) = 4; % A:4 B:4
    blocks(i,13) = 4; % A:4
    blocks(i,14) = 1;
end

for i = 16:19
    blocks(19,13) = .25
end

for i = [13,16]
    blocks(i,9:11) = [nHorizontalFore,nVertical,nz];
end

for i = 1:20
    blocks(i,5) = blocks(i,1)+32;
    blocks(i,6) = blocks(i,2)+32;
    blocks(i,7) = blocks(i,3)+32;
    blocks(i,8) = blocks(i,4)+32;
end

% OpenFoam block Index
for i = 1:20
    blocks(i,15) = i-1;
end
%blocks(9:20,12:14) = ones(12,3);
end

function [arcs] = arc(blocks,R)
arcs = zeros(32,5);
% blocks with arc edges 
block_arcs_Index = blocks(1:8,1:8);

next_arc = 1;
n = 0;
for  i = 1:8
    % Midpoint 
    theta = pi*(1/8*(i+n));
    n = n+1;
    % % vertecies in each arc
    arcs(next_arc+1,1:2) = block_arcs_Index(i,2:3);
    arcs(next_arc+1,3:5) = [R*cos(theta),R*sin(theta),-0.05];

    arcs(next_arc,1:2) = block_arcs_Index(i,[1,4]);
    arcs(next_arc,3:5) = [0.5*cos(theta),0.5*sin(theta),-0.05];

    arcs(next_arc+3,1:2) = block_arcs_Index(i,6:7);
    arcs(next_arc+3,3:5) = [R*cos(theta),R*sin(theta),0.05];

    arcs(next_arc+2,1:2) = block_arcs_Index(i,[5,8]);
    arcs(next_arc+2,3:5) = [0.5*cos(theta),0.5*sin(theta),0.05];

    next_arc = next_arc+4;
end


end

function [inlets] = inlet(blocks,v)
inlets = [];
left_end = min(v(:,1));
for i = 1:(length(blocks))
    if (v(blocks(i,1)+1,1) == left_end) & (v(blocks(i,4)+1,1) == left_end)
        inlets = [inlets; blocks(i,4) blocks(i,4)+32 blocks(i,1)+32 blocks(i,1)];
    end
end
end

function [outlets] = outlet(blocks,v)
outlets = [];
right_end = max(v(:,1));
for i = 1:(length(blocks))
    if (v(blocks(i,2)+1,1) == right_end) & (v(blocks(i,3)+1,1) == right_end)
        outlets = [outlets; blocks(i,3) blocks(i,3)+32 blocks(i,2)+32 blocks(i,2)];
    end
end
outlets = outlets([2;1;4;3],:);
end

function [tops] = top(blocks,v)

tops = [];
top_end = max(v(:,2));
for i = (length(blocks)):-1:1
    if (v(blocks(i,3)+1,2) == top_end) & (v(blocks(i,4)+1,2) == top_end)
        tops = [tops; blocks(i,4) blocks(i,4)+32 blocks(i,3)+32 blocks(i,3)];
    end
end
end

function [bottoms] = bottom(blocks,v)
bottoms = [];
bottom_end = min(v(:,2));
for i = 1:(length(blocks))
    if (v(blocks(i,1)+1,2) == bottom_end) & (v(blocks(i,2)+1,2) == bottom_end)
        bottoms = [bottoms; blocks(i,1) blocks(i,1)+32 blocks(i,2)+32 blocks(i,2)];
    end
end
end

function [cylinders] = cylinder(blocks,v)
cylinders = [];
rad = 0.5;
for i = 1:(length(blocks))
    if (norm(v(blocks(i,1)+1,1:2),2) - rad < 0.0001) & (norm(v(blocks(i,4)+1,1:2),2) - rad < 0.0001)
        cylinders = [cylinders; blocks(i,1) blocks(i,1)+32 blocks(i,4)+32 blocks(i,4)];
    end
end
end

function publish(blocks,v,arcs,inlets,outlets,tops,bottoms,cylinders)

fid = fopen('blockMeshDictB2','wt');
% write header
fprintf(fid,'// Mesh generated on 3/24/22\n');
fprintf(fid,'// Local directory: /Users/fbisetti/teaching/ASE347/SP-2021/OF-assignments/of2/simulations/mesh\n');
fprintf(fid,'// Username: fbisetti\n\n\n');

fprintf(fid,'// Lf (fore)  =  4.00000e+00\n');
fprintf(fid,'// Lw (wake)  =  6.00000e+00\n');
fprintf(fid,'// R  (outer) =  1.00000e+00\n');
fprintf(fid,'// H  (top/bottom) = 4.00000\n\n');

fprintf(fid,'FoamFile\n');
fprintf(fid,'{\n\tversion  2.0;\n\tformat   ascii;\n\tclass    dictionary;\n\tobject   blockMeshDict;\n}\n\n');

fprintf(fid,'convertToMeters 1.0;\n\n');

% write verticies
fprintf(fid,'vertices\n(\n');
fprintf(fid,'\t( %e %e %e ) // %1.0f\n',v(:,1:4)');
fprintf(fid,');\n\n');

% write blocks
fprintf(fid,'blocks\n(\n\n');
fprintf(fid,'\t// block %1.0f \n\thex (%2.0f %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f %2.0f) ( %2.0f %2.0f %2.0f) simpleGrading ( %e %e %3.1f)\n',[blocks(:,end) blocks(:,1:(end-1))]');
fprintf(fid,'\n);\n\n');

% write arcs
fprintf(fid,'edges\n(\n\n');
fprintf(fid,'arc %2.0f %2.0f ( %e %e %e )\n',arcs');
fprintf(fid,'\n\n);\n\n');

% write boundaries
fprintf(fid,'boundary\n(\n\n');

fprintf(fid,'\tinlet\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n');
fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',inlets');
fprintf(fid,'\t\t);\n\t}\n\n');

fprintf(fid,'\toutlet\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n');
fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',outlets');
fprintf(fid,'\t\t);\n\t}\n\n');

fprintf(fid,'\tcylinder\n\t{\n\t\ttype wall;\n\t\tfaces\n\t\t(\n');
fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',cylinders');
fprintf(fid,'\t\t);\n\t}\n\n');

fprintf(fid,'\ttop\n\t{\n\t\ttype symmetryPlane;\n\t\tfaces\n\t\t(\n');
fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',tops');
fprintf(fid,'\t\t);\n\t}\n\n');

fprintf(fid,'\tbottom\n\t{\n\t\ttype symmetryPlane;\n\t\tfaces\n\t\t(\n');
fprintf(fid,'\t\t\t(%2.0f %2.0f %2.0f %2.0f)\n',bottoms');
fprintf(fid,'\t\t);\n\t}\n\n');

fprintf(fid,');');

fclose(fid);

end
