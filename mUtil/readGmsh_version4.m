function [p, TEMP] = readGmsh_version4(filename)

fid=fopen([filename '.msh'],'r');

while 1
    temp=fgetl(fid);
    if strcmp(temp,'$Nodes') || strcmp(temp,'$NOD')
        break
    end
end

temp = fscanline(fid, '%d', 2);
Nentities = temp{1,1};
Nnd = temp{2,1};
p=zeros(Nnd,3);
disp('Reading points...');
for ii = 1:Nentities
    temp = fscanline(fid, '%d', 4);
    n = temp{4,1};
    id=fscanf(fid,'%d',n);
    tmp = fscanf(fid,'%f',3*n)';
    p(id,:) = reshape(tmp,3,n)';
end

disp('Reading Elements...');
while 1
    temp=fgetl(fid);
    if strcmp(temp,'$Elements') || strcmp(temp,'$ELM')
        break
    end
end

temp = fscanline(fid, '%d', 2);
Nentities = temp{1,1};
Nel = temp{2,1};
TEMP=nan(Nel,30);


Nnodes = Nodes_per_element();
not_tested = [3 8 9 10 16];


for ii = 1:Nentities
    temp = fscanline(fid, '%d', 4);
    elcode = temp{3,1};
    if ~isempty(find(not_tested == elcode, 1))
        warning(['Element type: ' num2str(elcode) ' has not been tested'])
    end
        
    
    n = temp{4,1};
    m = Nnodes(Nnodes(:,1) == elcode,2);
    tmp = fscanf(fid, '%f',m*n)';
    tmp = reshape(tmp, m, n)';
    TEMP(tmp(:,1), 1:m) = [double(elcode)*ones(n,1) double(tmp(:,2:end))];
end



fclose(fid);
end

function out = fscanline(fid, format, N)
    temp = [];
    while isempty(temp) || strcmp(temp, ' ')
        temp = fgetl(fid);
    end

    C = textscan(temp, format, N);
    for ii = 1:N
        out{ii,1} = C{1,1}(ii);
    end
end


function Nnodes = Nodes_per_element()
% Returns a map between element code and the number of numbers that have to
% be read
Nnodes = [15 2; ... % 1-node point.
          1 3; ... % 2-node line.
          2 4; ... % 3-node triangle.
          3 5; ... % 4-node quadrangle.
          8 4; ... % 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
          9 7; ... % 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
          10 10; ... % 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
          16 9; ... % 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
          ];
end

