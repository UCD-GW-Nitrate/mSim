function [p MSH]=read_2D_Gmsh(varargin)
% [p MSH]=read_2D_Gmsh(varargin)
%
% Reads the *.msh file generated from Gmsh
%
% INPUT
% filename : the name file without the extension *.msh
% display_count : when the meshes are very large this function can be time
%                 consuming. Because I'm very impatient and always want to
%                 know that my script is running I allow to display the
%                 loop counters by setting this to 1. 0 or empty does not
%                 display anything
% Use_Cpp       : Setting this to 1 will run some parts of the code in c++. 
%                 This is experimental and should be used only for very
%                 large meshes.
% 
%
% OUTPUT
% p   :   coordinates of points
% MSH : A structure with 3 rows. The first row corresponds to the zero
%       dimension elements (e.g. points) The second row to 1D elements
%       (lines) and the third row to 2D elements (triangles or
%       quadrilaterals). 
%       Each row has the field elem. The elem consist of one or more rows
%       where each row correspond to a different type of elements. If only
%       one type of element is used then elem has only one row. Each row of
%       elem has the following fields:
%       type is the type of element. for example 'triangle' 'quad'
%       id is the matrix which contains the connentivity of the points
%       the last field is dom and it is not used. 
%       
% To plot a triangular mesh in matlab or octave do the following
% First identify which row of the MSH(3,1).elem has type 'triangle'. If
% this is 1 then do :
% triplot(MSH(3,1).elem(1).id, p(:,1), p(:,2)).
% 
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 25-Mar_2013
% Department of Land Air and Water
% University of California Davis
%

if nargin == 3
    filename = varargin{1};
    display_count = varargin{2};
    Use_Cpp = varargin{3};
elseif nargin == 2
    filename = varargin{1};
    display_count = varargin{2};
    Use_Cpp = 0;
elseif nargin == 1
    filename = varargin{1};
    display_count = 0;
    Use_Cpp = 0;
else
    error('Wrong number of arguments')
end
if isempty(display_count)
    display_count = 0;
end
if isempty(Use_Cpp)
    Use_Cpp = 0;
end

if Use_Cpp
    [p, TEMP] = read_2D_mesh_cpp([filename '.msh']);
    TEMP(:,1)=[];
else
    
    fid=fopen([filename '.msh'],'r');
    
    while 1
        temp=fgetl(fid);
        if strcmp(temp,'$MeshFormat')
            temp = fscanf(fid, '%d', 3);
            GMSHVersion = floor(temp(1));
            fclose(fid);
            break
        end
    end
    
    if GMSHVersion == 4
        [p, TEMP] = readGmsh_version4(filename);
    else
        error('Older versions of GMSH are no longer supported')
    end
end
Np = size(p,1);

b=unique(TEMP(:,1));
cnt0D=1;cnt1D=1;cnt2D=1;cnt3D=1;
for i=1:size(b)
    if b(i)==15
        MSH(1,1).elem(cnt0D,1).type='BndNode';
        id= TEMP(:,1)==b(i);
        MSH(1,1).elem(cnt0D,1).id=TEMP(id,2);
        %MSH(1,1).elem(cnt0D,1).dom=TEMP_dom(id,1);
        cnt0D=cnt0D+1;
    elseif b(i)==1
        MSH(2,1).elem(cnt1D,1).type='line';
        id= TEMP(:,1)==b(i);
        MSH(2,1).elem(cnt1D,1).id=TEMP(id,2:3);
        %MSH(2,1).elem(cnt1D,1).dom=TEMP_dom(id,1);
        cnt1D=cnt1D+1;
    elseif b(i)==8
        MSH(2,1).elem(cnt1D,1).type='line';
        id= TEMP(:,1)==b(i);
        MSH(2,1).elem(cnt1D,1).id=TEMP(id,2:4);
        %MSH(2,1).elem(cnt1D,1).dom=TEMP_dom(id,1);
        cnt1D=cnt1D+1;
    elseif b(i)==2
        MSH(3,1).elem(cnt2D,1).type='triangle';
        id= TEMP(:,1)==b(i);
        MSH(3,1).elem(cnt2D,1).id=TEMP(id,2:4);
        %MSH(3,1).elem(cnt2D,1).dom=TEMP_dom(id,1);
        tf=iscw(p,MSH(3,1).elem(cnt2D,1).id);
        temp=MSH(3,1).elem(cnt2D,1).id(tf,:);
        temp=[temp(:,1) temp(:,3) temp(:,2)];
        MSH(3,1).elem(cnt2D,1).id(tf,:)=temp;
        cnt2D=cnt2D+1;
    elseif b(i)==9
        MSH(3,1).elem(cnt2D,1).type='triangle';
        id= TEMP(:,1)==b(i);
        MSH(3,1).elem(cnt2D,1).id=TEMP(id,2:7);
        %MSH(3,1).elem(cnt2D,1).dom=TEMP_dom(id,1);
        tf=iscw(p,MSH(3,1).elem(cnt2D,1).id);
        temp=MSH(3,1).elem(cnt2D,1).id(tf,:);
        temp=[temp(:,1) temp(:,3) temp(:,2) temp(:,6) temp(:,5) temp(:,4)];
        MSH(3,1).elem(cnt2D,1).id(tf,:)=temp;
        cnt2D=cnt2D+1;
    elseif b(i)==3
        MSH(3,1).elem(cnt2D,1).type='quad';
        id= TEMP(:,1)==b(i);
        MSH(3,1).elem(cnt2D,1).id=TEMP(id,2:5);
        %MSH(3,1).elem(cnt2D,1).dom=TEMP_dom(id,1);
        tf=iscw(p,MSH(3,1).elem(cnt2D,1).id);
        temp=MSH(3,1).elem(cnt2D,1).id(tf,:);
        temp=[temp(:,1) temp(:,4) temp(:,3) temp(:,2)];
        MSH(3,1).elem(cnt2D,1).id(tf,:)=temp;
        cnt2D=cnt2D+1;
    elseif b(i)==10
        MSH(3,1).elem(cnt2D,1).type='quad';
        id= TEMP(:,1)==b(i);
        MSH(3,1).elem(cnt2D,1).id=TEMP(id,2:10);
        %MSH(3,1).elem(cnt2D,1).dom=TEMP_dom(id,1);
        tf=iscw(p,MSH(3,1).elem(cnt2D,1).id);
        temp=MSH(3,1).elem(cnt2D,1).id(tf,:);
        temp=[temp(:,1) temp(:,4) temp(:,3) temp(:,2) temp(:,8) temp(:,7) temp(:,6) temp(:,5) temp(:,9)];
        MSH(3,1).elem(cnt2D,1).id(tf,:)=temp;
        cnt2D=cnt2D+1;
    elseif b(i)==16
        MSH(3,1).elem(cnt2D,1).type='quad_8';
        id= TEMP(:,1)==b(i);
        MSH(3,1).elem(cnt2D,1).id=TEMP(id,2:9);
        %MSH(3,1).elem(cnt2D,1).dom=TEMP_dom(id,1);
        tf=iscw(p,MSH(3,1).elem(cnt2D,1).id);
        temp=MSH(3,1).elem(cnt2D,1).id(tf,:);
        temp=[temp(:,1) temp(:,4) temp(:,3) temp(:,2) temp(:,8) temp(:,7) temp(:,6) temp(:,5)];
        MSH(3,1).elem(cnt2D,1).id(tf,:)=temp;
        cnt2D=cnt2D+1;
    end
end
%find if there are nodes without connections in the simulation dimension 
%delete them
% dlt=[];
% dim=size(MSH,1);
% Np = size(p,1);
% for i=1:Np
%     if display_count
%         [i Np] 
%     end
%     for j=1:size(MSH(dim,1).elem,1)
%         if isempty(find(MSH(dim,1).elem(j,1).id==i, 1))
%             dlt=[dlt;i];
%         end
%     end
% end

% This is a faster way to do the above loop but requires more memory
TestMat=sparse(Np,Np);
dim=size(MSH,1);
for k=1:size(MSH(dim,1).elem,1)
    for i = 1:size(MSH(dim,1).elem(k,1).id,2)
        for j = 1:size(MSH(dim,1).elem(k,1).id,2)
            TestMat = TestMat + sparse(MSH(dim,1).elem(k,1).id(:,i), MSH(dim,1).elem(k,1).id(:,j), ones(size(MSH(dim,1).elem(k,1).id,1),1),Np, Np);
        end
    end
end
TestMat = sum(TestMat,2);
dlt = find(TestMat == 0);


%delete the points and update the MSH structure
p(dlt,:)=[];
for i=length(dlt):-1:1
    if display_count
        i
    end
    for j=1:size(MSH,1)
        if isempty(MSH(j,1));continue;end
        %The nodes in the main mesh which do not have any connection with
        %other nodes in the main mesh may have connections on other lower
        %dimension elements. However these elements must be removed
        for k=1:size(MSH(j,1).elem,1)
            [ii, jj]=find(MSH(j,1).elem(k,1).id==dlt(i));
            MSH(j,1).elem(k,1).id(ii,:)=[];
            %MSH(j,1).elem(k,1).dom(ii,:)=[];
            ind=find(MSH(j,1).elem(k,1).id>dlt(i));
            MSH(j,1).elem(k,1).id(ind)=MSH(j,1).elem(k,1).id(ind)-1;
        end
    end
end


    

