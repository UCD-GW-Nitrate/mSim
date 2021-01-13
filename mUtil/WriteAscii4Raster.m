function WriteAscii4Raster(filename,TAB, xl, yl,csz, nodata)
% WriteAscii4Raster(filename,TAB,xl,yl,csz,nodata)
%
% Writes Grid data to an Ascii format suitable for ArcGis import
%
% Input
% filename : The name of the ascii file where the data wil be printed
% TAB      : Data matrix
% xl       : x coordinate of the left lower corner
% yl       : y coordinate of the left lower corner
% csz      : Cell size
% nodata   : Scalar that represents no data value. Do not use nan
%
% Output
% No output on the matlab workspace. Just prints into a file
%
% 
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 22-Mar_2014 
% Department of Land Air and Water
% University of California Davis

if isnan(csz)
    error('Do not use nan as no data value')
end

[nr nc]=size(TAB);
fid=fopen([filename '.asc'],'w');
if fid < 0
    error(['Cannot create the file "' [filename '.asc'] '"' ])
    return;
end
fprintf(fid,'NCOLS %g\n',nc);
fprintf(fid,'NROWS %g\n',nr);
fprintf(fid,'XLLCORNER %f\n',xl);
fprintf(fid,'YLLCORNER %f\n',yl);
fprintf(fid,'CELLSIZE %f\n',csz);
fprintf(fid,'NODATA_VALUE %g\n',nodata);
frmt = '%10.5f ';
for j=1:nc
    frmt = [frmt ' %10.5f'];
end
frmt = [frmt '\n'];
fprintf(fid,frmt ,TAB');

%for i=1:nr
%    for j=1:nc
%        fprintf(fid,'%10.5f ',TAB(i,j));
%    end
%    fprintf(fid,'\n');
%end
fclose(fid);