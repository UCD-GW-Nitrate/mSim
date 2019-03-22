function root_folder = msim_root()
% root_folder = msim_root()
%
% Returns the path of the root folder of the msim. This is the path where
% all the subfolders mFlow, mTrans, etc are located plus the file separator
%
% Version : 1.0
% Author : George Kourakos
% email: giorgk@gmail.com
% web : http://groundwater.ucdavis.edu/msim
% Date 22-Mar_2019 
% Department of Land Air and Water
% University of California Davis

root_folder = which('msim_root');
root_folder = root_folder(1:end - length('mUtil/msim_root.m'));

end

