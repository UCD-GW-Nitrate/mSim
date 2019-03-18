%% msim_polysplit
% 
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%
% This is equivalent to the matlab polysplit command. The only reason that
% we reimplement this, it's to avoid the dependency on the mapping
% toolbox. In addition as for now the Octave is missing this function
% therefore many functions of mSim couldn't run.
% Version : 1.0
%
% Author : George Kourakos
%
% email: giorgk@gmail.com
%
% web : http://groundwater.ucdavis.edu/msim
%
% Date 18-Mar-2014
%
% Department of Land Air and Water
%
% University of California Davis
%
%% Usage
% [Xs, Ys] = msim_polysplit(X, Y)
%
%%
%% Input
%
% _*X*_: [n x 1] A vector of x coordinates
%
% _*Y*_: [n x 1] A vector of y coordinates
%
%% Output
%
% _*Xs*_: A cell array of segments of the X coordinates split into Nans
%
% _*Ys*_: A cell array of segments of the Y coordinates split into Nans
%
%% Example
% Create 2 vectors X and Y with some nans at the same colums:
X = [8,9,1,9,NaN,1,NaN,NaN,NaN,10,NaN,10,NaN,5,8,1,4,9,8,NaN];
Y = [3,2,7,8,NaN,8,NaN,NaN,NaN, 4,NaN, 0,NaN,4,5,8,3,8,5,NaN];
[Xs, Ys] = msim_polysplit(X, Y)
%%
% 
% <msim_help_main.html | main>   <msim_help_demos.html | Tutorials> 
% <msim_function_categories.html | Functions> <http://www.subsurface.gr | website> |
%