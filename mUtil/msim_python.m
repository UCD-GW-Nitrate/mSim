function [result status] = msim_python(varargin)
% This file is a copy of the perl.m file provided by Matlab
% and we have just rename the function.
% @@@@@@@@@@@@ IMPORTANT @@@@@@@@@@@@@
% TO USE THIS SCRIPT YOU NEED TO DEFINE THE python path
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%
% Original file doc below:
%PERL Execute Perl command and return the result.
%   PERL(PERLFILE) calls perl script specified by the file PERLFILE
%   using appropriate perl executable.
%
%   PERL(PERLFILE,ARG1,ARG2,...) passes the arguments ARG1,ARG2,...
%   to the perl script file PERLFILE, and calls it by using appropriate
%   perl executable.
%
%   RESULT=PERL(...) outputs the result of attempted perl call.  If the
%   exit status of perl is not zero, an error will be returned.
%
%   [RESULT,STATUS] = PERL(...) outputs the result of the perl call, and
%   also saves its exit status into variable STATUS. 
% 
%   If the Perl executable is not available, it can be downloaded from:
%     http://www.cpan.org
%
%   See also SYSTEM, JAVA, MEX.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.4.8 $

%%%%%%% Change the following line with your python path
arcpy_path = 'c:\Python27\ArcGIS10.2\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmdString = '';

% Add input to arguments to operating system command to be executed.
% (If an argument refers to a file on the MATLAB path, use full file path.)
for i = 1:nargin
    thisArg = varargin{i};
    if isempty(thisArg) || ~ischar(thisArg)
        error('MATLAB:perl:InputsMustBeStrings', 'All input arguments must be valid strings.');
    end
    if i==1
        if exist(thisArg, 'file')==2
            % This is a valid file on the MATLAB path
            if isempty(dir(thisArg))
                % Not complete file specification
                % - file is not in current directory
                % - OR filename specified without extension
                % ==> get full file path
                thisArg = which(thisArg);
            end
        else
            % First input argument is PerlFile - it must be a valid file
            error('MATLAB:perl:FileNotFound', 'Unable to find Perl file: %s', thisArg);
        end
    end
  
  % Wrap thisArg in double quotes if it contains spaces
  if any(thisArg == ' ')
    thisArg = ['"', thisArg, '"'];
  end
  
  % Add argument to command string
  cmdString = [cmdString, ' ', thisArg];
end

% Execute Perl script
errTxtNoPerl = 'Unable to find Perl executable.';

if isempty(cmdString)
  error('MATLAB:perl:NoPerlCommand', 'No perl command specified');
elseif ispc
  % PC
  %perlCmd = fullfile('c:\Python25\');%With this command you can run python and acrgis 9.3 commands only
  %perlCmd = fullfile('c:\Python26\ArcGIS10.0\');arcpypath
  perlCmd = fullfile(arcpy_path);
  cmdString = ['python' cmdString];	 
  perlCmd = ['set PATH=',perlCmd, ';%PATH%&' cmdString];
  [status, result] = dos(perlCmd);
else
  % UNIX
  [status ignore] = unix('which perl'); %#ok
  if (status == 0)
    cmdString = ['perl', cmdString];
    [status, result] = unix(cmdString);
  else
    error('MATLAB:perl:NoExecutable', errTxtNoPerl);
  end
end

% Check for errors in shell command
if nargout < 2 && status~=0
  error('MATLAB:perl:ExecutionError', ...
        'System error: %sCommand executed: %s', result, cmdString);
end

