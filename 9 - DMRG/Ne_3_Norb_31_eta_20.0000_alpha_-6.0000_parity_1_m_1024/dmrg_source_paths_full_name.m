%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is called automatically when matlab started to set path variable        
% where matlab searches for DMRG components.                                    
%                                                    
% Default list of directories are in "./dmrg_sourcedirs.txt" and
% additional lists are in "pathdef.." files.
%
% [DMRG_SOURCEDIRS] = DMRG_SOURCE_PATHS( )
% [DMRG_SOURCEDIRS] = DMRG_SOURCE_PATHS( 'PATHDEF1', ... )
% [DMRG_SOURCEDIRS] = DMRG_SOURCE_PATHS( __ , SETPATH )
%
%   PATHDEFx       @char, names of additional pathdef files
%   SETPATH        @logical|char,  true/false switch, do not add the DMRG_SOURCEDIRS 
%       to the path automatically if false (default: true)
%   DMRG_SOURCEDIRS  @cellstr[], list of directories to add the path
%
% Note #1: Directory names in the list file sould be relatives, except
% that starts with "/" or "~".
%
% Note #2: Directory name (in the list file) ending with "/*" will be extended
% to the path with the directory itself and all of its subdirectories.
%
% Note #3: The logical SETPATH argument can be given as string 'true' or 'false' too
% (for command line calls). With this feature e.g. you can get a directory list from valid sources.
%                                                                              
% Note: path could not be obtained in main file qc_1_62 because this           
% command produces a run-time error if the code is compiled.                   
%                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change log:
%   20xx.xx.xx by O. Legeza -- created
%   2015.01.14 by T. Mosoni -- add "/*" and "~" handling, more pathdef files
%   2015.11.06 by T. Mosoni -- add SETPATH option
%   2016.04.18 by T. Mosoni -- disable line numbers on warning messages


function [dmrg_sourcedirs] = dmrg_source_paths(varargin)

localpwd='/home/wigner/DMRG_2022_06';

path(path, ['.',filesep]);  % reset search path to actual directory
setpath = true;
dmrg_sourcedirs = cell(0,1);
save_warn = warning('on','all');  warning('off','backtrace');  % warnings without stack info
fprintf('\n');

for fn = [['.',filesep,'dmrg_sourcedirs.txt'], varargin]  % default + given files
    
    if islogical(fn{:}) || isnumeric(fn{:}) % SETPATH option
        setpath = logical(fn{:});
    elseif any(strcmpi({'true','false'}, fn{:}))  % SETPATH as string (e.g. from command line)
        setpath = strcmpi('true',fn{:});
    else
        
        fprintf('Reading directory list file %s\n', fn{:});
        [fid, errmsg] = fopen(fn{:},'r');
        if fid < 0
            error('Can not open file %s -- %s', fn{:}, errmsg);
        else
            while true
                tline = fgetl(fid); %---reading in a line from the text file
                if ~ischar(tline),  break, end;  %---end of file reached
                
                tline = strtrim(tline);  % remove trailing and ending whitespaces
                if any(tline) && (tline(1) ~= '%') %---avoid comments
                    if (tline(1) ~= filesep) && (tline(1) ~= '~') %---extends to relative path
                       % tline = ['.',filesep, tline]; %#ok<AGROW>
                       tline = [localpwd,filesep, tline]
                    end%if
                    if (length(tline) > 2) && ...
                            (tline(end) == '*') && (tline(end-1) == filesep) %--- "/*" to get all subdirs
                        dmrg_sourcedirs = [dmrg_sourcedirs; ...
                            reshape( regexp(genpath(tline(1:end-2)),pathsep,'split'), [],1) ]; %#ok<AGROW>
                    else
                        dmrg_sourcedirs = [dmrg_sourcedirs; tline]; %#ok<AGROW>
                    end%if "/*"
                end%if tline
                
            end%while fid
            fclose(fid);
        end%if fid
        
    end%if islogical
end%for fn

dmrg_sourcedirs = unique(dmrg_sourcedirs,'stable');  % clear duplicates;
dmrg_sourcedirs = [dmrg_sourcedirs; localpwd];  % clear duplicates;

% dmrg_sourcedirs = dmrg_sourcedirs( ~cellfun(@isempty,dmrg_sourcedirs) );   % drop empty items

if setpath
    %     fprintf('Extending MATLAB path with directory %s\n', dmrg_sourcedirs{:});
    %     addpath(dmrg_sourcedirs{:});
    for tline = dmrg_sourcedirs'
        fprintf('Extending MATLAB path with directory %s\n', tline{:});
        addpath(tline{:});
    end%for tline
end%if setpath

warning(save_warn);   % restore warning mode
end%function


%{
-------- Original routine by O. Legeza -----------------
function [dmrg_sourcedirs] = dmrg_source_paths()

%fid = fopen('dmrg_sourcedirs.txt','r');
fid = fopen('./dmrg_sourcedirs.txt','r');
 path(path, './');
 nn = 0;
 while 1
     tline = fgetl(fid); %---reading in a line from the text file
     if ~ischar(tline), break, end
       if any(tline)
       if tline(1) ~= '%' %---avoid comments
          tline = ['./', tline];
          %path(path, tline);
          fprintf('\n Extending MATLAB path with directory %s', tline); 
          addpath(tline);
          nn = nn + 1;
       end%if
       end%if
 end
 fprintf('\n'); 
fclose(fid);

%---create cell object

dmrg_sourcedirs = cell(nn,1);
nn = 0;
%fid = fopen('/home/legeza/DMRG_PROGRAM/Sourcecode/dmrg_sourcedirs.txt','r');
fid = fopen('./dmrg_sourcedirs.txt','r');
 while 1
    tline = fgetl(fid); %---reading in a line from the text file
    if ~ischar(tline), break, end
      if any(tline)
      if tline(1) ~= '%' %---avoid comments
         nn = nn + 1; 
         dmrg_sourcedirs{nn} = tline;
      end%if
      end%if
 end%while
fclose(fid);


%path
-------- Original routine by O.Legeza -----------------
%}

