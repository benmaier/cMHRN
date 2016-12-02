% use C++ in the following
mex -setup C++;

% compile functions
mex -v CXXFLAGS='-O3 -std=c++11 -stdlib=libc++ -I./cMHRN -I./matlab' matlab/fast_mhrn.cpp matlab/CastResult.cpp cMHRN/mhrn.cpp cMHRN/Utilities.cpp 

% move compiled functions to new folder
mkdir matlabbuild
movefile('fast_mhrn*','./matlabbuild');

% add path to matlab environment via startup-file in user directory
up = userpath;
startuppath = [up(1:end-1),'/startup.m'];
libpath = [pwd,'/matlabbuild'];

fid = fopen(startuppath, 'at');  % append to possibly existing startup-file
fprintf(fid,['\naddpath(''',libpath,''')\n']);
fclose(fid);

% add path for this session
addpath(libpath);
