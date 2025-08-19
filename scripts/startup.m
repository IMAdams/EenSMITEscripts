%DipImage

% addpath('C:\Program Files\DIPimage 2.9\common\dipimage')
% dip_initialise
% dipsetpref('imagefilepath','C:\Program Files\DIPimage 2.9\images')
% addpath('C:\Program Files\DIPimage 2.9\common\dipimage\demos')



%set up SMITE
% run(fullfile(userpath,'smite','MATLAB','setupSMITE.m'))
% % addpath(genpath('C:\Users\imadams\Documents\MATLAB\smite'));
setupSMITE
warning('off','MATLAB:oldPfileVersion') 
diproot='C:\Program Files\DIPimage 2.9';
run([diproot '\dipstart.m'])
dipsetpref('DefaultMappingMode','lin');
dipsetpref('DefaultFigureHeight', 512,'DefaultFigureWidth', 1024 )
dipsetpref('TrueSize','off')
dipsetpref('ImageFilePath','C:\Program Files\DIPimage 2.9\images')
set(0,'DefaultFigurePaperType','US Letter');
dipfig -unlink

% diproot = 'C:\Program Files\DIPimage 2.9';
% run([diproot '\dipstart.m'])


%Setup Python- only need to run this once 
%install anaconda
% pe = pyenv('Version','C:\Users\imadams\AppData\Local\anaconda3\envs\een_matpy\python.exe') ;


% --- OLD ---


%HSM (requires DipImage)
%addpath([ dipshare 'HSM\tags\RC2\HSM_helpers']);
%addpath([ dipshare 'HSM\tags\RC2\hyperView']);
%addpath([ dipshare 'HSM\tags\RC2\Mex']);

%close and restart matlab. 
%you can now drag and drop *.hsi files into the command window. This will
%load a varialbe called 'hsiObj'. 

%Pop up HyperView:   
%hsiObj.hyperView
%addpath([ dipshare 'SPT\branches\Rev1']);









