%% Main code to execute QA EPI Stability.
more off; close all; clear; clc;
% The directory selected MUST be the parent directory. It means, that the
% main directory must contain the First, Second and Third EPI measurements,
% and in each folder must contain the Dicoms and/or Nii. 
main_dir = dir(uigetdir());
dir_list = {main_dir([main_dir.isdir] & ~strcmp({main_dir.name},'.') & ~strcmp({main_dir.name},'..')).name};
full_path_dir = fullfile(main_dir(1).folder, dir_list);

isNii = input('Is the data DICOM or Nii? [D/I] (default: D): ','s');

if isempty(isNii) || isNii == 'D'
    disp('Data is DICOM');
    params.isNii = 0;
else
    disp('Data is Nii');
    params.isNii = 1;
end

execute_defaults = input('Do you want to execute the defaults? [y/n] (default: y): ','s');

if isempty(execute_defaults) || execute_defaults == 'y'
    fprintf(['Defaults will be executed. It means: \n' ...
            '1. CheckROI will be deactivated \n' ...
            '2. Results are saved in each measurement folder, not in the code folder \n' ...
            '3. Number of slices is 16 \n' ...
            '4. Repetition time is 1.2 seconds \n' ...
            '5. All slices across all repetitions will be analysed. \n' ...
            '6. These parameters are predefined even if Nii data is being read. \n' ...
            'If any of these defaults do not match with your acquisition, please run the code without executing the defaults \n']);
    params.defaults = 1;
    pause(2);
else
    params.defaults = 0;
    disp('Defaults will be not executed, please check the command window during execution');
end


for dir_indx = 1:numel(full_path_dir)
    disp(['Analysing data (' num2str(dir_indx) '/' num2str(numel(full_path_dir)) ') located in: ' ...
          fullfile(full_path_dir{dir_indx})]);
    pause(2);
    params.dir_to_analyse = full_path_dir{dir_indx};
    FJmodified_QualityCheck_23v1204(params);
    close all; 
    more off;
end

CreateFigures(full_path_dir);