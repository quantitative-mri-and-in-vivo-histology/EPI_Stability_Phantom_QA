function FJmodified_QualityCheck_23v1204(params)

%	To perform a quantitation of snr, sfnr, stability and drift
%	includes weisskoff plot  MRM 36:643 (1996)
%
%	rev 0	3/3/00		original from noiseave and imgroi
%	rev 1	3/29/02		fix a few header things
%	rev 2	9/4/02		add weissnoise plot
%	rev 3	1/28/03		add phase drift plot 
%				.freq image is scaled 10x
%	rev 4	4/1/03		for fbirn

%	acq: 35	slice 64x64 TR 3.0 TE 30 (3T) 40 (1.5T) 200 frames
%	grecons -s 18 -e 18 Pxxxxx

%FJv23(25.01): Modifications goes as follows: 
%1. Explicit input/output folders were removed. They are based from where
%the code is being executed + the date and time of execution.
%2. Reads now the entire folder with DICOM images (assuming that comes from
%the same measurement). Based on the number of files, it will format
%accordingly.
%3. Small search (using the while function) to obtain the minimum row value
%to generate the corresponding evaluation contour.
%
%FJv23(05.05): Addition of a manual "ROI" selection in case the automatic
%contour is not capable to generate a correct ROI.
%FJv23(16.05): All ROIs become circular; therefore the steps that use rectangular ROIs
%will be changed.

%more off; close all; clear; clc;

if params.defaults
    checkROI = 0;
    DataIsNii = params.isNii;
    fpath_output = fullfile(params.dir_to_analyse,'Results');
    if ~exist(fpath_output,'dir')
        mkdir(fpath_output);
    end
    TR = 1.2;
    TotalSlices = 16;
    slice_eval_start = 1;
    slice_eval_stop = TotalSlices;
    Measurements = [1 250];
else
    checkROI = input('Check ROI per slice: [y/n] (default: y) ', 's');
    
    if isempty(checkROI) || strcmp(checkROI,'y')
        checkROI = 1;
    else
        checkROI = 0;
    end
    
    DataIsNii = input('Is data to be analysed Nii? [y/n] (default: n) ', 's');
    
    if isempty(DataIsNii) || strcmp(DataIsNii,'n')
        DataIsNii = 0;
    else
        DataIsNii = 1;
    end
    
    fpath_output_option = input('Where to save the outputs: DICOM Input (I) or Main Code (M) directory? [I/M]: ','s');

    switch strcmp(fpath_output_option,'I')
        case 1
            fpath_output = fullfile(params.dir_to_analyse,'Results');
            if ~exist(fpath_output,'dir')
                mkdir(fpath_output);
            end
        case 0
            fpath_output = fullfile(pwd,'Results');
            if ~exist(fpath_output,'dir')
                mkdir(fpath_output);
            end
    end
end

if ~isempty(dir(fullfile(params.dir_to_analyse,'loaded_images.mat'))) 
    load(fullfile(params.dir_to_analyse,'loaded_images.mat'));
    disp('Loaded images file exists. This will be used across this analysis');
    
    %load header information
    load(fullfile(params.dir_to_analyse,'header_information.mat'));
    disp('Header file exists. This information will be used across this analysis');
else
    if ~DataIsNii
        disp('Images have not been loaded before from DICOM - therefore the DICOMS will be loaded instead');
        
        % Info extracted from DICOM header if available
        DICOMS_dir = dir(params.dir_to_analyse);
        DICOMS_filelist = {DICOMS_dir(~[DICOMS_dir.isdir]).name};
        fname = DICOMS_filelist{1};
        header = dicominfo(fullfile(params.dir_to_analyse,fname));
        
        if contains(header.ImageType,'MOSAIC')
        %In case of mosaic, do this here:
            loaded_images = mosaic_images_to_data(fullfile(params.dir_to_analyse,DICOMS_filelist));
        else
        %In case of single slices, do this here:
            loaded_images = single_slice_to_data(fullfile(params.dir_to_analyse,DICOMS_filelist));
        end
        save(fullfile(params.dir_to_analyse,'loaded_images.mat'),'loaded_images');
        save(fullfile(params.dir_to_analyse,'header_information.mat'),'header');
    else
        disp('Images are from a Nii file, so this need to be chose before proceeding');
        nii_file = uigetfile();
        header = niftiinfo(fullfile(params.dir_to_analyse,nii_file));
        loaded_images = double(niftiread(header));
        save(fullfile(params.dir_to_analyse,'loaded_images.mat'),'loaded_images');
        save(fullfile(params.dir_to_analyse,'header_information.mat'),'header');
    end
end

% some defaults
I1 = 1;		% first image
I2 = size(loaded_images,4);%500;	% last image
numwin = 4;

NPIX = size(loaded_images,1);%header.AcquisitionMatrix(1);  % No. of Pix Read

if ~params.defaults
    if isfield(header,'RepetitionTime')
        TR = header.RepetitionTime/1000;		% rep time, s
    else
        TR = str2num(input('Repetition time not found. Probably an enhanced DICOM. Insert the TR in seconds: ', 's'));
    end
end

%ROI definitions
if(NPIX == 128)     % H? ??
  R = 30;			% ROI width
else
  R = 15;			% probably 64x64
end
r1 = 1;
r2 = R;

theta = 0:0.01:2*pi;
Xcir = round(R/2*cos(theta) + R/2);
Ycir = round(R/2*sin(theta) + R/2);

%  set up input and ROI mask
mask = poly2mask(Xcir,Ycir,R,R);%ones(R);
npx = sum(mask(:)==1);
%img = zeros(header.AcquisitionMatrix(4),header.AcquisitionMatrix(1),1);
    
% Mosaic Dimensionen ********Nur wenn phase H/F?********
%Mosaic_rows = header.Height / header.AcquisitionMatrix(4);
%Mosaic_columns = header.Width / header.AcquisitionMatrix(1);

% Anzahl Rep.
if ~params.defaults
str = sprintf('Indicate range of measurements to be analysed [min max] (default = [%d %d]) = ', I1, I2);
foo = input(str);
if(~isempty(foo))
    Measurements(1) = foo(1); Measurements(2) = foo(2);
else
    Measurements(1) = I1;
    Measurements(2) = I2;
end

% Schichtanzahl
%slices = 20;
%-FJv23(25.01): Ok, I am not completely sure if this is the same across all
%DICOMs in Siemens, but it seems that the number of slices is given by a
%private variable.
%if ~isfield(header,'Private_0019_100a')
%    slices = double(input('DICOM flag for slice number is not incorporated. Please add it manually: '));
    %slices = Mosaic_columns * Mosaic_rows;
%else
%    slices = double(header.Private_0019_100a);
%end
TotalSlices = size(loaded_images,3);
str = sprintf('Indicate the number of slices in your measurement (default = %d) = ', TotalSlices);
foo = input(str);
if(~isempty(foo))
  TotalSlices = foo(1);
end

% Schichtauswahl f???r eval
slice_eval_start = 1; slice_eval_stop = TotalSlices;
if ~(slice_eval_start == slice_eval_stop) %In case of 1 slice scenario, then it will not require user input 
    str = sprintf('Indicate the number of slices to be analysed [min max] (cr =[%d %d]) = ', slice_eval_start, slice_eval_stop);
    foo = input(str);
    if(~isempty(foo))
        slice_eval_start = foo(1); slice_eval_stop = foo(2);
    end
end
end

N = Measurements(2) - Measurements(1) + 1;                % num time frames
M = r2 - r1 + 1;                % num ROI's

%
% All images are already loaded in loaded_images variable. Therefore, many
% steps are skipped or already defined based on the loaded_images variable.
%
for slice_eval = slice_eval_start:slice_eval_stop
    
    %Einzelbildkoordinaten im Mosaic
    %imgX = fix((slice_eval -1) / double(Mosaic_columns)) * header.AcquisitionMatrix(4) + 1;
    %imgY = mod(slice_eval - 1,Mosaic_columns) * header.AcquisitionMatrix(1) +1;
        
% init ROIs    
    roir = zeros(N,M);

%  begin loop through images

    %Iodd  = zeros(header.AcquisitionMatrix(4),header.AcquisitionMatrix(1),1); % Nur f???r ein Einzelbild im Mosaic
    %Ieven = zeros(header.AcquisitionMatrix(4),header.AcquisitionMatrix(1),1);
    %Syy   = zeros(header.AcquisitionMatrix(4),header.AcquisitionMatrix(1),1); 
    %Syt   = zeros(header.AcquisitionMatrix(4),header.AcquisitionMatrix(1),1);
    Iodd = zeros(size(loaded_images,[1,2]));
    Ieven = Iodd;
    Syy = Iodd;
    Syt = Iodd;
    St = 0;
    Stt = 0;
    S0 = 0;

% Phantomposition im ersten Bild finden
    %dummy = double(dicomread(fullfile(fpath,fname))); % this is the fname from the beginning.
    %FitI = dummy(imgX:imgX+header.AcquisitionMatrix(4)-1,imgY:imgY+header.AcquisitionMatrix(1)-1)/4096; % Nur Einzelbild skaliert auf 0-1
    FitI = loaded_images(:,:,slice_eval,1)/max(max(loaded_images(:,:,slice_eval,1)));
    threshold = graythresh(FitI);
    BW = im2bw(FitI,threshold);
    dim = size(BW);
    %-FJv23(26.01): Ok, this value is weird... or it is not clear what it
    %is for. Therefore, I will just go for a small weird approach:
%    col = 86;
    col = 86;%double(header.AcquisitionMatrix(1));%80;
    row = [];
    while isempty(row)
        row = min(find(BW(:,col)));
        col = col - 1;
    end
    col = col + 1; %-FJv23(26.01): Since it was forced to remove one unit before leaving the while function
                   % then here must be incorporated. Otherwise, the contour
                   % gets empty.
    connectivity = 8;
    num_points   = 180;
    contour = bwtraceboundary(BW, [row, col], 'N', connectivity, num_points,'clockwise');
    x = contour(:,2);
    y = contour(:,1);

    % solve for parameters a, b, and c in the least-squares sense by
    % using the backslash operator
    abc=[x y ones(length(x),1)]\(-(x.^2+y.^2));
    a = abc(1); b = abc(2); c = abc(3);

    % calculate the location of the center and the radius
    xc = -a/2;
    yc = -b/2;
    radius  =  sqrt((xc^2+yc^2)-c);
    
    xc=round(xc);
    yc=round(yc);
    
    Xfit = radius*cos(theta) + xc;
    Yfit = radius*sin(theta) + yc;
      
    % ROI Koordinaten
     %X1 = xc - fix(R/2);
     %X2 = X1 + R -1; 
     %Y1 = yc - fix(R/2); 
     %Y2 = Y1 + R - 1; 
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of "is this ROI ok or should I make it myself?"
    if checkROI      
        figure;
        imshow(FitI); hold on;
        plot(Xfit, Yfit,'LineWidth',2);
        hold off;
        
        str_option = input('Is the ROI selecting properly the phantom? [y/n]: ','s');
        close;
        
        while strcmp(str_option,'n')
            % Draw a circle using the drawcircle function
            figure;
            imshow(FitI); hold on;
            h = drawcircle;

            % Manually select any point on the circumference of the circle
            [xcirc, ycirc] = ginput(1);
            % And update information for further analysis            
            xc = round(h.Center(1));
            yc = round(h.Center(2)); % Get center from drawcircle output
            radius = round(sqrt((xcirc-xc)^2 + (ycirc-yc)^2)); % Calculate radius
            %X1 = round(xc - fix(R/2));
            %X2 = round(X1 + R - 1); 
            %Y1 = round(yc - fix(R/2)); 
            %Y2 = round(Y1 + R - 1); 
            str_option = input('Is the new ROI selecting properly the phantom? [y/n]: ','s');
            Xfit = radius*cos(theta) + xc;
            Yfit = radius*sin(theta) + yc;
            close;
        end    
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    circle_roi(:,:,slice_eval) = poly2mask(Xfit,Yfit,size(loaded_images,1), size(loaded_images,2));
    
    for j = Measurements(1):Measurements(2)
        %fname_second = DICOMS_filelist{j};
        %dummy = double(dicomread(fullfile(fpath,fname_second)));
        %I = dummy(imgX:imgX+header.AcquisitionMatrix(4)-1,imgY:imgY+header.AcquisitionMatrix(1)-1); % Nur Einzelbild
        I = loaded_images(:,:,slice_eval,j);
        if(mod(j,2)==1)
          Iodd = Iodd + I;
        else
          Ieven = Ieven + I;
        end
        
        Syt = Syt + I*j;
        Syy = Syy + I.*I;

        S0 = S0 + 1;
        St = St + j;
        Stt = Stt + j*j;
        roi(S0) = sum(I(squeeze(circle_roi(:,:,slice_eval))))/npx;
        
        for r = r1:r2                 % each roi size
          sub_x = r*cos(theta) + xc;
          sub_y = r*sin(theta) + yc;
          sub_circle_roi = poly2mask(sub_x,sub_y,size(loaded_images,1), size(loaded_images,2));
          %x1 = xc - fix(r/2);   
          %x2 = x1 + r - 1; % Reihe ende unten
          %y1 = yc - fix(r/2); 
          %y2 = y1 + r - 1; % Spalte ende rechts
          %sub = I(x1:x2,y1:y2);
          roir(j-Measurements(1)+1, r) = mean(I(sub_circle_roi));
        end
    end

%  write out diff image
    [~, ProtocolName] = fileparts(params.dir_to_analyse);
    
    Isub = Iodd - Ieven;
    %sub = Isub(X1:X2,Y1:Y2);
    varI = var(Isub(squeeze(circle_roi(:,:,slice_eval))));
    
    diff_volume(:,:,slice_eval) = Isub;
    diff_fname=fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_' ProtocolName '_DiffIm_Slice_' num2str(slice_eval) '.nave']);
    fout = fopen(diff_fname, 'w');
    fwrite(fout, Isub, 'short');
    fprintf('\nwrite file %s\n', diff_fname);
    fclose(fout);

    %  write out ave image
    Sy = Iodd + Ieven; % MosaicBild
    Iave = Sy/N;
    %sub = Iave(X1:X2,Y1:Y2);
    meanI = mean(Iave(squeeze(circle_roi(:,:,slice_eval))));
    
    mean_volume(:,:,slice_eval) = Iave;
    ave_fname=fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_' ProtocolName '_AveIm_Slice_' num2str(slice_eval) '.ave']);
    fout = fopen(ave_fname, 'w');
    fwrite(fout, Iave, 'short');
    fprintf('write file %s\n', ave_fname);
    fclose(fout);

    % find trend line at + b Einzelbild

    D = (Stt*S0 - St*St);
    a = (Syt*S0 - St*Sy)/D;
    b = (Stt*Sy - St*Syt)/D;

    % make sd image

    Var = Syy + a.*a*Stt +b.*b*S0 + 2*a.*b*St - 2*a.*Syt - 2*b.*Sy;
    Isd = sqrt(Var/(N-1));
    
    sd_volume(:,:,slice_eval) = Isd;
    sd_fname=fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_' ProtocolName '_Sd_Slice_' num2str(slice_eval) '.sd']);
    fout = fopen(sd_fname, 'w');
    fwrite(fout, 10*Isd, 'short');
    fprintf('write file %s\n', sd_fname);
    fclose(fout);

    % make sfnr image

    sfnr = Iave./(Isd + eps);
    %img(:) = sfnr;
    %sub = sfnr(X1:X2,Y1:Y2);
    sfnrI = mean(sfnr(squeeze(circle_roi(:,:,slice_eval))));

    SFNR_volume(:,:,slice_eval) = sfnr;
    sfnr_fname=fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_' ProtocolName '_snfr_Slice_' num2str(slice_eval) '.sfnr']);
    fout = fopen(sfnr_fname, 'w');
    fwrite(fout, 10*sfnr, 'short');
    fprintf('write file %s\n', sfnr_fname);
    fclose(fout);

    snr = meanI/sqrt(varI/N);
    fprintf('\nmean, SNR, SFNR = %5.1f  %5.1f  %5.1f\n', meanI, snr, sfnrI);
    %
    figure (slice_eval);

    %  Do fluctation analysis

    x=(1:N);
    p=polyfit(x,roi,2);
    yfit = polyval(p, x);
    y = roi - yfit;

    subplot(numwin,1,1)
    plot(x,roi,x,yfit);
    xlabel('frame num');
    ylabel('Raw signal');
    grid
    m=mean(roi);
    sd=std(y);
    drift = (yfit(N)-yfit(1))/m;
    title(sprintf('%s-slice%02d-%s   percent fluct (trend removed), drift= %5.2f %5.2f',ProtocolName,slice_eval,header.PerformedProcedureStepStartDate, 100*sd/m, 100*drift));

    fprintf('std, percent fluc, drift = %5.2f  %6.2f %6.2f \n', sd, 100*sd/m, 100*drift);

    z = fft(y);
    fs = 1/TR;
    nf = N/2+1;
    f = 0.5*(1:nf)*fs/nf;
    subplot(numwin,1,2);plot(f, abs(z(1:nf)));grid
    ylabel('spectrum');
    xlabel('frequency, Hz');
    ax = axis;

    text(ax(2)*.2, ax(4)*.8, sprintf('mean, SNR, SFNR = %5.1f  %5.1f  %5.1f', meanI, snr, sfnrI));


    %  now do analysis for each roi size
    t = (1:N);
    for r = r1:r2
      y = roir(:, r)';
      yfit = polyval(polyfit(t, y, 2), t);  % 2nd order trend
      F(r) = std(y - yfit)/mean(yfit);
    end
    rr = (r1:r2);
    F = 100*F;              % percent
    fcalc = F(1)./rr;
    rdc = F(1)/F(r2);	% decorrelation distance

    % write log file
    %-FJv23(06.03): This is converted to a table for easier accesibility
    logname = fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_LogFile.log']);
    logout = fopen(logname, 'a');
    fprintf(logout,'Slice %04d %6.2f %6.2f %5.1f  %5.1f  %5.1f %3.1f \n',slice_eval, 100*sd/m, 100*drift, meanI, snr, sfnrI,rdc );
    fclose(logout);
    
    Array_for_table(slice_eval - slice_eval_start+1,:) = [slice_eval,100*sd/m, 100*drift, meanI, snr, sfnrI,rdc]; 

    % plot
    subplot(numwin,1,3);
    loglog(rr, F, '-x', rr, fcalc, '--');
    grid
    xlabel('ROI full width, pixels');
    ylabel('Relative std, %');
    axis([r1 r2 .01 1]);
    text(6, 0.5, 'solid: meas   dashed: calc');
    text(6, 0.25, sprintf('rdc = %3.1f pixels',rdc));

    subplot(numwin,1,4);
    imshow(FitI);
    hold on;
    plot(contour(:,2),contour(:,1),'g','LineWidth',1);
    % display the calculated center
    plot(xc,yc,'yx','LineWidth',2);

    % plot the entire circle
    %theta = 0:0.01:2*pi;
    % use parametric representation of the circle to obtain coordinates
    % of points on the circle
    %Xfit = radius*cos(theta) + xc;
    %Yfit = radius*sin(theta) + yc;

    plot(Xfit, Yfit,'y','LineWidth',2);

    % plot the ROI
    %plot(xc-fix(R/2):xc+fix(R/2)-1,yc-fix(R/2):yc-fix(R/2),'y','Linewidth',2);
    %plot(xc+fix(R/2)-1:xc+fix(R/2)-1,yc-fix(R/2):yc+fix(R/2)-1,'y','Linewidth',2);
    %plot(xc-fix(R/2):xc+fix(R/2)-1,yc+fix(R/2)-1:yc+fix(R/2)-1,'y','Linewidth',2);
    %plot(xc-fix(R/2):xc-fix(R/2),yc-fix(R/2):yc+fix(R/2)-1,'y','Linewidth',2); 
    %hold off;
    
    %circle_roi(:,:,slice_eval) = poly2mask(Xfit,Yfit,size(loaded_images,1), size(loaded_images,2));
    fig_fname=fullfile(fpath_output,[header.PerformedProcedureStepStartDate '_' ProtocolName '_ImageResults_Slice_' num2str(slice_eval) '.jpg']);
    print ('-djpeg',fig_fname);
end % slice_eval 

%Save table here
tab_varNames = ["Slice","Perc. Fluctuation","Drift","Mean signal","SNR","SFNR","Rdc"];

if exist(fullfile(fpath_output,['TableResults_' ProtocolName '.mat']))
    disp('Table results for this measurement already exists. It will be rewritten.');
    load(fullfile(fpath_output,['TableResults_' ProtocolName '.mat']));
    TableAux = array2table(Array_for_table,'VariableNames',tab_varNames);
    TableResults = TableAux;
    save(fullfile(fpath_output,['TableResults_' ProtocolName '.mat']),'TableResults');
else
    TableResults = array2table(Array_for_table,'VariableNames',tab_varNames);
    save(fullfile(fpath_output,['TableResults_' ProtocolName '.mat']),'TableResults');
end

%Save volumes here (? - check header)
niftiwrite(diff_volume,fullfile(fpath_output,[datestr(today,'yymmdd') '_DifferenceVolume.nii']));
niftiwrite(mean_volume,fullfile(fpath_output,[datestr(today,'yymmdd') '_MeanVolume.nii']));
niftiwrite(sd_volume,fullfile(fpath_output,[datestr(today,'yymmdd') '_SdVolume.nii']));
niftiwrite(SFNR_volume,fullfile(fpath_output,[datestr(today,'yymmdd') '_SFNRVolume.nii']));
niftiwrite(double(circle_roi),fullfile(fpath_output,[datestr(today,'yymmdd') '_CircleROI.nii']));

more on;

end
