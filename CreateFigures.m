function CreateFigures(path_measurements)

for meas_indx = 1:numel(path_measurements)
    load(fullfile(path_measurements{meas_indx},'loaded_images.mat'));
    ROImask = dir(fullfile(path_measurements{meas_indx},'Results','*_CircleROI.nii'));
    
    if numel(ROImask) > 1
        ROImask_volume = niftiread(fullfile(ROImask(end).folder,ROImask(end).name));
    else
        ROImask_volume = niftiread(fullfile(ROImask.folder,ROImask.name));
    end
        
    figure(meas_indx); 
    for i = 1:16
        subplot(4,4,i), imagesc(loaded_images(:,:,i,1)/800,[0,1]); title(['Slice #' num2str(i)]); hold on;       
        subplot(4,4,i), imcontour(ROImask_volume(:,:,i),1,'-r');
    end
    
    clear loaded_images ROImask ROImask_volume 

    figure(numel(path_measurements)+1);
    TableFile = dir(fullfile(path_measurements{meas_indx},'Results','TableResults_*.mat'));
    load(fullfile(path_measurements{meas_indx},'Results', TableFile.name));
    subplot(2,2,1), plot(TableResults.Slice, TableResults.("Perc. Fluctuation"),'o--'); hold on; 
    subplot(2,2,2), plot(TableResults.Slice, TableResults.Drift,'o--'); hold on;
    subplot(2,2,3), plot(TableResults.Slice, TableResults.SFNR,'o--');hold on; 
    subplot(2,2,4), plot(TableResults.Slice, TableResults.Rdc,'o--'); hold on;

    for i = 1:4 
        subplot(2,2,i), grid minor; box off; xlabel('Slice number/position'); 
    end
    subplot(2,2,1), title('Percentage fluctuation (in units)');
    subplot(2,2,2), title('Drift');
    subplot(2,2,3), title('Signal-to-fluctuation-noise ratio (SFNR)');
    subplot(2,2,4), title('Radius of decorrelation');
    
    clear TableResults

end