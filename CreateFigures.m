function CreateFigures(path_measurements)

for meas_indx = 1:numel(path_measurements)
    load(fullfile(path_measurements{meas_indx},'loaded_images.mat'));

    figure(meas_indx); 
    for i = 1:16
        subplot(4,4,i), imagesc(loaded_images(:,:,i,1)/800,[0,1]); title(['Slice #' num2str(i)]);
    end
    
    clear loaded_images

    figure(numel(path_measurements)+1);
    TableFile = dir(fullfile(path_measurements{meas_indx},'Results','TableResults_*.mat'));
    load(fullfile(path_measurements{meas_indx},'Results', TableFile.name));
    subplot(2,2,1), plot(TableResults.("Perc. Fluctuation"),'o--'); hold on; subplot(2,2,2), plot(TableResults.Drift,'o--'); hold on;
    subplot(2,2,3), plot(TableResults.SFNR,'o--');hold on; subplot(2,2,4), plot(TableResults.Rdc,'o--'); hold on;

    for i = 1:4 
        subplot(2,2,i), grid minor; box off; xlabel('Slice number/position'); 
    end
    subplot(2,2,1), title('Percentage fluctuation (in units)');
    subplot(2,2,2), title('Drift');
    subplot(2,2,3), title('Signal-to-fluctuation-noise ratio (SFNR)');
    subplot(2,2,4), title('Radius of decorrelation');

end