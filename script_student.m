 clear;
close;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));
addpath(genpath('codeProject1'));

load('channel_location_16_10-20_mi');

folderName =  'folder_runs_ak6';

params_spectrogram.mlength    = 1;
params_spectrogram.wlength    = 0.5;
params_spectrogram.pshift     = 0.25;                  
params_spectrogram.wshift     = 0.0625;              
params_spectrogram.selchans   = 1:16;     
params_spectrogram.freq = 4:0.1:40;

% File Processing
subSessionFolders = dir([folderName filesep '*.gdf']);
nFile = numel(subSessionFolders);


for iFile = 1:nFile
    disp('****************************')
    fprintf('******* Run %d/%d ******* \n',iFile,nFile)
    fileName = fullfile(folderName,subSessionFolders(iFile).name);
    
    %% Extract data
    % Load the data and put it into a structure as before
    [signal,header] = sload(fileName);
    session.data=signal;
    session.rate=512;
    session.event.name=header.EVENT.TYP;
    session.event.position=header.EVENT.POS;
    % here you put your structure function for sessioninstead of create_your_own_structure()
    s{iFile}  = session;
end
%% BUTTER FILTER

[b,a]=butter(2,[5 40]/session.rate/2); %bandpass
for j=1:numel(s)
    data_filter{j}=s{j};
    
    for i=1:size(chanlocs16,2)
        
        data=s{j}.data(:,i);
        data_filter{j}.data(:,i)=filter(b,a,data);
    end
    clear data;
   
end
%% CAR FILTER

for j=1:numel(data_filter)
 medium_channels=mean(data_filter{j}.data');
    signal_car{j}=data_filter{j};
    for i=1:size(data_filter{j}.data,1)
        signal_car{j}.data(i,:)=data_filter{j}.data(i,:)-medium_channels(1,i);
        
    end
end
%% Extract PSD
s = preprocess_spectrogram(s,params_spectrogram);
% return session with data in 3D (a new dimension for frequency!)


%% Epoching
% do the epoching you will need to adapt and create a new fuction maybe (keep the last one!)
% + concantenate runs 
%epochs_PSD_ONSET  = epoching_function();

for i=1:numel(s)
struct_epoch=epoch_window(s{i},Align_Event,timebEvent,timeaEvent,windowlength,overlap)
end