%% SCP, concatenation files

clear;
close;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));


load('channel_location_16_10-20_mi');


filename_6a = 'ak6_run1_offlineMIterm_20181603155146.gdf';
[s, h]= sload(filename_6a);
session.data=(s)';
session.fs=h.SampleRate;
session.channels={chanlocs16.labels};
session.Event_type=h.EVENT.TYP;
session.Event_pos=h.EVENT.POS;

filename_6b = 'ak6_run2_offlineMIterm_20181603160414.gdf';
[s, h]= sload(filename_6b);
session.data=cat(2,session.data,s');
session.Event_type=cat(1,session.Event_type,h.EVENT.TYP);
session.Event_pos=cat(1,session.Event_pos,h.EVENT.POS);

% filename_6c = 'ak6_run3_offlineMIterm_20181603162521.gdf';
% [s, h]= sload(filename_6c);
% session.data=cat(2,session.data,s');
% session.data=cat(2,session.data,s');
% 
% filename_6d = 'ak6_run4_offlineMIterm_20181603160052.gdf';
% [s, h]= sload(filename);
% session.data=cat(2,session.data,s');



session.data=session.data(1:16,:); % we delete the last channels 17


%% filter butter
[b,a]=butter(2,[0.1 3]/session.fs/2); %bandpass
data_filter=session.data;
for i=1:size(chanlocs16,2)
    data=session.data(i,:);
    data_filter(i,:)=filter(b,a,data);
end 


%% car filter

medium_channels=mean(data_filter);
signal_car=zeros(size(data_filter',1),size(data_filter',2));


 for i=1:size(data_filter',1)
    signal_car(i,:)=data_filter(:,i)'-medium_channels(1,i);

 end

epoch_baseline=epoch_struct(session,200,0,3);
epoch_MI=epoch_struct(session,400,0,3);

    
for i=1:size(epoch_MI.data,3)
    
    c(:,:,i)=epoch_MI.data(:,:,i)';
end
   
c2=mean(c,3);

 %% plottopo
 plottopo(c2,chanlocs16);
 
 %% DECODER 
 
 EpochTraining.data=cat(3, epoch_baseline.data(:,:,1:30), epoch_MI.data(:,:,1:30));
 EpochTraining.labels=cat(1, epoch_baseline.Event_type(1:30), epoch_MI.Event_type(1:30));
 
 EpochTesting.data=cat(3, epoch_baseline.data(:,:,31:60), epoch_MI.data(:,:,31:60));
 EpochTesting.labels=cat(1, epoch_baseline.Event_type(31:60), epoch_MI.Event_type(31:60));
 
 for i=1:size(chanlocs16,2)
   [PwelchTrainingOnechannel{i},freq_1]=pwelch_decoder(i,epoch_baseline,[4:2:40],epoch_baseline.fs); 
   [PwelchTestingOnechannel{i},freq_2]=pwelch_decoder(i,epoch_MI,[4:2:40],epoch_MI.fs); 
   
%    subplot(4,4,i)
%    plot(freq_1(1:55),10*log10(pwelch_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_MI_onechannel{i}(1:55)))
%    xlabel('Frequency [Hz]');
%    ylabel('PSD [dB]');
%    title(sprintf('spectral density for channels %s',epoch_baseline.channels{1,i}));
%    legend('baseline','MI');
 end 
 
% TAGLIARE FROM 4 to 40 Hz
% interset function to cut at the right 

TrainingData=zeros(60, 19*16);
TestingData=zeros(60, 19*16);

for i=1:60
    
    a=[];
    for j=1:16
        
    TrainingData(i,1:(19*j))=cat(2,a,cell2mat(PwelchTrainingOnechannel{1,i}{1,j}));
    a=TrainingData(i,1:(19*j));
        
    end
    
end
 
 
 
 
 
 
 
 
 
 