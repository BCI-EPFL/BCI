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

 filename_6c = 'ak6_run3_offlineMIterm_20181603162521.gdf';
 [s, h]= sload(filename_6c);
session.data=cat(2,session.data,s');
session.Event_type=cat(1,session.Event_type,h.EVENT.TYP);
session.Event_pos=cat(1,session.Event_pos,h.EVENT.POS);

filename_6d = 'ak6_run4_offlineMIterm_20181603164052.gdf';
session.data=cat(2,session.data,s');
session.Event_type=cat(1,session.Event_type,h.EVENT.TYP);
session.Event_pos=cat(1,session.Event_pos,h.EVENT.POS);


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

session.data=signal_car';
epoch_baseline=epoch_struct(session,200,0,3);
epoch_MI=epoch_struct(session,400,0,3);
epoch_MI_term=epoch_struct(session,555,0,3);
    
for i=1:size(epoch_MI.data,3)
    
    c(:,:,i)=epoch_MI.data(:,:,i)';
end
   
medium_average=mean(c,3);

 %% plottopo
 plottopo(medium_average,chanlocs16);
 
 %% DECODER 
 
 EpochTraining.BaseMI.data=cat(3, epoch_baseline.data(:,:,1:60), epoch_MI.data(:,:,1:60));
 EpochTraining.BaseMI.labels=cat(1, epoch_baseline.Event_type(1:60), epoch_MI.Event_type(1:60));
 
 EpochTesting.BaseMI.data=cat(3, epoch_baseline.data(:,:,61:120), epoch_MI.data(:,:,61:120));
 EpochTesting.BaseMI.labels=cat(1, epoch_baseline.Event_type(61:120), epoch_MI.Event_type(61:120));
 
 for i=1:size(chanlocs16,2)
   [PwelchTrainingOnechannel{i},freq]=pwelch_decoder(i,EpochTraining.BaseMI,[4:2:40],epoch_baseline.fs); 
   [PwelchTestingOnechannel{i},freq]=pwelch_decoder(i,EpochTesting.BaseMI,[4:2:40],epoch_MI.fs); 
   
%    subplot(4,4,i)
%    plot(freq_1(1:55),10*log10(pwelch_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_MI_onechannel{i}(1:55)))
%    xlabel('Frequency [Hz]');
%    ylabel('PSD [dB]');
%    title(sprintf('spectral density for channels %s',epoch_baseline.channels{1,i}));
%    legend('baseline','MI');
 end 
 

% interset function to cut at the right 

TrainingData.BaseMI=zeros(size(epoch_baseline.data,3), size(freq,2)*size(chanlocs16,2));
TestingData.BaseMI=zeros(size(epoch_baseline.data,3), size(freq,2)*size(chanlocs16,2));

for i=1:size(epoch_baseline.data,3)
    
    a=[];
    for j=1:size(chanlocs16,2)
        
    TrainingData.BaseMI(i,1:(size(freq,2)*j))=cat(2,a,(PwelchTrainingOnechannel{1,j}{1,i}));
    a=TrainingData.BaseMI(i,1:(size(freq,2)*j));
        
    end
    
end
 
 
[ind, power_feat] = rankfeat(TrainingData.BaseMI, EpochTraining.BaseMI.labels,  'fisher');
figure(1); clf
plot(power_feat, '-')
xlabel('feature index')
ylabel('Fisher score')
grid minor

PowerTraining.BaseMI=zeros(1,size(power_feat,2));
for i=1:size(PowerTraining.BaseMI,2)
    PowerTraining.BaseMI(i)=power_feat(ind==i); 
end
PowerTraining.BaseMI=(reshape(PowerTraining.BaseMI,[size(freq,2),size(chanlocs16,2)]))';

figure
imagesc([4 40],[1 16],PowerTraining.BaseMI)
title('Fisher scores - baseline vs MI')
 

%decoder Baseline-MI_termination

EpochTraining.BaseMIterm.data=cat(3, epoch_baseline.data(:,:,1:60), epoch_MI_term.data(:,:,1:60));
EpochTraining.BaseMIterm.labels=cat(1, epoch_baseline.Event_type(1:60), epoch_MI_term.Event_type(1:60));
 
 EpochTesting.BaseMIterm.data=cat(3, epoch_baseline.data(:,:,61:120), epoch_MI_term.data(:,:,61:120));
 EpochTesting.BaseMIterm.labels=cat(1, epoch_baseline.Event_type(61:120), epoch_MI_term.Event_type(61:120));
 
 
 clear PwelchTrainingOnechannel;
 clear freq;
 for i=1:size(chanlocs16,2)
   [PwelchTrainingOnechannel{i},freq]=pwelch_decoder(i,EpochTraining.BaseMIterm,[4:2:40],epoch_baseline.fs); 
   [PwelchTestingOnechannel{i},freq]=pwelch_decoder(i,EpochTesting.BaseMIterm,[4:2:40],epoch_MI_term.fs); 
  
 end 
 

TrainingData.BaseMIterm=zeros(size(epoch_baseline.data,3), size(freq,2)*size(chanlocs16,2));
TestingData.BaseMIterm=zeros(size(epoch_baseline.data,3), size(freq,2)*size(chanlocs16,2));

for i=1:size(epoch_baseline.data,3)
    
    a=[];
    for j=1:size(chanlocs16,2)
        
    TrainingData.BaseMIterm(i,1:(size(freq,2)*j))=cat(2,a,(PwelchTrainingOnechannel{1,j}{1,i}));
    a=TrainingData.BaseMIterm(i,1:(size(freq,2)*j));
        
    end
    
end
 
 
[ind, power_feat] = rankfeat(TrainingData.BaseMIterm, EpochTraining.BaseMIterm.labels,  'fisher');
figure(1); clf
plot(power_feat, '-')
xlabel('feature index')
ylabel('Fisher score')
grid minor

PowerTraining.BaseMIterm=zeros(1,size(power_feat,2));
for i=1:size(PowerTraining.BaseMIterm,2)
    PowerTraining.BaseMIterm(i)=power_feat(ind==i); 
end
PowerTraining.BaseMIterm=(reshape(PowerTraining.BaseMIterm,[size(freq,2),size(chanlocs16,2)]))';

figure
imagesc([4 40],[1 16],PowerTraining.BaseMIterm)
title('Fisher scores - baseline vs MI')
 
