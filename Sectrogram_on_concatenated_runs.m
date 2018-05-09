clear;
close;
clc

x=input('Enter subject: 1 for Elisabetta, 2 for Marco, 3 for Giammarco \n');

%% INPUTS
%Defining the windown for epoching

%--Baseline
TimeBeforeEventBaseline=-3;
TimeAfterEventBaseline=3;
%--MI
TimeBeforeEventMI=-2;
TimeAfterEventMI=6;
%--MI termination
TimeBeforeEventMItermination=-3;
TimeAfterEventMItermination=3;

%--Cyclic frequency for spectrogram
Cyclic_freq=[5:0.1:40]; %inital:resolution:final


%% Loading paths and files

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));
addpath(genpath('folder_runs_ak6'));

load('channel_location_16_10-20_mi');

FilesEli={'ak6_run1_offlineMIterm_20181603155146.gdf', 'ak6_run2_offlineMIterm_20181603160414.gdf', 'ak6_run3_offlineMIterm_20181603162521.gdf','ak6_run4_offlineMIterm_20181603164052.gdf'};
FilesMarco={'ak4_run1_offlineMIterm_20182103093936.gdf', 'ak4_run2_offlineMIterm_20182103095132.gdf', 'ak4_run3_offlineMIterm_20182103101438.gdf', 'ak4_run4_offlineMIterm_20182103102802.gdf'};
FilesGiam={'ak5_run1_offlineMIterm_20182003154933.gdf', 'ak5_run2_offlineMIterm_20182003160110.gdf','ak5_run3_offlineMIterm_20182003163100.gdf', 'ak5_run4_offlineMIterm_20182003164443.gdf'};

AllFiles=[FilesEli; FilesMarco; FilesGiam];


for j=1:length(FilesEli)

filename = (AllFiles{x,j});
[s, h]= sload(filename);


% Creating the main structure from data

session.fs=h.SampleRate;
session.data=(s)';
session.channels={chanlocs16.labels};
session.Event_type=h.EVENT.TYP;
session.Event_pos=h.EVENT.POS;

% spatial filtering on the raw data.CAR

    medium_channels=mean(s');
    signal_car=zeros(size(s,1),size(s,2));
    for i=1:size(s,1)
        signal_car(i,:)=s(i,:)-medium_channels(1,i);
        
    end

    session_filt_CAR=session;
    session_filt_CAR.data=signal_car';
    
% p_welch on car data

    filt_epoch_baseline_CAR=epoch_struct(session_filt_CAR,200,TimeBeforeEventBaseline,TimeAfterEventBaseline);
    filt_epoch_MI_CAR=epoch_struct(session_filt_CAR,400,TimeBeforeEventMI,TimeAfterEventMI);
    filt_epoch_MI_CAR_termination=epoch_struct(session_filt_CAR,555,TimeBeforeEventMItermination,TimeAfterEventMItermination);
    
    for i=1:size(chanlocs16,2)
        [pwelch_car_bas_onechannel{i},freq_1]=pwelch_for_each_channel(i,filt_epoch_baseline_CAR,500,session_filt_CAR.fs);
        [pwelch_car_MI_onechannel{i},freq_2]=pwelch_for_each_channel(i,filt_epoch_MI_CAR,500,session_filt_CAR.fs);
        
        if j==1
        
        subplot(4,4,i)
        plot(freq_1(1:55),10*log10(pwelch_car_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_car_MI_onechannel{i}(1:55)))
        xlabel('Frequency [Hz]');
        ylabel('PSD [dB]');
        title(sprintf('SPD comparison for channels %s',filt_epoch_baseline_CAR.channels{1,i}));
        legend('baseline','MI');
        
        end
    end

% Spectrogram CAR MI termination

for i=1:16
   
   channels={'FZ','FC3','FC1','FCz','FC2','FC4','C3','C1', 'Cz', 'C2', 'C4' ,'CP3' ,'CP1' ,'CPZ', 'CP2' ,'CP4'}; 
 
   [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_CAR, filt_epoch_MI_CAR_termination, i, filt_epoch_baseline_CAR.fs, filt_epoch_baseline_CAR.fs-32, Cyclic_freq);
   
   if j==1
    Matrix.(channels{i})(1:size(spect_for_one_channel,1), 1:size(spect_for_one_channel,2))=0;
   end
   
   if i==1
       time=t;
   end
       
   Matrix.(channels{i})=Matrix.(channels{i}) + spect_for_one_channel;
   flag=flag+1;
end

end

figure;

for i=1:16
   
   Matrix.(channels{i})=Matrix.(channels{i})./length(FilesEli);
   subplot(4,4,i);
   time=t(1,:)+ filt_epoch_MI_CAR_termination.time;
   imagesc('XData',time,'YData',f,'CData', 10*log(Matrix.(channels{i})));
   axis tight;
   colorbar;
   % set(h,'ylim',[-20 20]);
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   title((sprintf('Channel %s',filt_epoch_baseline_CAR.channels{1,i})));
  
end


% suptitle('Spectrogram CAR filtered data - MI termination');


