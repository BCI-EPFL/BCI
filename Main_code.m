clear;
close;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));

load('channel_location_16_10-20_mi');

filename = 's_run1_offlineMIterm_20180703154501.gdf';
[s, h]= sload(filename);



session.fs=h.SampleRate;
session.data=(s)';
session.channels={chanlocs16.labels};
session.Event_type=h.EVENT.TYP;
session.Event_pos=h.EVENT.POS;

%% pwelch
epoch_baseline=epoch_struct(session,200,0,3);
epoch_MI=epoch_struct(session,400,0,3);

for i=1:16
   [pwelch_bas_onechannel,freq]=pwelch_for_each_channel(i,epoch_baseline,500,epoch_baseline.fs); 
   [pwelch_MI_onechannel,freq]=pwelch_for_each_channel(i,epoch_MI,500,epoch_MI.fs); 
   figure
   plot(freq,10*log10(pwelch_bas_onechannel),freq,10*log10(pwelch_MI_onechannel))
   xlabel('Frequency [Hz]');
   ylabel('Power Spectral Density [dB]');
end    

%% temporal filtering on the raw data
[b,a]=butter(2,[5 40]/512/2); %bandpass
data_filter=zeros(16,330752);

for i=1:size(chanlocs16,2)
    data=session.data(i,:);
    data_filter(i,:)=filter(b,a,data);
end 

session_filt=session;
session_filt.data=data_filter;
session_filt.data(17,:)=zeros(1,length(session_filt.data(2,:)));
filt_epoch_baseline=epoch_struct(session_filt,200,0,3);
filt_epoch_MI=epoch_struct(session_filt,400,0,3);

for i=1:16
   [filt_pwelch_bas_onechannel,freq]=pwelch_for_each_channel(i,filt_epoch_baseline,500,filt_epoch_baseline.fs); 
   [filt_pwelch_MI_onechannel,freq]=pwelch_for_each_channel(i,filt_epoch_MI,500,filt_epoch_MI.fs); 
   figure
   plot(freq,10*log10(pwelch_bas_onechannel),freq,10*log10(pwelch_MI_onechannel),freq,10*log10(filt_pwelch_bas_onechannel),freq,10*log10(filt_pwelch_MI_onechannel))
   xlabel('Frequency [Hz]');
   ylabel('Power Spectral Density [dB]');
end    


%% sdp after temporal filtering

session.data=data_filter;
baseline=epoch_struct(session, 200, 0,3);
MI=epoch_struct(session, 400, 0,3);


for i=1:size(chanlocs16,2)
[pwr_b, fre_b]=pwelch(squeeze(epoch_baseline.data(:,i,:)),0.5*session.fs,0.5*0.5*session.fs,500,session.fs);

plot(fre_b,10*log(pwr_b));%10log
hold on
legend('baseline channel 1','baseline channel 2','baseline channel 3','baseline channel 4','baseline channel 5','baseline channel 6','baseline channel 7','baseline channel 8','baseline channel 9','baseline channel 10','baseline channel 11','baseline channel 12','baseline channel 13','baseline channel 14','baseline channel 15','baseline channel 16');
xlabel('frequencies');
ylabel('spectral density');
title('spectral density');
end 
%%
for i=1:size(chanlocs16,2)
[pwr_MI,fre_MI]=pwelch(squeeze(MI.data(:,i,:)),0.5*session.fs,0.5*0.5*session.fs,500,session.fs);
plot(fre_MI,pwr_MI);
hold on 
legend('motor imagery channel 1','motor imagery channel 2','motor imagery channel 3','motor imagery channel 4','motor imagery channel 5','motor imagery channel 6','motor imagery channel 7','motor imagery channel 8','motor imagery channel 9','motor imagery channel 10','motor imagery channel 11','motor imagery channel 12','motor imagery channel 13','motor imagery channel 14','motor imagery channel 15','motor imagery channel 16');
xlabel('frequencies');
ylabel('spectral density');
title('spectral density');
end 


%% spatial filtering on the raw data.CAR
medium=mean(mean(s));
signal_car=zeros(size(s,1),size(s,2));
for i=1:size(s,2)-1
    signal_car(:,i)=s(:,1)-medium;
end

plot(signal_car(:,9))
hold on
plot(s(:,9))
title ('car filter and raw signal');
labels('CAR signal','raw signal');

%%laplacian filter
laplacian=s(:,1:16)*lap;
subplot(3,1,1)
plot(s(:,9))
subplot(3,1,2)
plot(laplacian(:,9))
subplot(3,1,3)
plot(signal_car(:,9))
  

%topoplot






















