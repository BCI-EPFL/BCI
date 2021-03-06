clear;
close;
clc

%% INPUTS
%Defining the windown for epoching

%--Baseline
TimeBeforeEventBaseline=0;
TimeAfterEventBaseline=3;
%--MI
TimeBeforeEventMI=-2;
TimeAfterEventMI=6;
%--MI termination
TimeBeforeEventMItermination=-3;
TimeAfterEventMItermination=3;

%--Cyclic frequency for spectrogram
Cyclic_freq=[5:0.1:40]; %inital:resolution:final

%% SELECT OPTIONS BY UNCOMMENTING 
Action=[0, 0, 0, 0, 0];

% Action(1)=1;              % Pwelch on raw data

% Action(1)=2;              % Temporal filtering on raw data

% Action(1)=3;              % CAR on raw data

% Action=[3, 1, 0, 0, 0]    % Pwelch on CAR data

% Action=[3, 1, 2, 0, 0]    % Plot_comparison CAR_pwelch VS RAW_pwelch

% Action(1)=4;              % Laplacian filter

% Action=[1, 0, 0, 1, 0]    % Spectrogram MI init/term (on raw Pwelch)

 Action=[3, 1, 0, 2, 0]    % Spectrogram CAR MI init/term

% Action=[4, 0, 0, 3, 0]    % Spectrogram Laplacian MI init/term

% Action=[1, 0, 0, 0, 1]    % Topoplot raw

% Action=[3, 1, 0, 0, 2]    % Topoplot CAR

% Action=[4, 0, 0, 0, 3]    % Topoplot Laplacian

%% Loadin paths and files

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));
addpath(genpath('folder_runs_ak6'));

load('channel_location_16_10-20_mi');

filename = 'ak6_run2_offlineMIterm_20181603160414.gdf';
[s, h]= sload(filename);

%% Creating the main structure from data

session.fs=h.SampleRate;
session.data=(s)';
session.channels={chanlocs16.labels};
session.Event_type=h.EVENT.TYP;
session.Event_pos=h.EVENT.POS;

%% Pwelch on raw data
if (Action(1)==1)
    epoch_baseline=epoch_struct(session,200,TimeBeforeEventBaseline,TimeAfterEventBaseline);
    epoch_MI=epoch_struct(session,400,TimeBeforeEventMI,TimeAfterEventMI);
    epoch_MI_termination=epoch_struct(session,555,TimeBeforeEventMItermination,TimeAfterEventMItermination);
    
    for i=1:size(chanlocs16,2)
        [pwelch_bas_onechannel{i},freq_1]=pwelch_for_each_channel(i,epoch_baseline,500,epoch_baseline.fs);
        [pwelch_MI_onechannel{i},freq_2]=pwelch_for_each_channel(i,epoch_MI,500,epoch_MI.fs);
        
        subplot(4,4,i)
        plot(freq_1(1:55),10*log10(pwelch_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_MI_onechannel{i}(1:55)))
        xlabel('Frequency [Hz]');
        ylabel('PSD [dB]');
        title(sprintf('spectral density for channels %s',epoch_baseline.channels{1,i}));
        legend('baseline','MI');
    end
end
%% temporal filtering on the raw data
if Action(1)==2
    [b,a]=butter(2,[5 40]/session.fs/2); %bandpass
    data_filter=session.data;
    for i=1:size(chanlocs16,2)
        data=session.data(i,:);
        data_filter(i,:)=filter(b,a,data);
    end
    
    session_filt=session;
    session_filt.data=data_filter;
    filt_epoch_baseline=epoch_struct(session_filt,200,TimeBeforeEventBaseline,TimeAfterEventBaseline);
    filt_epoch_MI=epoch_struct(session_filt,400,TimeBeforeEventMI,TimeAfterEventMI);
    filt_epoch_MI_termination=epoch_struct(session_filt,555,TimeBeforeEventMItermination,TimeAfterEventMItermination);
    
    for i=1:size(chanlocs16,2)
        [filt_pwelch_bas_onechannel,freq_1]=pwelch_for_each_channel(i,filt_epoch_baseline,500,filt_epoch_baseline.fs);
        [filt_pwelch_MI_onechannel,freq_2]=pwelch_for_each_channel(i,filt_epoch_MI,500,filt_epoch_MI.fs);
        subplot(4,4,i)
        plot(freq_1(1:55),10*log10(pwelch_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_MI_onechannel{i}(1:55)),freq_1(1:55),10*log10(filt_pwelch_bas_onechannel(1:55)),freq_2(1:55),10*log10(filt_pwelch_MI_onechannel(1:55)));
        xlabel('Frequency [Hz]');
        ylabel('Power Spectral Density [dB]');
        legend('Raw data baseline','Raw data MI','Filtered baseline','Filtered MI');
        legend('Location','best');
        title(sprintf('Comparison pwelch butter VS pwelch raw for channels %s',epoch_baseline.channels{1,i}));
    end
end
%% spatial filtering on the raw data.CAR
if Action(1)==3
    medium_channels=mean(s');
    signal_car=zeros(size(s,1),size(s,2));
    for i=1:size(s,1)
        signal_car(i,:)=s(i,:)-medium_channels(1,i);
        
    end
    
    % plot(s(:,9))
    % hold on
    % plot(signal_car(:,9)')
    % title (sprintf('car filter and raw signal for the channel %d',9));
    % xlabel('CAR signal');
    % ylabel('raw signal');
    
    session_filt_CAR=session;
    session_filt_CAR.data=signal_car';
    
end
%% p_welch on car data
if Action(2)==1
    filt_epoch_baseline_CAR=epoch_struct(session_filt_CAR,200,TimeBeforeEventBaseline,TimeAfterEventBaseline);
    filt_epoch_MI_CAR=epoch_struct(session_filt_CAR,400,TimeBeforeEventMI,TimeAfterEventMI);
    filt_epoch_MI_CAR_termination=epoch_struct(session_filt_CAR,555,TimeBeforeEventMItermination,TimeAfterEventMItermination);
    
    for i=1:size(chanlocs16,2)
        [pwelch_car_bas_onechannel{i},freq_1]=pwelch_for_each_channel(i,filt_epoch_baseline_CAR,500,session_filt_CAR.fs);
        [pwelch_car_MI_onechannel{i},freq_2]=pwelch_for_each_channel(i,filt_epoch_MI_CAR,500,session_filt_CAR.fs);
        subplot(4,4,i)
        plot(freq_1(1:55),10*log10(pwelch_car_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_car_MI_onechannel{i}(1:55)))
        xlabel('Frequency [Hz]');
        ylabel('PSD [dB]');
        title(sprintf('SPD comparison for channels %s',filt_epoch_baseline_CAR.channels{1,i}));
        legend('baseline','MI');
    end
end
%% plot_comparison CAR_pwelch VS RAW_pwelch
if Action(3)==1;
for i=1:size(chanlocs16,2)
   [pwelch_car_bas_onechannel{i},freq_1]=pwelch_for_each_channel(i,filt_epoch_baseline_CAR,500,session_filt_CAR.fs); 
   [pwelch_car_MI_onechannel{i},freq_2]=pwelch_for_each_channel(i,filt_epoch_MI_CAR,500,session_filt_CAR.fs); 
   
   [pwelch_bas_onechannel{i},freq_1]=pwelch_for_each_channel(i,epoch_baseline,500,epoch_baseline.fs); 
   [pwelch_MI_onechannel{i},freq_2]=pwelch_for_each_channel(i,epoch_MI,500,epoch_MI.fs); 
   
   subplot(4,4,i)
   plot(freq_1(1:55),10*log10(pwelch_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_MI_onechannel{i}(1:55)))
   hold on
   plot(freq_1(1:55),10*log10(pwelch_car_bas_onechannel{i}(1:55)),freq_2(1:55),10*log10(pwelch_car_MI_onechannel{i}(1:55)))
   
   xlabel('Frequency [Hz]');
   ylabel('PSD [dB]');
   title(sprintf('SPD comparison for channels %s',filt_epoch_baseline_CAR.channels{1,i}));
   lgd=legend('bas raw','MI raw','bas car','MI car');
   lgd.FontSize=5.5;
end  
end

%% laplacian filter
if Action(1)==4
load('laplacian_16_10-20_mi.mat');
signal_laplacian = s(:,1:16)*lap; 

% subplot(3,1,1)
% plot(s(:,9))
% subplot(3,1,2)
% plot(signal_laplacian(:,9))
% subplot(3,1,3)
% plot(signal_car(:,9))

session_filt_lap=session;
session_filt_lap.data=signal_laplacian';

filt_epoch_baseline_lap=epoch_struct(session_filt_lap,200,TimeBeforeEventBaseline,TimeAfterEventBaseline);
filt_epoch_MI_lap=epoch_struct(session_filt_lap,400,TimeBeforeEventMI,TimeAfterEventMI);
filt_epoch_MI_lap_termination=epoch_struct(session_filt_lap,555,TimeBeforeEventMItermination,TimeAfterEventMItermination);
end

%% Spectrogram MI initiation
if Action(4)==1
Cyclic_freq=[5:0.1:40];

figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(epoch_baseline, epoch_MI, i, epoch_baseline.fs, epoch_baseline.fs-32, Cyclic_freq);
   subplot(4,4,i)
   t=t(1,:)+epoch_MI.time;
   imagesc('XData',t,'YData',f, 'CData',10*log(spect_for_one_channel)); 
   axis tight
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   h=colorbar;
   set(h,'ylim',[-5 5]);
   title((sprintf('Channel %s',epoch_baseline.channels{1,i})));
end
suptitle('Spectrogram raw data - MI initiation')
end

%% Spectrogram MI termination
if Action(4)==1
Cyclic_freq=[5:0.1:40];

figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(epoch_baseline, epoch_MI_termination, i, epoch_baseline.fs, epoch_baseline.fs-32, Cyclic_freq);
   subplot(4,4,i)
   t=t(1,:)+epoch_MI_termination.time;
   imagesc('XData',t,'YData',f, 'CData',10*log(spect_for_one_channel)); 
   axis tight
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   h=colorbar;
   set(h,'ylim',[-5 5]);
   title((sprintf('Channel %s',epoch_baseline.channels{1,i})));
end
suptitle('Spectrogram raw data - MI termination')
end

%% Spectrogram CAR MI initiation
if Action(4)==2
figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_CAR, filt_epoch_MI_CAR, i, filt_epoch_baseline_CAR.fs, filt_epoch_baseline_CAR.fs-32, Cyclic_freq);
   subplot(4,4,i)
    t=t(1,:)+ filt_epoch_MI_CAR.time;
   imagesc('XData',t,'YData',f,'CData', 10*log(spect_for_one_channel)); 
   axis tight
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   h=colorbar;
   set(h,'ylim',[-5 5]);
   title((sprintf('Channel %s',filt_epoch_baseline_CAR.channels{1,i})));
  
end

suptitle('Spectrogram CAR filtered data - MI initiation')
end
%% Spectrogram CAR MI termination
if Action(4)==2
figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_CAR, filt_epoch_MI_CAR_termination, i, filt_epoch_baseline_CAR.fs, filt_epoch_baseline_CAR.fs-32, Cyclic_freq);
   subplot(4,4,i);
   t=t(1,:)+ filt_epoch_MI_CAR_termination.time;
   imagesc('XData',t,'YData',f,'CData', 10*log(spect_for_one_channel));
   axis tight
   h=colorbar;
   set(h,'ylim',[-5 5]);
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   title((sprintf('Channel %s',filt_epoch_baseline_CAR.channels{1,i})));
   
end
suptitle('Spectrogram CAR filtered data - MI termination')
end

%% Spectrogram Laplacian MI initiation
if Action(4)==3
figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_lap, filt_epoch_MI_lap, i, filt_epoch_baseline_lap.fs, filt_epoch_baseline_lap.fs-32, Cyclic_freq);
   subplot(4,4,i)
   t=t(1,:)+ filt_epoch_MI_lap.time;
   imagesc('XData',t,'YData',f,'CData', 10*log(spect_for_one_channel)); % in order to put in line the ferquencies and teh time
   axis tight
   h=colorbar;
   set(h,'ylim',[-5 5]);
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   title((sprintf('Channel %s',filt_epoch_baseline_lap.channels{1,i})));
end
suptitle('Spectrogram Laplacian filtered data - MI initiation')
end
%% Spectrogram Laplacian MI termination
if Action(4)==3
figure
for i=1:16
    
   [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_lap, filt_epoch_MI_lap_termination, i, filt_epoch_baseline_lap.fs, filt_epoch_baseline_lap.fs-32, Cyclic_freq);
   subplot(4,4,i)
    t=t(1,:)+ filt_epoch_MI_lap_termination.time;
   imagesc('XData',t,'YData',f,'CData', 10*log(spect_for_one_channel)); % in order to put in line the ferquencies and teh time
   axis tight
   h=colorbar;
   set(h,'ylim',[-5 5]);
   xlabel('time[s]');
   ylabel('frequency[Hz]');
   title((sprintf('Channel %s',filt_epoch_baseline_lap.channels{1,i})));
end
suptitle('Spectrogram Laplacian filtered data - MI termination')
end
%% topoplot_raw data_fre:8-12;
if Action(5)==1
spect_for_one_channel_top=zeros(16,113); % 113 is the size of the time
Cyclic_freq2=[8:0.1:12];
for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(epoch_baseline, epoch_MI, i, epoch_baseline.fs, epoch_baseline.fs-32, Cyclic_freq2);
    for j=1:113
    
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end
end

    figure
    for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
    end
    suptitle('Topoplot raw data - alpha band');
end

%% topoplot raw data_fre:13-25;
if Action(5)==1
Cyclic_freq2=[13:0.1:25];
spect_for_one_channel_top=zeros(16,113);
for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(epoch_baseline, epoch_MI, i, epoch_baseline.fs, epoch_baseline.fs-32, Cyclic_freq2);
    for j=1:113
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end
end

    figure;
    for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
    end
    suptitle('Topoplot raw data - beta band');
end    
%% topoplot_car_fre:8-12;
if Action(5)==2
Cyclic_freq2=[8:0.1:12];
spect_for_one_channel_top=zeros(16,1);

for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_CAR, filt_epoch_MI_CAR, i, filt_epoch_baseline_CAR.fs, filt_epoch_baseline_CAR.fs-32, Cyclic_freq);
    for j=1:size(t,2)
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end 
end
    
    figure
    for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
    end
    suptitle('Topoplot CAR filtered data - alpha band');
end
%% topoplot_car_fre:13-30;
if Action(5)==2
Cyclic_freq2=[13:0.1:30];
spect_for_one_channel_top=zeros(16,1);

for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_CAR, filt_epoch_MI_CAR, i, filt_epoch_baseline_CAR.fs, filt_epoch_baseline_CAR.fs-32, Cyclic_freq);
    for j=1:size(t,2)
        
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end
end
    
   figure
   for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
   end
   suptitle('Topoplot CAR filtered data - beta band'); 
end

%% topoplot laplacian_fre:8-12;
if Action(5)==3
spect_for_one_channel_top=zeros(16,1);
Cyclic_freq2=[8:0.1:12];

for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_lap, filt_epoch_MI_lap, i,filt_epoch_baseline_lap.fs, filt_epoch_baseline_lap.fs-32, Cyclic_freq);
    
    for j=1:size(t,2)
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end
end
       
    figure
    for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
    end
    suptitle('Topoplot Laplacian filtered data - alpha band');
end

%% topoplot laplacian_fre:13-25
if Action(5)==3
spect_for_one_channel_top=zeros(16,1);
Cyclic_freq2=[13:0.1:25];

for i=1:16
  
    [spect_for_one_channel,t, f]=Spectrogram_function(filt_epoch_baseline_lap, filt_epoch_MI_lap, i,filt_epoch_baseline_lap.fs, filt_epoch_baseline_lap.fs-32, Cyclic_freq);
    for j=1:size(t,2)
         spect_for_one_channel_top(i,j)=mean(10*log(spect_for_one_channel(:,j)));
         
    end
end
  
    figure
    for j=1:9
        subplot(3,3,j)
        topoplot(spect_for_one_channel_top(:,1+13*(j-1)),chanlocs16,'style','both','electrodes','ptslabels','chaninfo', session.channels);
        title(sprintf('time : %.2f ',t(1+13*(j-1))));
    end
    suptitle('Topoplot Laplacian filtered data - beta band');
end

%%
% 
