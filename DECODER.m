clear;
close;
clc;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));
addpath(genpath('codeProject1'));

load('channel_location_16_10-20_mi');

folderName =  'mi614';

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

%% CAR FILTER
signal_car=[];

for j=1:numel(s)
 medium_channels=mean(s{j}.data');
    signal_car{j}=s{j};
    for i=1:size(s{j}.data,1)
        signal_car{j}.data(i,:)=s{j}.data(i,:)-medium_channels(1,i);
        
    end
end
%% Extract PSD
sessionPSD = preprocess_spectrogram(signal_car,params_spectrogram);
% return session with data in 3D (a new dimension for frequency!)

%elisabetta: all the events starts 1 second before the real one (backwards)
%at the beginning of each window.

%% Concatenation of the data and Epoching
% do the epoching you will need to adapt and create a new fuction maybe (keep the last one!)
% + concantenate runs 
%epochs_PSD_ONSET  = epoching_function();

runconc.data=[];
runconc.event.name=[];
runconc.event.position=[];
runconc.freq=sessionPSD{1}.freq;
for i=1:numel(sessionPSD)
    runconc.data=cat(3,runconc.data,sessionPSD{i}.data);% concatenation of all the 4 runs_data
    runconc.event.name=cat(1,runconc.event.name,sessionPSD{i}.event.name);
    runconc.event.position=cat(1,runconc.event.position,sessionPSD{i}.event.position);%+size(runconc.event.position,3));
    
end

%create epoching 
% epoch_baseline=epoch_window(runconc,200,0,2,params_spectrogram.mlength,params_spectrogram.wshift);
% epoch_MI=epoch_window(runconc,400,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
% epoch_MI_term=epoch_window(runconc,555,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_baseline=epoch_window(runconc,200,0,2,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,555,-2,0,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,0.5,2.5,params_spectrogram.mlength,params_spectrogram.wshift);




%% MI initiation vs MI Termination

thresholdCross=0.75;

EpochTraining.MITerm.data=cat(3, epoch_MI_term.samples(:,:,1:floor(thresholdCross*size(epoch_MI_term.samples,3))), epoch_MI.samples(:,:,1:floor(thresholdCross*size(epoch_MI.samples,3))));
EpochTraining.MITerm.labels=cat(1,555*ones(floor(thresholdCross*size(epoch_MI_term.labels,1)),1), 400*ones((floor(thresholdCross*size(epoch_MI.labels,1))),1));
 
EpochTesting.MITerm.data=cat(3, epoch_MI_term.samples(:,:,(floor(thresholdCross*size(epoch_MI_term.samples,3))+1):end), epoch_MI.samples(:,:,(floor(thresholdCross*size(epoch_MI.samples,3))+1):end));
EpochTesting.MITerm.labels=cat(1, 555*ones(ceil((1-thresholdCross)*size(epoch_MI_term.labels,1)),1), 400*(ones((ceil((1-thresholdCross)*size(epoch_MI.labels,1))),1)));

NoTrainingSamplesPerClass=thresholdCross*size(epoch_MI.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(epoch_MI.duration*10);

% folds separation
Folds.MI.data=[];
Folds.MI.labels=[];
k=0;
cont=1;
trials=0;
for j=1:10
    for i=1:TrialsPerFold
        Folds.MI.data=cat(3,Folds.MI.data,EpochTraining.MITerm.data(:,:,((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.data(:,:,((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        % -1 is to not consider twice the same number.
        Folds.MI.labels=cat(1,Folds.MI.labels,EpochTraining.MITerm.labels(((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.labels(((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        k=1;
        cont=i*epoch_MI.duration;
    end

trials=cont;
cont=1;
k=0;
Folds.MITerm.data{j}= Folds.MI.data;
Folds.MITerm.labels{j}=Folds.MI.labels;
Folds.MI.data=[];
Folds.MI.labels=[];
end
 

%% Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)
a=Folds; 
Folds.MITerm=rmfield(Folds.MITerm,'data');
for i=1:10
    for j=1:size(a.MITerm.data{1,i},3)
        Folds.MITerm.data{1,i}(j,:)=reshape(a.MITerm.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end

%---
%% Proper CV
for i=1:10
   %definition of the folds
   TrainingFolds.MITerm.data=cat(1,Folds.MITerm.data{[1:(i-1),(i+1):10]});
   TrainingFolds.MITerm.labels=cat(1,Folds.MITerm.labels{[1:(i-1),(i+1):10]});
   ValidationFold.MITerm.data=Folds.MITerm.data{i};
   ValidationFold.MITerm.labels=Folds.MITerm.labels{i};
   
   %normalization
    [TrainingFolds.MITerm.data, mu.MITerm, sigma.MITerm]=zscore(TrainingFolds.MITerm.data);
    ValidationFold.MITerm.data=(ValidationFold.MITerm.data-mu.MITerm)./sigma.MITerm;
%    
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.MITerm(i,:), power_feat.MITerm(i,:)] = rankfeat(TrainingFolds.MITerm.data, TrainingFolds.MITerm.labels,  'fisher');

    Classifier={'linear', 'diaglinear','diagquadratic'}; %
   %loop over the number of features
   for j=1:numel(Classifier)
    for Nsel=1:20%304
       classifier.MITerm=fitcdiscr(TrainingFolds.MITerm.data(:,ind.MITerm(i,1:Nsel)),TrainingFolds.MITerm.labels,'discrimtype', Classifier{1,j});
       [yhat.MITerm,PosteriorProb.MITerm{i},~]=predict(classifier.MITerm,ValidationFold.MITerm.data(:,ind.MITerm(i,1:Nsel)));
       ClassError.MITerm{j}(i,Nsel)=classerror(ValidationFold.MITerm.labels,yhat.MITerm);
       fprintf('iteration number =%d /10 , classifier=%d, features number :%d/304\n',i,j,Nsel);
       
   end
    end
end

% %ROC CURVE
% 
% [X,Y] = perfcurve(ValidationFold.MITerm.labels,PosteriorProb.MITerm{10}(:,2),555);
% plot(X,Y);

%error--> selection Classifier and number of features
for j=1:numel(Classifier)
    MeanClassError.MITerm.noPCA{j}=mean(ClassError.MITerm{j});
    plot(1:20,MeanClassError.MITerm.noPCA{j});
   
    title('Classifier and class error ');
    legend('Classifier: Linear','Classifier: DiagLinear','diagquadratic');
    xlabel('features');
    ylabel('class error');
    hold on 
end

% %% average outside CV
% averageFisher.MITerm=zeros(size(power_feat.MITerm,2),1);
% for i=1:size(power_feat.MITerm,2)
%     for j=1:10
%         averageFisher.MITerm(i)=averageFisher.MITerm(i)+power_feat.MITerm(j,ind.MITerm(j,:)==i);
%     end
%     averageFisher.MITerm(i)=averageFisher.MITerm(i)/10;
% end
% 
% %reshaping average Fisher's scores from line to matrix to plot them
% averageFisher.MITerm=(reshape(averageFisher.MITerm,[19 16]))';
% 
% %plot of Fisher's scores.
% figure
% imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.MITerm)
% axis tight
% title('Fisher scores - MI initiation vs MI termination');
% xlabel('Frequency [Hz]');
% ylabel('channels');
% set(gca,'yTick',1:16,'YTickLabel', {chanlocs16.labels});
% h=colorbar;

%% Train on all data--> feature selection based on the number Nsel
EpochTraining.MITerm.dataCorrected=[];
EpochTraining.MITerm.labelsCorrected=[];
k=0;
cont=1;
trials=0;

    for i=1:TrialsPerFold*10
        EpochTraining.MITerm.dataCorrected=cat(3,EpochTraining.MITerm.dataCorrected,EpochTraining.MITerm.data(:,:,((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.data(:,:,((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        % -1 is to not consider twice the same number.
        EpochTraining.MITerm.labelsCorrected=cat(1,EpochTraining.MITerm.labelsCorrected,EpochTraining.MITerm.labels(((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.labels(((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        k=1;
        cont=i*epoch_MI.duration;
    end
    
    
for j=1:size(EpochTraining.MITerm.dataCorrected,3)
        EpochTraining.MITerm.dataNorm(j,:)=reshape(EpochTraining.MITerm.dataCorrected(:,:,j)',[1,16*19]); 
end



EpochTest.MITerm.dataCorrected=[];
EpochTest.MITerm.labelsCorrected=[];
k=0;
cont=1;
trials=0;

    for i=1:30
        EpochTest.MITerm.dataCorrected=cat(3,EpochTest.MITerm.dataCorrected,EpochTesting.MITerm.data(:,:,((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.data(:,:,((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        % -1 is to not consider twice the same number.
        EpochTest.MITerm.labelsCorrected=cat(1,EpochTest.MITerm.labelsCorrected,EpochTesting.MITerm.labels(((cont+trials+k):(epoch_MI.duration*i+trials))),EpochTraining.MITerm.labels(((cont+trials+k+NoTrainingSamplesPerClass):(epoch_MI.duration*i+NoTrainingSamplesPerClass+trials))));
        k=1;
        cont=i*epoch_MI.duration;
    end

for j=1:size(EpochTest.MITerm.dataCorrected,3)%--> to do better in order like training
        EpochTesting.MITerm.dataNorm(j,:)=reshape(EpochTest.MITerm.dataCorrected(:,:,j)',[1,16*19]); 
end 


[EpochTraining.MITerm.dataNorm, mu.MINorm, sigma.MINorm]=zscore(EpochTraining.MITerm.dataNorm);% Mu is the mean and Sigma is the standard devition
EpochTesting.MITerm.dataNorm=(EpochTesting.MITerm.dataNorm-mu.MINorm)./sigma.MINorm;
   
[ind.MITermNorm, power_feat.MITermNorm] = rankfeat(EpochTraining.MITerm.dataNorm, EpochTraining.MITerm.labelsCorrected,  'fisher');



for i=1:size(power_feat.MITermNorm,2)
averageonFisher(1,i)=power_feat.MITermNorm(ind.MITermNorm==i);
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisherReshape=(reshape(averageonFisher,[19 16]))';

%plot of Fisher's scores.
figure;
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisherReshape)
axis tight
title('Fisher scores - MI initiation vs MI termination');
xlabel('Frequency [Hz]');
ylabel('channels');
set(gca,'yTick',1:16,'YTickLabel', {chanlocs16.labels});
h=colorbar;


%% ROC CURVE


Nsel=20;

classifier.MITermNorm=fitcdiscr(EpochTraining.MITerm.dataNorm(:,ind.MITermNorm(1,1:Nsel)),EpochTraining.MITerm.labelsCorrected,'discrimtype', 'linear');
[yhat.MITermNorm,PosteriorProb.MITermNorm,~]=predict(classifier.MITermNorm,EpochTesting.MITerm.dataNorm(:,ind.MITermNorm(1,1:Nsel)));

[X,Y] = perfcurve(EpochTest.MITerm.labelsCorrected,PosteriorProb.MITermNorm(:,2),555);
plot(X,Y);
% random curve is the bisettrice, 
