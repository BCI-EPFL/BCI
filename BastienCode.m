clear;
close;
clc;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));

global SubjectID
global chanlocs16
global params_spectrogram

load('channel_location_16_10-20_mi');

x=input('Enter subject: 1 for Marco (ak4), 2 for Giammarco (ak5), 3 for Elisabetta (ak6) \n');


if x==1
    SubjectID="ak4";
    folderName='mi614';
    thresholdCross=2/3;
elseif x==2
    SubjectID="ak5";
    folderName='folder_runs_ak5_Giammarco';
    thresholdCross=0.75;
elseif x==3
    SubjectID="ak6";
    folderName='folder_runs_ak6';
    thresholdCross=0.75;
else 
    disp('Error in files assignment');
end

params_spectrogram.mlength    = 1;
params_spectrogram.wlength    = 0.5;
params_spectrogram.pshift     = 0.25;                  
params_spectrogram.wshift     = 0.0625;              
params_spectrogram.selchans   = 1:16;     
params_spectrogram.freq = 4:2:40;

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

%% CAR filtering
signal_car=[];

nRun = numel(s);
signal_car= s;

for iRun = 1:nRun 
    mean_channels = mean(s{iRun}.data(:,1:16),2);
    signal_car{iRun}.data(:,1:16) = signal_car{iRun}.data(:,1:16)- mean_channels*ones(1,16);
end

%% Extract PSD
sessionPSD = preprocess_spectrogram(signal_car,params_spectrogram);

%% Concatenation of the data and Epoching
% do the epoching you will need to adapt and create a new fuction maybe (keep the last one!)
% + concantenate runs 
%epochs_PSD_ONSET  = epoching_function();

runconc.data=[];
runconc.event.name=[];
runconc.event.position=[];
runconc.freq=sessionPSD{1}.freq;
for iFold=1:numel(sessionPSD)
    runconc.data=cat(3,runconc.data,sessionPSD{iFold}.data);% concatenation of all the 4 runs_data
    runconc.event.name=cat(1,runconc.event.name,sessionPSD{iFold}.event.name);
    runconc.event.position=cat(1,runconc.event.position,sessionPSD{iFold}.event.position+size(runconc.event.position,3));
    
end

epoch_baseline=epoch_window(runconc,200,0,2,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,555,-2,0,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,0.5,2.5,params_spectrogram.mlength,params_spectrogram.wshift);

%% MI initiation vs MI Termination

EpochTraining.MITerm.data=cat(3, epoch_MI.samples(:,:,1:floor(thresholdCross*size(epoch_MI.samples,3))), epoch_MI_term.samples(:,:,1:floor(thresholdCross*size(epoch_MI_term.samples,3))));
EpochTraining.MITerm.labels=cat(1,400*ones(floor(thresholdCross*size(epoch_MI.labels,1)),1), epoch_MI_term.labels(1:floor(thresholdCross*size(epoch_MI_term.labels,1))));
 
EpochTesting.MITerm.data=cat(3, epoch_MI.samples(:,:,(floor(thresholdCross*size(epoch_MI.samples,3))+1):end), epoch_MI_term.samples(:,:,(floor(thresholdCross*size(epoch_MI_term.samples,3))+1):end));
EpochTesting.MITerm.labels=cat(1, 400*ones(floor((1-thresholdCross)*size(epoch_MI.labels,1)),1), epoch_MI_term.labels((floor(thresholdCross*size(epoch_MI_term.labels,1))+1):end));

NoTrainingSamplesPerClass=thresholdCross*size(epoch_MI.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(epoch_MI.duration*10);

% folds separation
Folds.MI.data=[];
Folds.MI.labels=[];
cont=1;
for j=1:10
    for i=1:TrialsPerFold
        Folds.MI.data=cat(3,Folds.MI.data,EpochTraining.MITerm.data(:,:,(cont:(cont+(epoch_MI.duration-1)))),EpochTraining.MITerm.data(:,:,((cont+NoTrainingSamplesPerClass):(cont+(epoch_MI.duration-1)+NoTrainingSamplesPerClass))));
        % -1 is to not consider twice the same number.
        Folds.MI.labels=cat(1,Folds.MI.labels,EpochTraining.MITerm.labels(cont:(cont+(epoch_MI.duration-1))),EpochTraining.MITerm.labels((cont+NoTrainingSamplesPerClass):(cont+(epoch_MI.duration-1)+NoTrainingSamplesPerClass)));
        cont=cont+epoch_MI.duration;
    end
Folds.MITerm.data{j}= Folds.MI.data;
Folds.MITerm.labels{j}=Folds.MI.labels;
Folds.MI.data=[];
Folds.MI.labels=[];
end

%% Now we need to transform every sample from a matrix to a line 
a=Folds; 
Folds.MITerm=rmfield(Folds.MITerm,'data');
for i=1:10
    for j=1:size(a.MITerm.data{1,i},3)
        Folds.MITerm.data{1,i}(j,:)=reshape(a.MITerm.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end
 
%% Proper CV 
featuresList = 1:25;
nFold = 10;

[hyperparameters,minClassError,BestModel.ClassError,BestModel.AUC]=cross_validation(Folds, featuresList,nFold);

%% Training on the whole training 75% (66%) data set with the parameters selection

for j=1:size(EpochTraining.MITerm.data,3)
    dataTraining(j,:)=reshape(EpochTraining.MITerm.data(:,:,j)',[1,16*19]); 
end    
[dataTraining, Mu, Sigma]=zscore(dataTraining);
[indexTraining, powerTraining] = rankfeat(dataTraining, EpochTraining.MITerm.labels,  'fisher');

classifier=fitcdiscr(dataTraining(:,indexTraining(1:hyperparameters.Nsel)),EpochTraining.MITerm.labels,'discrimtype',  hyperparameters.ClassifierType);

%% Capability evaluation

%reshaping the 25%(33%) Test set
for j=1:size(EpochTesting.MITerm.data,3)
    dataTesting(j,:)=reshape(EpochTesting.MITerm.data(:,:,j)',[1,16*19]); 
end  

%normalization
dataTesting=(dataTesting-ones(size(dataTesting,1),1)*Mu)./(ones(size(dataTesting,1),1)*Sigma);

%predict
[yhatTesting,PosteriorProbTesting,~]=predict(classifier,dataTesting(:,indexTraining(1:hyperparameters.Nsel)));

%class error
Metrics.ClassError=classerror(EpochTesting.MITerm.labels,yhatTesting);
    
%ROC curve
[X,Y,T,Metrics.AUC] = perfcurve(EpochTesting.MITerm.labels,PosteriorProbTesting(:,2),555);
figure
plot(X,Y);
xlabel('False positive rate') 
ylabel('True positive rate')
title(sprintf('ROC curve -  %s', SubjectID))

%confusion matrix
Metrics.ConfMat=confusionmat(EpochTesting.MITerm.labels,yhatTesting);

figure
histogram(dataTesting(EpochTesting.MITerm.labels==400,indexTraining(1)));
hold on
histogram(dataTesting(EpochTesting.MITerm.labels==555,indexTraining(1)));

figure
histogram(dataTraining(EpochTesting.MITerm.labels==400,indexTraining(1)));
hold on
histogram(dataTraining(EpochTesting.MITerm.labels==555,indexTraining(1)));

figure
scatter(dataTesting(EpochTesting.MITerm.labels==400,indexTraining(1)),dataTesting(EpochTesting.MITerm.labels==400,indexTraining(2)));
hold on
scatter(dataTesting(EpochTesting.MITerm.labels==555,indexTraining(1)),dataTesting(EpochTesting.MITerm.labels==400,indexTraining(2)));
line

%% ONLINE
[Trials,SmoothedTotal]=real_online(s{1,length(s)},Mu,Sigma,classifier,indexTraining(1:hyperparameters.Nsel),0.96);