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
s = preprocess_spectrogram(signal_car,params_spectrogram);
% return session with data in 3D (a new dimension for frequency!)


%% Epoching
% do the epoching you will need to adapt and create a new fuction maybe (keep the last one!)
% + concantenate runs 
%epochs_PSD_ONSET  = epoching_function();

runconc.data=[];
runconc.event.name=[];
runconc.event.position=[];
runconc.freq=[];
for i=1:numel(s)
    runconc.data=cat(3,runconc.data,s{i}.data);% concatenation of all the 4 runs_data
    runconc.event.name=cat(1,runconc.event.name,s{i}.event.name);
    runconc.event.position=cat(1,runconc.event.position,s{i}.event.position+size(runconc.event.position,3));
    runconc.freq=s{1}.freq;
end

%create epoching 
epoch_baseline=epoch_window(runconc,200,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,400,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
 

%% cross validation for each trial to avoid the proble of time-chronological event 
 
thresholdCross=0.75; %25% of the data is test set

%-----
%BASELINE vs MI

EpochTraining.BaseMI.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MI.samples(:,:,1:thresholdCross*size(epoch_MI.samples,3)));
EpochTraining.BaseMI.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MI.labels(1:thresholdCross*size(epoch_MI.labels,1)));
 
EpochTesting.BaseMI.data=cat(3, epoch_baseline.samples(:,:,(thresholdCross*size(epoch_baseline.samples,3)+1):end), epoch_MI.samples(:,:,(thresholdCross*size(epoch_MI.samples,3)+1):end));
EpochTesting.BaseMI.labels=cat(1, epoch_baseline.labels((thresholdCross*size(epoch_baseline.samples,1)+1):end), epoch_MI.labels((thresholdCross*size(epoch_MI.samples,1)+1):end));

NoTrainingSamplesPerClass=thresholdCross*size(epoch_baseline.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(33*10);  %33=no. of windows per trial, 10=folds cv;
for i=1:10
    contStart=(TrialsPerFold*(i-1)*33+1);
    contEnd=TrialsPerFold*i*33;
    Folds.BaseMI.data{i}=cat(3,EpochTraining.BaseMI.data(:,:,contStart:contEnd),EpochTraining.BaseMI.data(:,:,(contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
    Folds.BaseMI.labels{i}=cat(1,EpochTraining.BaseMI.labels(contStart:contEnd),EpochTraining.BaseMI.labels((contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
end

%Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)
a=Folds; 
Folds.BaseMI=rmfield(Folds.BaseMI,'data');
for i=1:10
    for j=1:size(a.BaseMI.data{1,i},3)
        Folds.BaseMI.data{1,i}(j,:)=reshape(a.BaseMI.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end

%---
%Proper CV
for i=1:10
   %definition of the folds
   TrainingFolds.BaseMI.data=cat(1,Folds.BaseMI.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMI.labels=cat(1,Folds.BaseMI.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMI.data=Folds.BaseMI.data{i};
   ValidationFold.BaseMI.labels=Folds.BaseMI.labels{i};
   
   %normalization
   [TrainingFolds.BaseMI.data, mu.BaseMI, sigma.BaseMI]=zscore(TrainingFolds.BaseMI.data);
   ValidationFold.BaseMI.data=(ValidationFold.BaseMI.data-mu.BaseMI)./sigma.BaseMI;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.BaseMI(i,:), power_feat.BaseMI(i,:)] = rankfeat(TrainingFolds.BaseMI.data, TrainingFolds.BaseMI.labels,  'fisher');

   %loop over the number of features
   for Nsel=[1:14,15:5:50,75:25:304]
       classifier.BaseMI=fitcdiscr(TrainingFolds.BaseMI.data(:,ind.BaseMI(1:Nsel)),TrainingFolds.BaseMI.labels,'discrimtype','linear');
       [yhat.BaseMI,PosteriorProb.BaseMI,~]=predict(classifier.BaseMI,ValidationFold.BaseMI.data(:,ind.BaseMI(1:Nsel)));
       ClassError.BaseMI(i,Nsel)=classerror(ValidationFold.BaseMI.labels,yhat.BaseMI);
       Nsel
   end
end

MeanClassError.BaseMI=mean(ClassError.BaseMI);

%average outside CV
averageFisher.BaseMI=zeros(size(power_feat.BaseMI,2),1);
for i=1:size(power_feat.BaseMI,2)
    for j=1:10
        averageFisher.BaseMI(i)=averageFisher.BaseMI(i)+power_feat.BaseMI(j,ind.BaseMI(j,:)==i);
    end
    averageFisher.BaseMI(i)=averageFisher.BaseMI(i)/10;
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMI=reshape(averageFisher.BaseMI,[19 16])';

%plot of Fisher's scores
figure
imagesc([4 40],[1 16],averageFisher.BaseMI)
title('Fisher scores - baseline vs MI')


%-----
%BASELINE vs MI Termination

EpochTraining.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MI_term.samples(:,:,1:thresholdCross*size(epoch_MI_term.samples,3)));
EpochTraining.BaseMITerm.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MI_term.labels(1:thresholdCross*size(epoch_MI_term.labels,1)));
 
EpochTesting.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,(thresholdCross*size(epoch_baseline.samples,3)+1):end), epoch_MI_term.samples(:,:,(thresholdCross*size(epoch_MI_term.samples,3)+1):end));
EpochTesting.BaseMITerm.labels=cat(1, epoch_baseline.labels((thresholdCross*size(epoch_baseline.samples,1)+1):end), epoch_MI_term.labels((thresholdCross*size(epoch_MI_term.samples,1)+1):end));

NoTrainingSamplesPerClass=thresholdCross*size(epoch_baseline.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(33*10);  %33=no. of windows per trial, 10=folds cv;
for i=1:10
    contStart=(TrialsPerFold*(i-1)*33+1);
    contEnd=TrialsPerFold*i*33;
    Folds.BaseMITerm.data{i}=cat(3,EpochTraining.BaseMITerm.data(:,:,contStart:contEnd),EpochTraining.BaseMITerm.data(:,:,(contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
    Folds.BaseMITerm.labels{i}=cat(1,EpochTraining.BaseMITerm.labels(contStart:contEnd),EpochTraining.BaseMITerm.labels((contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
end

%Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)
a=Folds; 
Folds.BaseMITerm=rmfield(Folds.BaseMITerm,'data');
for i=1:10
    for j=1:size(a.BaseMITerm.data{1,i},3)
        Folds.BaseMITerm.data{1,i}(j,:)=reshape(a.BaseMITerm.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end

%---
%Proper CV
for i=1:10
   %definition of the folds
   TrainingFolds.BaseMITerm.data=cat(1,Folds.BaseMITerm.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMITerm.labels=cat(1,Folds.BaseMITerm.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMITerm.data=Folds.BaseMITerm.data{i};
   ValidationFold.BaseMITerm.labels=Folds.BaseMITerm.labels{i};
   
   %normalization
   [TrainingFolds.BaseMITerm.data, mu.BaseMITerm, sigma.BaseMITerm]=zscore(TrainingFolds.BaseMITerm.data);
   ValidationFold.BaseMITerm.data=(ValidationFold.BaseMITerm.data-mu.BaseMITerm)./sigma.BaseMITerm;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.BaseMITerm(i,:), power_feat.BaseMITerm(i,:)] = rankfeat(TrainingFolds.BaseMITerm.data, TrainingFolds.BaseMITerm.labels,  'fisher');

   %loop over the number of features
   for Nsel=[1:14,15:5:50,75:25:304]
       classifier.BaseMITerm=fitcdiscr(TrainingFolds.BaseMITerm.data(:,ind.BaseMITerm(1:Nsel)),TrainingFolds.BaseMITerm.labels,'discrimtype','linear');
       [yhat.BaseMITerm,PosteriorProb.BaseMITerm,~]=predict(classifier.BaseMITerm,ValidationFold.BaseMITerm.data(:,ind.BaseMITerm(1:Nsel)));
       ClassError.BaseMITerm(i,Nsel)=classerror(ValidationFold.BaseMITerm.labels,yhat.BaseMITerm);
       Nsel
   end
end

MeanClassError.BaseMITerm=mean(ClassError.BaseMITerm);

%average outside CV
averageFisher.BaseMITerm=zeros(size(power_feat.BaseMITerm,2),1);
for i=1:size(power_feat.BaseMITerm,2)
    for j=1:10
        averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)+power_feat.BaseMITerm(j,ind.BaseMITerm(j,:)==i);
    end
    averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)/10;
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMITerm=reshape(averageFisher.BaseMITerm,[19 16])';

%plot of Fisher's scores.
figure
imagesc([4 40],[1 16],averageFisher.BaseMITerm)
title('Fisher scores - baseline vs MI termination')
