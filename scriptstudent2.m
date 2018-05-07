clear;
close;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));
addpath(genpath('codeProject1'));%--?

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

data_filter=[];
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
signal_car=[];

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

%elisabetta: all the events starts 1 second before the real one (backwards)
%at the beginning of each window.

%% Concatenation of the data and Epoching
% do the epoching you will need to adapt and create a new fuction maybe (keep the last one!)
% + concantenate runs 
%epochs_PSD_ONSET  = epoching_function();

runconc.data=[];
runconc.event.name=[];
runconc.event.position=[];
runconc.freq=s{1}.freq;
for i=1:numel(s)
    runconc.data=cat(3,runconc.data,s{i}.data);% concatenation of all the 4 runs_data
    runconc.event.name=cat(1,runconc.event.name,s{i}.event.name);
    runconc.event.position=cat(1,runconc.event.position,s{i}.event.position+size(runconc.event.position,3));
    
end

%create epoching 
epoch_baseline=epoch_window(runconc,200,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,400,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
 %% cross validation for each trial to avoid the problem of time-chronological event 
 
thresholdCross=0.75; %25% of the data is test set for online test

%-----
%BASELINE vs MI (cat of 75% of the baseline and MI data, labels for
%training and test set) 

EpochTraining.BaseMI.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MI.samples(:,:,1:thresholdCross*size(epoch_MI.samples,3)));
EpochTraining.BaseMI.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MI.labels(1:thresholdCross*size(epoch_MI.labels,1)));
 
EpochTesting.BaseMI.data=cat(3, epoch_baseline.samples(:,:,(thresholdCross*size(epoch_baseline.samples,3)+1):end), epoch_MI.samples(:,:,(thresholdCross*size(epoch_MI.samples,3)+1):end));
EpochTesting.BaseMI.labels=cat(1, epoch_baseline.labels((thresholdCross*size(epoch_baseline.samples,1)+1):end), epoch_MI.labels((thresholdCross*size(epoch_MI.samples,1)+1):end));


%% QUESTA è LA SEZIONE A CUI MI RIFERIVO:
%%put each 33 sample for baseline and MI events in a way to be following--> .

NoTrainingSamplesPerClass=thresholdCross*size(epoch_baseline.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(epoch_baseline.duration*10);

Folds.Base.data=[];
Folds.Base.labels=[];
k=0;
cont=1;
trials=0;
for j=1:10
for i=1:TrialsPerFold
    Folds.Base.data=cat(3,Folds.Base.data,EpochTraining.BaseMI.data(:,:,(cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMI.data(:,:,(cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    % -1 is to not consider twice the same number.
    Folds.Base.labels=cat(1,Folds.Base.labels,EpochTraining.BaseMI.labels((cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMI.labels((cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    k=1;
    cont=i*epoch_baseline.duration;
end

trials=cont;
cont=1;
k=0;
Folds.BaseMi.data{j}= Folds.Base.data;
Folds.BaseMi.labels{j}=Folds.Base.labels;
Folds.Base.data=[];
Folds.Base.labels=[];
end
 

%% marco 
% NoTrainingSamplesPerClass=thresholdCross*size(epoch_baseline.samples,3);
% TrialsPerFold=NoTrainingSamplesPerClass/(epoch_baseline.duration*10);  %33=no. of windows per each trial, 10=folds cv;
% for i=1:10
%     contStart=(TrialsPerFold*(i-1)*33+1);
%     contEnd=TrialsPerFold*i*33;
%     Folds.BaseMI.data{i}=cat(3,EpochTraining.BaseMI.data(:,:,contStart:contEnd),EpochTraining.BaseMI.data(:,:,(contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
%     Folds.BaseMI.labels{i}=cat(1,EpochTraining.BaseMI.labels(contStart:contEnd),EpochTraining.BaseMI.labels((contStart+NoTrainingSamplesPerClass):(contEnd+NoTrainingSamplesPerClass)));
% end

%% Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)

a=Folds; 
Folds.BaseMi=rmfield(Folds.BaseMi,'data');
for i=1:10
    for j=1:size(a.BaseMi.data{1,i},3)
         
        Folds.BaseMi.data{1,i}(j,:)=reshape(a.BaseMi.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end

%---
%% Proper CV, with rank feat:Marco plus all classifier
for i=1:10
   %definition of the folds
   TrainingFolds.BaseMi.data=cat(1,Folds.BaseMi.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMi.labels=cat(1,Folds.BaseMi.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMi.data=Folds.BaseMi.data{i};
   ValidationFold.BaseMi.labels=Folds.BaseMi.labels{i};
   
   %normalization
   [TrainingFolds.BaseMi.data, mu.BaseMi, sigma.BaseMi]=zscore(TrainingFolds.BaseMi.data);% Mu is the mean and Sigma is the standard devition
   ValidationFold.BaseMi.data=(ValidationFold.BaseMi.data-mu.BaseMi)./sigma.BaseMi;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.BaseMi(i,:), power_feat.BaseMi(i,:)] = rankfeat(TrainingFolds.BaseMi.data, TrainingFolds.BaseMi.labels,  'fisher');
   
   Classifier={'linear', 'diaglinear','diagquadratic'}; %the quadratic can not be performed
   
   for j=1:numel(Classifier)
       for Nsel=1:size(TrainingFolds.BaseMi.data,2)
          classifier.BaseMi=fitcdiscr(TrainingFolds.BaseMi.data(:,ind.BaseMi(1:Nsel)),TrainingFolds.BaseMi.labels,'discrimtype',Classifier{1,j});
          [yhat.BaseMi,PosteriorProb.BaseMi,~]=predict(classifier.BaseMi,ValidationFold.BaseMi.data(:,ind.BaseMi(1:Nsel)));
          ClassError.BaseMi{j}(i,Nsel)=classerror(ValidationFold.BaseMi.labels,yhat.BaseMi);
       %Nsel
       end
   end
end

%find the classifier and the number of features with rankfeat that i'm
%going to use after with pca

for j=1:numel(Classifier)
MeanClassError.BaseMi{j}=mean(ClassError.BaseMi{j});
[ClassError.MinValue{j},ClassError.MinInd{j}]=min(MeanClassError.BaseMi{j});
end


plot(1:304,mean(ClassError.BaseMi{1}));
title('Classifier and class error ');
legend('Classifier:linear');
xlabel('features');
ylabel('class error');

Nselected=200;   %from the graph, error of 0.1

%average outside CV, put in order all the features from rank feat for the
%final plot
averageFisher.BaseMi=zeros(size(power_feat.BaseMi,2),1);
for i=1:size(power_feat.BaseMi,2) % each features that we wanna find
    for j=1:10 %for each fold
        averageFisher.BaseMi(i)=averageFisher.BaseMi(i)+power_feat.BaseMi(j,ind.BaseMi(j,:)==i);
    end
    averageFisher.BaseMi(i)=averageFisher.BaseMi(i)/10;
end



%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMi=reshape(averageFisher.BaseMi,[19 16])'; % reshape in a way to have 16 channels for rows and 19 channles for columns

%plot of Fisher's scores
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMi);
title('Fisher scores - baseline vs MI');
axis tight
xlabel('frequencies');
ylabel('channels');
h=colorbar;

%% PCA and rank-feat: Elisabetta


%let's apply pca on the number of features reduces:N selected

for i=1:10
   %definition of the folds
   TrainingFolds.BaseMi.data=cat(1,Folds.BaseMi.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMi.labels=cat(1,Folds.BaseMi.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMi.data=Folds.BaseMi.data{i};
   ValidationFold.BaseMi.labels=Folds.BaseMi.labels{i};
   
   %normalization
   [TrainingFolds.BaseMi.data, mu.BaseMi, sigma.BaseMi]=zscore(TrainingFolds.BaseMi.data);% Mu is the mean and Sigma is the standard devition
   ValidationFold.BaseMi.data=(ValidationFold.BaseMi.data-mu.BaseMi)./sigma.BaseMi;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   %[ind.BaseMi(i,:), power_feat.BaseMi(i,:)] = rankfeat(TrainingFolds.BaseMi.data, TrainingFolds.BaseMi.labels,  'fisher');
   [coeff_train, score_train, variance_train] = pca(TrainingFolds.BaseMi.data(:,ind.BaseMi(i,1:Nselected)));
   [coeff_val, score_val, variance_val] = pca(ValidationFold.BaseMi.data(:,ind.BaseMi(i,1:Nselected)));
   
   
   %cross validation for choosing the PCA number of features,
   %projection
   TrainingFolds.BaseMi.dataPCA=TrainingFolds.BaseMi.data(:,ind.BaseMi(1:200))*(coeff_train*coeff_train');
   ValidationFold.BaseMi.dataPCA=ValidationFold.BaseMi.data(:,ind.BaseMi(1:200))*(coeff_val*coeff_val');
   
   
   Classifier={'linear', 'diaglinear','diagquadratic'};% quadratic no possible
   for j=1:numel(Classifier)
       for Nsel=1:200 
          classifier.BaseMiPCA=fitcdiscr(TrainingFolds.BaseMi.dataPCA(:,(1:Nsel)),TrainingFolds.BaseMi.labels,'discrimtype',Classifier{1,j});
          [yhat.BaseMiPCA,PosteriorProb.BaseMiPCA,~]=predict(classifier.BaseMiPCA,ValidationFold.BaseMi.dataPCA(:,(1:Nsel)));
          ClassError.BaseMiPCA{j}(i,Nsel)=classerror(ValidationFold.BaseMi.labels,yhat.BaseMiPCA);
       %Nsel
       end
   end
end

% mean of the class error to chose the type of classifier with PCA and the number of features--> 

for j=1:numel(Classifier)
MeanClassError.BaseMiPCA{j}=mean(ClassError.BaseMiPCA{j});
[ClassError.MinValuePCA{j},ClassError.MinIndPCA{j}]=min(MeanClassError.BaseMiPCA{j});
end

plot(1:200,mean(ClassError.BaseMiPCA{1}));
title('Classifier and class error ');
legend('Classifier:linear');
xlabel('features');
ylabel('class error');

Nselected=185;

%--> THERE IS NO DIFFERENCE WITH PCA AND WITHOUT PCA----%

%average outside CV, put in order all the features from rank feat for the
%final plot
averageFisher.BaseMi=zeros(size(power_feat.BaseMi,2),1);
for i=1:size(power_feat.BaseMi,2) % each features that we wanna find
    for j=1:10 %for each fold
        averageFisher.BaseMi(i)=averageFisher.BaseMi(i)+power_feat.BaseMi(j,ind.BaseMi(j,:)==i);
    end
    averageFisher.BaseMi(i)=averageFisher.BaseMi(i)/10;
end



%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMi=reshape(averageFisher.BaseMi,[19 16])'; % reshape in a way to have 16 channels for rows and 19 channles for columns

%plot of Fisher's scores
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMi);
title('Fisher scores - baseline vs MI');
axis tight
xlabel('frequencies');
ylabel('channels');
h=colorbar;


%% third case: PCA BEFORE RANK-FEAT

for i=1:10
   %definition of the folds
   TrainingFolds.BaseMi.data=cat(1,Folds.BaseMi.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMi.labels=cat(1,Folds.BaseMi.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMi.data=Folds.BaseMi.data{i};
   ValidationFold.BaseMi.labels=Folds.BaseMi.labels{i};
   
   %normalization
   [TrainingFolds.BaseMi.data, mu.BaseMi, sigma.BaseMi]=zscore(TrainingFolds.BaseMi.data);% Mu is the mean and Sigma is the standard devition
   ValidationFold.BaseMi.data=(ValidationFold.BaseMi.data-mu.BaseMi)./sigma.BaseMi;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)

   
   [coeff_train, score_train, variance_train] = pca(TrainingFolds.BaseMi.data(:,:));
   [coeff_val, score_val, variance_val] = pca(ValidationFold.BaseMi.data(:,:));
   TrainingFolds.BaseMi.dataPCArank=TrainingFolds.BaseMi.data*(coeff_train*coeff_train.');
   ValidationFold.BaseMi.dataPCArank=ValidationFold.BaseMi.data*(coeff_val*coeff_val.');
   [ind.BaseMiPCArank(i,:), power_feat.BaseMiPCArank(i,:)] = rankfeat(TrainingFolds.BaseMi.dataPCArank, TrainingFolds.BaseMi.labels,'fisher');


  Classifier={'linear', 'diaglinear','diagquadratic'};% quadratic no possible
   for j=1:numel(Classifier)
       for Nsel=1:304 
          classifier.BaseMiPCArank=fitcdiscr(TrainingFolds.BaseMi.dataPCArank(:,(1:Nsel)),TrainingFolds.BaseMi.labels,'discrimtype',Classifier{1,j});
          [yhat.BaseMiPCArank,PosteriorProb.BaseMiPCArank,~]=predict(classifier.BaseMiPCArank,ValidationFold.BaseMi.dataPCArank(:,(1:Nsel)));
          ClassError.BaseMiPCA{j}(i,Nsel)=classerror(ValidationFold.BaseMi.labels,yhat.BaseMiPCArank);
       %Nsel
       end
   end
end

% mean of the class error to chose the type of classifier with PCA and the number of features--> 

for j=1:numel(Classifier)
MeanClassError.BaseMiPCA{j}=mean(ClassError.BaseMiPCA{j});
[ClassError.MinValuePCA{j},ClassError.MinIndPCA{j}]=min(MeanClassError.BaseMiPCA{j});
end

plot(1:304,mean(ClassError.BaseMiPCA{1}));% min error 0.0396, linear
title('Classifier and class error ');
legend('Classifier:linear');
xlabel('features');
ylabel('class error');

Nselected=150;

%--> THERE IS NO DIFFERENCE WITH PCA AND WITHOUT PCA----%

%average outside CV, put in order all the features from rank feat for the
%final plot
averageFisher.BaseMiPCArank=zeros(size(power_feat.BaseMiPCArank,2),1);
for i=1:size(power_feat.BaseMiPCArank,2) % each features that we wanna find
    for j=1:10 %for each fold
        averageFisher.BaseMiPCArank(i)=averageFisher.BaseMiPCArank(i)+power_feat.BaseMiPCArank(j,ind.BaseMiPCArank(j,:)==i);
    end
    averageFisher.BaseMiPCArank(i)=averageFisher.BaseMiPCArank(i)/10;
end



%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMiPCArank=reshape(averageFisher.BaseMiPCArank,[19 16])'; % reshape in a way to have 16 channels for rows and 19 channles for columns

%plot of Fisher's scores
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMiPCArank);
title('Fisher scores - baseline vs MI');
axis tight
xlabel('frequencies');
ylabel('channels');
h=colorbar;


%% -----
%%BASELINE vs MI Termination

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
title('Fisher scores - baseline vs MI termination');