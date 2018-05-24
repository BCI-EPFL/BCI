clear;
close;
clc;

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
    runconc.event.position=cat(1,runconc.event.position,sessionPSD{i}.event.position+size(runconc.event.position,3));
    
end

%create epoching 
% epoch_baseline=epoch_window(runconc,200,0,2,params_spectrogram.mlength,params_spectrogram.wshift);
% epoch_MI=epoch_window(runconc,400,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
% epoch_MI_term=epoch_window(runconc,555,0,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_baseline=epoch_window(runconc,200,2,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,400,3,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,3,params_spectrogram.mlength,params_spectrogram.wshift);

%% cross validation for each trial to avoid the problem of time-chronological event 
 
thresholdCross=0.75; %25% of the data is test set for online test

%-----
%BASELINE vs MI (cat of 75% of the baseline and MI data, labels for
%training and test set) 

EpochTraining.BaseMI.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MI.samples(:,:,1:thresholdCross*size(epoch_MI.samples,3)));
EpochTraining.BaseMI.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MI.labels(1:thresholdCross*size(epoch_MI.labels,1)));

EpochTesting.BaseMI.data=cat(3, epoch_baseline.samples(:,:,((thresholdCross*size(epoch_baseline.samples,3)+1):end)), epoch_MI.samples(:,:,(thresholdCross*size(epoch_MI.samples,3)+1):end));
EpochTesting.BaseMI.labels=cat(1,epoch_baseline.labels((thresholdCross*size(epoch_baseline.labels,1)+1):end), epoch_MI.labels((thresholdCross*size(epoch_MI.labels,1)+1):end));

%% marco 



NoTrainingSamplesPerClassBas=thresholdCross*size(epoch_baseline.samples,3);
NoTrainingSamplesPerClassMI=thresholdCross*size(epoch_MI.samples,3);

TrialsPerFold=NoTrainingSamplesPerClassBas/(epoch_baseline.duration*10);  

for i=1:10
    contStartBas=(TrialsPerFold*(i-1)*epoch_baseline.duration+1); 
    contEndBas=TrialsPerFold*i*epoch_baseline.duration;
    
    contStartMI=(TrialsPerFold*(i-1)*epoch_MI.duration+1); 
    contEndMI=TrialsPerFold*i*epoch_MI.duration;
    
    Folds.BaseMI.data{i}=cat(3,EpochTraining.BaseMI.data(:,:,contStartBas:contEndBas),EpochTraining.BaseMI.data(:,:,(contStartMI+NoTrainingSamplesPerClassBas):(contEndMI+NoTrainingSamplesPerClassBas)));
    Folds.BaseMI.labels{i}=cat(1,EpochTraining.BaseMI.labels(contStartBas:contEndBas),EpochTraining.BaseMI.labels((contStartMI+NoTrainingSamplesPerClassBas):(contEndMI+NoTrainingSamplesPerClassBas)));
end

%% Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)


%in this way i'll have folds.BaseMI.data organized in a such a way that
%(16*19)


a=Folds; 
Folds.BaseMI=rmfield(Folds.BaseMI,'data');
for i=1:10
    for j=1:size(a.BaseMI.data{1,i},3)
         
        Folds.BaseMI.data{1,i}(j,:)=reshape(a.BaseMI.data{1,i}(:,:,j),[1,16*19]); % each channel for all the frequencies
    end    
end

%---
%% Proper CV, with rank feat:Marco plus all classifier, all features are considered 
for i=1:10
   %definition of the folds
   TrainingFolds.BaseMI.data=cat(1,Folds.BaseMI.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMI.labels=cat(1,Folds.BaseMI.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMI.data=Folds.BaseMI.data{i};
   ValidationFold.BaseMI.labels=Folds.BaseMI.labels{i};
   
   %normalization
   [TrainingFolds.BaseMI.data, mu.BaseMI, sigma.BaseMI]=zscore(TrainingFolds.BaseMI.data);% Mu is the mean and Sigma is the standard devition
   ValidationFold.BaseMI.data=(ValidationFold.BaseMI.data-mu.BaseMI)./sigma.BaseMI;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.BaseMI.noPCA(i,:), power_feat.BaseMI.noPCA(i,:)] = rankfeat(TrainingFolds.BaseMI.data, TrainingFolds.BaseMI.labels,  'fisher');
   
   Classifier={'linear', 'diaglinear','diagquadratic'}; %the quadratic can not be performed
   
   for j=1:numel(Classifier)
       for Nsel=1:20 %304
          classifier.BaseMI.noPCA=fitcdiscr(TrainingFolds.BaseMI.data(:,sort(ind.BaseMI.noPCA(i,1:Nsel))),TrainingFolds.BaseMI.labels,'discrimtype',Classifier{1,j});
          [yhat.BaseMI.noPCA,PosteriorProb.BaseMI.noPCA,~]=predict(classifier.BaseMI.noPCA,ValidationFold.BaseMI.data(:,sort(ind.BaseMI.noPCA(i,1:Nsel))));
          ClassError.BaseMI.noPCA{j}(i,Nsel)=classerror(ValidationFold.BaseMI.labels,yhat.BaseMI.noPCA);
       Nsel
       end
   end
end

%% Choice of the classifier and the number of features with rankfeat 
%going to use after with pca

for j=1:numel(Classifier)
    MeanClassError.BaseMI.noPCA{j}=mean(ClassError.BaseMI.noPCA{j});
    plot(1:20,MeanClassError.BaseMI.noPCA{j});
   
    title('Classifier and class error ');
    legend('Classifier: Linear','Classifier: DiagLinear','diagquadratic');
    xlabel('features');
    ylabel('class error');
    hold on 
end

%from the plot the linear shows higher perfomance, let's consider 150
%features good one with error 0.15

Nselected.noPCA=4;   %for Elisabetta

%% average outside CV, put in order all the features from rank feat for the
%final plot
averageFisher.BaseMI.noPCA=zeros(size(power_feat.BaseMI.noPCA,2),1);
for i=1:size(power_feat.BaseMI.noPCA,2) % each features that we wanna find
    for j=1:10 %for each fold
        averageFisher.BaseMI.noPCA(i)=averageFisher.BaseMI.noPCA(i)+power_feat.BaseMI.noPCA(j,ind.BaseMI.noPCA(j,:)==i);
    end
    averageFisher.BaseMI.noPCA(i)=averageFisher.BaseMI.noPCA(i)/10;
end



%% reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMI.noPCA=reshape(averageFisher.BaseMI.noPCA,[19 16]); % reshape in a way to have 16 channels for rows and 19 channles for columns

%plot of Fisher's scores
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMI.noPCA);
title('Fisher scores - baseline vs MI');
axis tight
xlabel('frequencies');
ylabel('channels');
h=colorbar;

% %% third case: PCA BEFORE RANK-FEAT
% 
% for i=1:10
%    %definition of the folds
%    TrainingFolds.BaseMI.data=cat(1,Folds.BaseMI.data{[1:(i-1),(i+1):10]});
%    TrainingFolds.BaseMI.labels=cat(1,Folds.BaseMI.labels{[1:(i-1),(i+1):10]});
%    ValidationFold.BaseMI.data=Folds.BaseMI.data{i};
%    ValidationFold.BaseMI.labels=Folds.BaseMI.labels{i};
%    
%    %normalization
%    [TrainingFolds.BaseMI.data, mu.BaseMI, sigma.BaseMI]=zscore(TrainingFolds.BaseMI.data);% Mu is the mean and Sigma is the standard devition
%    ValidationFold.BaseMI.data=(ValidationFold.BaseMI.data-mu.BaseMI)./sigma.BaseMI;
%    
%    %Fisher's score: we need to save it for every feature for every
%    %iteration (we will make the average outside the CV loop)
% 
%    
%    [coeff_train, score_train, variance_train] = pca(TrainingFolds.BaseMI.data(:,:));
%    %[coeff_val, score_val, variance_val] = pca(ValidationFold.BaseMI.data(:,:));
%    score_val=ValidationFold.BaseMI.data*coeff_train;
%    
% %    TrainingFolds.BaseMI.dataPCArank=TrainingFolds.BaseMI.data*(coeff_train*coeff_train.');
% %    ValidationFold.BaseMI.dataPCArank=ValidationFold.BaseMI.data*(coeff_val*coeff_val.');
%   [ind.BaseMI.PCA(i,:), power_feat.BaseMI.PCA(i,:)] = rankfeat(score_train, TrainingFolds.BaseMI.labels,'fisher');
% 
% 
%   Classifier={'linear', 'diaglinear','diagquadratic','quadratic'};% quadratic no possible
%    for j=1:numel(Classifier)
%        for Nsel=1:50 %until 304 
%           classifier.BaseMI.PCA=fitcdiscr(score_train(:,(1:Nsel)),TrainingFolds.BaseMI.labels,'discrimtype',Classifier{1,j});
%           [yhat.BaseMI.PCA,PosteriorProb.BaseMI.PCA,~]=predict(classifier.BaseMI.PCA,score_val(:,(1:Nsel)));
%           ClassError.BaseMI.PCA{j}(i,Nsel)=classerror(ValidationFold.BaseMI.labels,yhat.BaseMI.PCA);
%        Nsel
%        end
%    end
% end
% 
% %% mean of the class error to chose the type of classifier with PCA and the number of features--> 
% 
% for j=1:numel(Classifier)
%     MeanClassError.BaseMI.PCA{j}=mean(ClassError.BaseMI.PCA{j});
%     plot(1:50,MeanClassError.BaseMI.PCA{j});
%    
%     title('Classifier and class error ');
%     legend('Classifier: Linear','Classifier: DiagLinear','Classifier: DiagQuadratic','Classifier: Quadratic');
%     xlabel('features');
%     ylabel('class error');
%     hold on 
% end
% 
% Nselected.PCA=15; %linear is still the best one, error is half smaller (0.05)
% 
% 
% %% average outside CV, put in order all the features from rank feat for the
% %final plot
% averageFisher.BaseMI.PCA=zeros(size(power_feat.BaseMI.PCA,2),1);
% for i=1:size(power_feat.BaseMI.PCA,2) % each features that we wanna find
%     for j=1:10 %for each fold
%         averageFisher.BaseMI.PCA(i)=averageFisher.BaseMI.PCA(i)+power_feat.BaseMI.PCA(j,ind.BaseMI.PCA(j,:)==i);
%     end
%     averageFisher.BaseMI.PCA(i)=averageFisher.BaseMI.PCA(i)/10;
% end
% 
% 
% 
% %reshaping average Fisher's scores from line to matrix to plot them
% averageFisher.BaseMI.PCA=reshape(averageFisher.BaseMI.PCA,[19 16])'; % reshape in a way to have 16 channels for rows and 19 channles for columns
% 
% %plot of Fisher's scores
% figure
% imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMI.PCA);
% title('Fisher scores - baseline vs MI');
% axis tight
% xlabel('frequencies');
% ylabel('channels');
% h=colorbar;

%no difference in the final spectrogram


%% BASELINE vs MI Termination
 EpochTraining.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MI_term.samples(:,:,1:thresholdCross*size(epoch_MI_term.samples,3)));
 EpochTraining.BaseMITerm.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MI_term.labels(1:thresholdCross*size(epoch_MI_term.labels,1)));
 
EpochTesting.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,(thresholdCross*size(epoch_baseline.samples,3)+1):end), epoch_MI_term.samples(:,:,(thresholdCross*size(epoch_MI_term.samples,3)+1):end));
EpochTesting.BaseMITerm.labels=cat(1, epoch_baseline.labels((thresholdCross*size(epoch_baseline.samples,1)+1):end), epoch_MI_term.labels((thresholdCross*size(epoch_MI_term.samples,1)+1):end));

%
NoTrainingSamplesPerClassBas=thresholdCross*size(epoch_baseline.samples,3);
NoTrainingSamplesPerClassMI=thresholdCross*size(epoch_MI_term.samples,3);

TrialsPerFold=NoTrainingSamplesPerClassBas/(epoch_baseline.duration*10);  

for i=1:10
    contStartBas=(TrialsPerFold*(i-1)*epoch_baseline.duration+1); 
    contEndBas=TrialsPerFold*i*epoch_baseline.duration;
    
    contStartMI=(TrialsPerFold*(i-1)*epoch_MI_term.duration+1); 
    contEndMI=TrialsPerFold*i*epoch_MI_term.duration;
    
    Folds.BaseMITerm.data{i}=cat(3,EpochTraining.BaseMITerm.data(:,:,contStartBas:contEndBas),EpochTraining.BaseMITerm.data(:,:,(contStartMI+NoTrainingSamplesPerClassBas):(contEndMI+NoTrainingSamplesPerClassBas)));
    Folds.BaseMITerm.labels{i}=cat(1,EpochTraining.BaseMITerm.labels(contStartBas:contEndBas),EpochTraining.BaseMITerm.labels((contStartMI+NoTrainingSamplesPerClassBas):(contEndMI+NoTrainingSamplesPerClassBas)));
end


%Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)
a=Folds; 
Folds.BaseMITerm=rmfield(Folds.BaseMITerm,'data');
for i=1:10
    for j=1:size(a.BaseMITerm.data{1,i},3)
        Folds.BaseMITerm.data{1,i}(j,:)=reshape(a.BaseMITerm.data{1,i}(:,:,j),[1,16*19]); 
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
    Classifier={'linear', 'diaglinear','diagquadratic'}; %the quadratic can not be performed
   
   for j=1:numel(Classifier)
    for Nsel=1:20
       classifier.BaseMITerm=fitcdiscr(TrainingFolds.BaseMITerm.data(:,ind.BaseMITerm(i,1:Nsel)),TrainingFolds.BaseMITerm.labels,'discrimtype',Classifier{1,j});
       [yhat.BaseMITerm,PosteriorProb.BaseMITerm,~]=predict(classifier.BaseMITerm,ValidationFold.BaseMITerm.data(:,ind.BaseMITerm(i,1:Nsel)));
       ClassError.BaseMITerm{j}(i,Nsel)=classerror(ValidationFold.BaseMITerm.labels,yhat.BaseMITerm);
       Nsel
   end
    end
end

%error
for j=1:numel(Classifier)
    MeanClassError.BaseMITerm.noPCA{j}=mean(ClassError.BaseMITerm{j});
    plot(1:20,MeanClassError.BaseMITerm.noPCA{j});
   
    title('Classifier and class error ');
    legend('Classifier: Linear','Classifier: DiagLinear','diagquadratic');
    xlabel('features');
    ylabel('class error');
    hold on 
end


%average outside CV
averageFisher.BaseMITerm=zeros(size(power_feat.BaseMITerm,2),1);
for i=1:size(power_feat.BaseMITerm,2)
    for j=1:10
        averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)+power_feat.BaseMITerm(j,ind.BaseMITerm(j,:)==i);
    end
    averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)/10;
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMITerm=reshape(averageFisher.BaseMITerm,[19 16]);

%plot of Fisher's scores.
figure
imagesc([4 40],[1 16],averageFisher.BaseMITerm)
title('Fisher scores - baseline vs MI termination');

%% MI initiation vs MI Termination

EpochTraining.BaseMITerm.data=cat(3, epoch_MI.samples(:,:,1:thresholdCross*size(epoch_MI.samples,3)), epoch_MI_term.samples(:,:,1:thresholdCross*size(epoch_MI_term.samples,3)));
EpochTraining.BaseMITerm.labels=cat(1,epoch_MI.labels(1:thresholdCross*size(epoch_MI.labels,1)), epoch_MI_term.labels(1:thresholdCross*size(epoch_MI_term.labels,1)));
 
EpochTesting.BaseMITerm.data=cat(3, epoch_MI.samples(:,:,(thresholdCross*size(epoch_MI.samples,3)+1):end), epoch_MI_term.samples(:,:,(thresholdCross*size(epoch_MI_term.samples,3)+1):end));
EpochTesting.BaseMITerm.labels=cat(1, epoch_MI.labels((thresholdCross*size(epoch_MI.labels,1)+1):end), epoch_MI_term.labels((thresholdCross*size(epoch_MI_term.labels,1)+1):end));

NoTrainingSamplesPerClass=thresholdCross*size(epoch_MI.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(epoch_MI.duration*10);

% folds separation
Folds.Base.data=[];
Folds.Base.labels=[];
k=0;
cont=1;
trials=0;
for j=1:10
for i=1:TrialsPerFold
    Folds.Base.data=cat(3,Folds.Base.data,EpochTraining.BaseMITerm.data(:,:,(cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMITerm.data(:,:,(cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    % -1 is to not consider twice the same number.
    Folds.Base.labels=cat(1,Folds.Base.labels,EpochTraining.BaseMITerm.labels((cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMITerm.labels((cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    k=1;
    cont=i*epoch_baseline.duration;
end

trials=cont;
cont=1;
k=0;
Folds.BaseMITerm.data{j}= Folds.Base.data;
Folds.BaseMITerm.labels{j}=Folds.Base.labels;
Folds.Base.data=[];
Folds.Base.labels=[];
end
 

%% Now we need to transform every sample from a matrix to a line (because
%rankfeat wants lines)
a=Folds; 
Folds.BaseMITerm=rmfield(Folds.BaseMITerm,'data');
for i=1:10
    for j=1:size(a.BaseMITerm.data{1,i},3)
        Folds.BaseMITerm.data{1,i}(j,:)=reshape(a.BaseMITerm.data{1,i}(:,:,j),[1,16*19]); 
    end    
end

%---
%% Proper CV
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

    Classifier={'linear', 'diaglinear','diagquadratic'}; %
   %loop over the number of features
   for j=1:numel(Classifier)
    for Nsel=1:20
       classifier.BaseMITerm=fitcdiscr(TrainingFolds.BaseMITerm.data(:,ind.BaseMITerm(1:Nsel)),TrainingFolds.BaseMITerm.labels,'discrimtype', Classifier{1,j});
       [yhat.BaseMITerm,PosteriorProb.BaseMITerm,~]=predict(classifier.BaseMITerm,ValidationFold.BaseMITerm.data(:,ind.BaseMITerm(1:Nsel)));
       ClassError.BaseMITerm{j}(i,Nsel)=classerror(ValidationFold.BaseMITerm.labels,yhat.BaseMITerm);
       Nsel
   end
    end
end

%error
for j=1:numel(Classifier)
    MeanClassError.BaseMITerm.noPCA{j}=mean(ClassError.BaseMITerm{j});
    plot(1:20,MeanClassError.BaseMITerm.noPCA{j});
   
    title('Classifier and class error ');
    legend('Classifier: Linear','Classifier: DiagLinear','diagquadratic');
    xlabel('features');
    ylabel('class error');
    hold on 
end

%% average outside CV
averageFisher.BaseMITerm=zeros(size(power_feat.BaseMITerm,2),1);
for i=1:size(power_feat.BaseMITerm,2)
    for j=1:10
        averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)+power_feat.BaseMITerm(j,ind.BaseMITerm(j,:)==i);
    end
    averageFisher.BaseMITerm(i)=averageFisher.BaseMITerm(i)/10;
end

%reshaping average Fisher's scores from line to matrix to plot them
averageFisher.BaseMITerm=reshape(averageFisher.BaseMITerm,[19 16]);

%plot of Fisher's scores.
figure
imagesc('XData',[4 40],'YData',[1 16],'CData',averageFisher.BaseMITerm)
axis tight
title('Fisher scores - baseline vs MI termination');
ylabel('channels');
h=colorbar;


%% TERMINATION

 
 
thresholdCross=0.75; %25% of the data is test set for online test

%-----
%BASELINE vs MI (cat of 75% of the baseline and MI data, labels for
%training and test set) 

EpochTraining.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,1:thresholdCross*size(epoch_baseline.samples,3)), epoch_MITerm.samples(:,:,1:thresholdCross*size(epoch_MITerm.samples,3)));
EpochTraining.BaseMITerm.labels=cat(1,epoch_baseline.labels(1:thresholdCross*size(epoch_baseline.labels,1)), epoch_MITerm.labels(1:thresholdCross*size(epoch_MITerm.labels,1)));
 
EpochTesting.BaseMITerm.data=cat(3, epoch_baseline.samples(:,:,(thresholdCross*size(epoch_baseline.samples,3)+1):end), epoch_MITerm.samples(:,:,(thresholdCross*size(epoch_MITerm.samples,3)+1):end));
EpochTesting.BaseMITerm.labels=cat(1, epoch_baseline.labels((thresholdCross*size(epoch_baseline.samples,1)+1):end), epoch_MITerm.labels((thresholdCross*size(epoch_MITerm.samples,1)+1):end));


%% QUESTA è LA SEZIONE A CUI MI RIFERIVO:
%%put each 33 sample for baseline and MI events in a way to be the following--> .

NoTrainingSamplesPerClass=thresholdCross*size(epoch_baseline.samples,3);
TrialsPerFold=NoTrainingSamplesPerClass/(epoch_baseline.duration*10);

Folds.Base.data=[];
Folds.Base.labels=[];
k=0;
cont=1;
trials=0;
for j=1:10
for i=1:TrialsPerFold
    Folds.Base.data=cat(3,Folds.Base.data,EpochTraining.BaseMITerm.data(:,:,(cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMI.data(:,:,(cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    % -1 is to not consider twice the same number.
    Folds.Base.labels=cat(1,Folds.Base.labels,EpochTraining.BaseMITerm.labels((cont+trials+k:epoch_baseline.duration*i+trials)),EpochTraining.BaseMI.labels((cont+trials+k+NoTrainingSamplesPerClass:(epoch_baseline.duration*i+NoTrainingSamplesPerClass+trials))));
    k=1;
    cont=i*epoch_baseline.duration;
end

trials=cont;
cont=1;
k=0;
Folds.BaseMiTerm.data{j}= Folds.Base.data;
Folds.BaseMiTerm.labels{j}=Folds.Base.labels;
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
Folds.BaseMiTerm=rmfield(Folds.BaseMiTerm,'data');
for i=1:10
    for j=1:size(a.BaseMiTerm.data{1,i},3)
         
        Folds.BaseMiTerm.data{1,i}(j,:)=reshape(a.BaseMiTerm.data{1,i}(:,:,j)',[1,16*19]); 
    end    
end

%---
%% Proper CV, with rank feat:Marco plus all classifier, all features are considered 
for i=1:10
   %definition of the folds
   TrainingFolds.BaseMiTerm.data=cat(1,Folds.BaseMiTerm.data{[1:(i-1),(i+1):10]});
   TrainingFolds.BaseMiTerm.labels=cat(1,Folds.BaseMiTerm.labels{[1:(i-1),(i+1):10]});
   ValidationFold.BaseMiTerm.data=Folds.BaseMiTerm.data{i};
   ValidationFold.BaseMiTerm.labels=Folds.BaseMiTerm.labels{i};
   
   %normalization
   [TrainingFolds.BaseMiTerm.data, mu.BaseMi, sigma.BaseMi]=zscore(TrainingFolds.BaseMiTerm.data);% Mu is the mean and Sigma is the standard devition
   ValidationFold.BaseMiTerm.data=(ValidationFold.BaseMiTerm.data-mu.BaseMi)./sigma.BaseMi;
   
   %Fisher's score: we need to save it for every feature for every
   %iteration (we will make the average outside the CV loop)
   [ind.BaseMiTerm(i,:), power_feat.BaseMiTerm(i,:)] = rankfeat(TrainingFolds.BaseMiTerm.data, TrainingFolds.BaseMiTerm.labels,  'fisher');
   
   Classifier={'linear', 'diaglinear','diagquadratic'}; %the quadratic can not be performed
   
   for j=1:numel(Classifier)
       for Nsel=1:50 %size(TrainingFolds.BaseMi.data,2)
          classifier.BaseMiTerm=fitcdiscr(TrainingFolds.BaseMiTerm.data(:,ind.BaseMiTerm(1:Nsel)),TrainingFolds.BaseMiTerm.labels,'discrimtype',Classifier{1,j});
          [yhat.BaseMiTerm,PosteriorProb.BaseMiTerm,~]=predict(classifier.BaseMiTerm,ValidationFold.BaseMiTerm.data(:,ind.BaseMiTerm(1:Nsel)));
          ClassError.BaseMiTerm{j}(i,Nsel)=classerror(ValidationFold.BaseMiTerm.labels,yhat.BaseMiTerm);
          fprintf('iteration number =%d /10, features number :%d/304\n',i,Nsel);
       end
   end
end

%% Choice of the classifier and the number of features with rankfeat 
%going to use after with pca

for j=1:numel(Classifier)
    MeanClassError.BaseMiTerm{j}=mean(ClassError.BaseMiTerm{j});
    plot(1:50,MeanClassError.BaseMi{j});
   
    title('Classifier and class error ');
    legend('Classifier: Linear','Classifier: DiagLinear','Classifier: DiagQuadratic');
    xlabel('features');
    ylabel('class error');
    hold on 
end

%from the plot the linear shows higher perfomance, let's consider 40
%features good one with error 0.2

Nselected=40;   %for Elisabetta

%% average outside CV, put in order all the features from rank feat for the
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