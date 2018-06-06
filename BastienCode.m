clear;
close;
clc;

addpath(genpath('biosig'));
addpath(genpath('folder_runs'));
addpath(genpath('data'));
addpath(genpath('eeglab13_4_4b'));


load('channel_location_16_10-20_mi');

folderName =  'folder_runs_ak5_Giammarco';

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

%% CAR FILTER
signal_car=[];

% for iClassifier=1:numel(s)
%  medium_channels=mean(s{iClassifier}.data');
%     signal_car{iClassifier}=s{iClassifier};
%     for iFold=1:size(s{iClassifier}.data,1)
%         signal_car{iClassifier}.data(iFold,:)=s{iClassifier}.data(iFold,:)-medium_channels(1,iFold);
%         
%     end
% end

nRun = numel(s);
signal_car= s;

for iRun = 1:nRun 
    mean_channels = mean(s{iRun}.data(:,1:16),2);
    signal_car{iRun}.data(:,1:16) = signal_car{iRun}.data(:,1:16)- mean_channels*ones(1,16);
end

%% Extract PSD

sessionPSD = preprocess_spectrogram(signal_car,params_spectrogram);
% Aux=signal_car;
% signal_car=[];
% for iRun=1:nRun
%     for i=1:length(Aux)
%     j = 1;
%         for t = 512:32:size(Aux{1,iRun}.data,1) %we remove the last window (not of 1 second)
%             for idxChannels = 1:16
%                 signal_car{iRun}.Windows(:,idxChannels,j)=Aux{1,iRun}.data(t-512+1:t,idxChannels);
%             end 
%             j= j+1;
%             j
%         end
%     end
% end
% 
% sessionPSD=[];
% for iRun=2:nRun %1:nRun
%     sessionPSD{iRun}.data=zeros(16,19,size(signal_car{iRun}.Windows,3));
%     for j=1:size(signal_car{iRun}.Windows,3)
%         for idxChannels=1:16
%             [psd,freqgrid] = pwelch(signal_car{iRun}.Windows(:,idxChannels,j),0.5*512,0.4375*512,f,512);
%             psd=log(psd);
% %             [freqs, idfreqs] = intersect(freqgrid,f);
% %             psd = psd(idfreqs);
%             sessionPSD{iRun}.data(idxChannels,:,j)=psd;
%         end
%     j    
%     end
% end
sessionPSD = preprocess_spectrogram(signal_car,params_spectrogram);
%sessionPSD2 = preprocess_spectrogram(Aux,params_spectrogram);
% for i=1:4
% sessionPSD{1,i}.rate=sessionPSD2{1,i}.rate;
% sessionPSD{1,i}.event=sessionPSD2{1,i}.event;
% sessionPSD{1,i}.freq=sessionPSD2{1,i}.freq;
% end
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
for iFold=1:numel(sessionPSD)
    runconc.data=cat(3,runconc.data,sessionPSD{iFold}.data);% concatenation of all the 4 runs_data
    runconc.event.name=cat(1,runconc.event.name,sessionPSD{iFold}.event.name);
    runconc.event.position=cat(1,runconc.event.position,sessionPSD{iFold}.event.position+size(runconc.event.position,3));
    
end

epoch_baseline=epoch_window(runconc,200,0,2,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI=epoch_window(runconc,555,-2,0,params_spectrogram.mlength,params_spectrogram.wshift);
epoch_MI_term=epoch_window(runconc,555,0.5,2.5,params_spectrogram.mlength,params_spectrogram.wshift);

%% MI initiation vs MI Termination

thresholdCross=0.75;

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

% hyperparameters
Classifier={'linear', 'diaglinear','diagquadratic','quadratic'};
featuresList = 1:20;
nFold = 10;

Pfeature = [];
class_error = [];
yPredicted = [];
PosteriorProb = [];
X=[];
Y=[];
Xtot=[];
Ytot=[];
AUCtot=[];

for iFold=1:nFold
    fprintf('Fold %d \n', iFold);
    
    %definition of the folds
    xTrain{iFold}=cat(1,Folds.MITerm.data{[1:(iFold-1),(iFold+1):10]});
    yTrain{iFold} =cat(1,Folds.MITerm.labels{[1:(iFold-1),(iFold+1):10]});
    xTest=Folds.MITerm.data{iFold};
    yTest{iFold} =Folds.MITerm.labels{iFold};
    
    %normalization
    [xTrain{iFold}, mu, sigma]=zscore(xTrain{iFold});
    xTest=(xTest-ones(size(xTest,1),1)*mu)./(ones(size(xTest,1),1)*sigma);
    
    %Fisher's score: we need to save it for every feature for every
    %iteration (we will make the average outside the CV loop)
    [indexPower, power_feat] = rankfeat(xTrain{iFold}, yTrain{iFold},  'fisher');
    [~,orderedInd] = sort(indexPower,'ascend'); %from 1 to 304
    Pfeature(iFold,:) = power_feat(orderedInd);
    
    %loop over the number of features
    for iClassifier =1:numel(Classifier)
        for iFeature= featuresList
            
            featuresSelected = xTrain{iFold}(:,indexPower(1:iFeature));
            model =fitcdiscr(featuresSelected,yTrain{iFold},'discrimtype', Classifier{iClassifier});
            
            [yPredicted{iFold},PosteriorProb{iFold},~]=predict(model,xTest(:,indexPower(1:iFeature)));
            class_error(iFold,iClassifier,iFeature)=classerror(yTest{iFold},yPredicted{iFold});
            [X{iFold},Y{iFold},T{iFold},AUC{iFold}] = perfcurve(yTest{iFold},PosteriorProb{iFold}(:,2),555);
           
            
        end
        XforRoc{iClassifier,iFold}=X{iFold};
        YforRoc{iClassifier,iFold}=Y{iFold};
        AUCforRoc{iClassifier,iFold}=AUC{iFold};
        Xtot=[];
        Ytot=[];
        AUCtot=[];
    end
end

%error
MinClassError=0.5;
figure
for iClassifier=1:numel(Classifier)
    class_error4Classifier = squeeze(class_error(:,iClassifier,:));
    
    meanClassError = mean(class_error4Classifier,1);
    [m,i]=min(meanClassError);
    if m<MinClassError
        MinClassError=m;
        hyperparameters.Nsel=i;
        hyperparameters.ClassifierType=Classifier{iClassifier};
    end
    % stdError(iClassifier,:)=std(meanClassError); %for each classifier
    stdError(iClassifier,:)=std(class_error4Classifier); %for each classifier
    %errorbar(featuresList,meanClassError,stdError(iClassifier,:));
    plot(featuresList,meanClassError)
    hold on
end
title('Classifier and class error ');
legend('Classifier: Linear','Classifier: DiagLinear','Classifier:Diagquadratic','Classifier:Quadratic');
xlabel('Features');
ylabel('Class error');


%% perfcurve
XRoc=[];
YRoc=[];

    for i=1:(iFold)
 XRoc=cat(2,XRoc,XforRoc{1,i});
 YRoc=cat(2,YRoc,YforRoc{1,i});   
    end
XRocC{1}=XRoc;
YRocC{1}=YRoc; XRoc=[];YRoc=[];

   for i=1:(iFold)
 XRoc=cat(2,XRoc,XforRoc{4,i});
 YRoc=cat(2,YRoc,YforRoc{4,i});   
    end
XRocC{4}=XRoc;
YRocC{4}=YRoc;

%it's not possible to perform the diaglinear or diagquadratic.
figure
plot(mean(XRocC{1},2),mean(YRocC{1},2));

hold on
plot(mean(XRocC{4},2),mean(YRocC{4},2));
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve ')
legend('Classifier: Linear','Classifier:Quadratic');




%% average outside CV

averageFisher = mean(Pfeature,1); 

%reshaping average Fisher's scores from line to matrix to plot them
fisherScore = reshape(averageFisher' , 19, 16);

%plot of Fisher's scores.
figure
imagesc('XData',params_spectrogram.freq,'YData',1:16,'CData',fisherScore')
set(gca, 'YTick', 1:16,'YTickLabels', {chanlocs16.labels})
xlabel('Frequencies (Hz)')
axis tight;
ylabel('Channels')
title('Fisher discriminant power of features')
h=colorbar;

%% train on all the data set with the parameters selection

for j=1:size(EpochTraining.MITerm.data,3)
    dataTraining(j,:)=reshape(EpochTraining.MITerm.data(:,:,j)',[1,16*19]); 
end    
[dataTraining, Mu, Sigma]=zscore(dataTraining);
classifier =fitcdiscr(dataTraining(:,indexPower(1:hyperparameters.Nsel)),EpochTraining.MITerm.labels,'discrimtype',  hyperparameters.ClassifierType);

%Online plot
[Trials,SmoothedTotal]=real_online(s{1,length(s)},Mu,Sigma,classifier,indexPower(1:hyperparameters.Nsel));