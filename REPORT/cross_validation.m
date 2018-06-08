function [hyperparameters,MinClassError,ClassError,AUC]=cross_validation(Folds,featuresList,nFold)
%{
Inputs
-Folds: trial-based folds. Need to be organized like this:
        Folds.MITerm.data (vector of nFold cells with inside each cell
        a matrix with n°trials x n°features) and Folds.MITerm (vector of
        nFold cells with inside each cell a vector of n°trials)
-featuresList: range for Nsel (e.g. featuresList=1:25)
-nFold: n° of folds (e.g. 10)
-------
Outputs
-hyperparameters: hyperparameters (Nsel, type of classifier) selected from
                  the best model (i.e. the model with lowest class error)
-MinClassError: class error (averaged over the nFold folds) corresponding
                to the best model
-ClassError: struct containing the mean (over the folds) of the class error
             for the best model (i.e. the model with the lowest average
             class error) and the std
-AUC: struct containing the mean (over the folds) value of the AUC for the 
      best model and the std
%}

global SubjectID
global chanlocs16
global params_spectrogram

Classifier={'linear', 'diaglinear','diagquadratic','quadratic'};


Pfeature = [];
class_error = [];
yPredicted = [];
PosteriorProb = [];
X=[];
Y=[];

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
            
            [X{iFeature,iClassifier,iFold},Y{iFeature,iClassifier,iFold},T{iFeature,iClassifier,iFold},AUC{iFeature,iClassifier,iFold}] = perfcurve(yTest{iFold},PosteriorProb{iFold}(:,2),555);
            
        end
    end
end

%error
MinClassError=0.5;
figure
for iClassifier=1:numel(Classifier)
    class_error4Classifier = squeeze(class_error(:,iClassifier,:));
    
    meanClassError = mean(class_error4Classifier,1);
    stdError=std(class_error4Classifier); %for each classifier  
    [m,i]=min(meanClassError);
    if m<MinClassError
        MinClassError=m;
        StdClassError=stdError(i);
        hyperparameters.Nsel=i;
        hyperparameters.ClassifierType=Classifier{iClassifier};
    end
    plot(featuresList,meanClassError)
    hold on
end
title(sprintf('Classifier and class error -  %s', SubjectID));
legend('Classifier: Linear','Classifier: DiagLinear','Classifier:Diagquadratic','Classifier:Quadratic');
xlabel('Features');
ylabel('Class error');

ClassError.mean=MinClassError;
ClassError.std=StdClassError;

%% perfcurve

if strcmp(hyperparameters.ClassifierType,'linear')||strcmp(hyperparameters.ClassifierType,'quadratic')
    iClassifier=find(strcmp(Classifier, hyperparameters.ClassifierType));
    iFeature=hyperparameters.Nsel;
    XRoc=[];
    YRoc=[];
    AUCRoc=[];
    for iFold=1:10
        XRoc=cat(1,XRoc,(X{iFeature,iClassifier,iFold})');
        YRoc=cat(1,YRoc,(Y{iFeature,iClassifier,iFold})');
        AUCRoc=cat(1,AUCRoc,(AUC{iFeature,iClassifier,iFold})');
    end
    XRocMean=mean(XRoc);
    YRocMean=mean(YRoc);
    YRocStd=std(YRoc);
    AUCMean=mean(AUCRoc);
    AUCStd=std(AUCRoc);
    AUC=[];
    AUC.Mean=AUCMean;
    AUC.Std=AUCStd;
    
    figure
    plotshaded(XRocMean,[YRocMean-YRocStd; YRocMean; YRocMean+YRocStd],'b')
    set(gca,'YLim',[0 1])
    xlabel('False positive rate') 
    ylabel('True positive rate')
    title(sprintf('ROC curve averaged over the folds -  %s', SubjectID))
else
    disp('mean ROC curve not possible for diaglinear/diagquadratic');
    AUC=[];
end

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
title(sprintf('Fisher scores map -  %s', SubjectID))
h=colorbar;



end