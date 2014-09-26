function [calfilename] = RandomPredictVoiceConfusionMatrices_Cal(MatfilePath)
Bootmax=1000;%Number of random matrices to compute per cell
MAT=load(MatfilePath);
    
%% Retrieve the best Individual confusion matrix of that cell
Winsize=MAT.winSize;
MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
optWin_CT=Winsize(MAXWinCT);
if length(optWin_CT)>1
    optWin_CT=optWin_CT(1);%by default, take the first one
end
VocConfMat = MAT.confusionMatrix{find(Winsize==optWin_CT)};
% convert the joint probability matrix to an event matrix
SF_ConfMat = VocConfMat.*MAT.neventMatrix(find(Winsize==optWin_CT));
% retrieve the vocalization type order
VoiceTypeSel = MAT.VoiceTypeSel;
% retrieve the number of different categories
StimTypeCM=unique(VoiceTypeSel);
% Make sure stimType Background is at the end of the category list
BG_Ind = find(strcmp(StimTypeCM, 'BG'));
StimTypeCM = [StimTypeCM(1:(BG_Ind-1)); StimTypeCM((BG_Ind + 1):end); StimTypeCM(BG_Ind)];
NstimTypeCM=length(StimTypeCM);

%% find the vector of indices for spike trains obtained during sound (suppress those obtain during background) in the confusion matrix SF_ConfMat
IndBackground = find(strcmp(VoiceTypeSel, 'BG'));
IndVocOnly = setdiff(1:length(VoiceTypeSel), IndBackground);

%% Calculate bootmax random matrices with Background intact and bootmax matrices with Background shuffled
CM_IV_RandBG = cell(Bootmax, 1);
List_VocRandBG = cell(Bootmax, 1);
Nb_VocPerCatBG = cell(Bootmax, 1);
CM_IV_Rand = cell(Bootmax, 1);
List_VocRand = cell(Bootmax, 1);
Nb_VocPerCat = cell(Bootmax, 1);
for bb=1:Bootmax
    bb
    rng('shuffle'); % seeds the random number generator based on the current time
    VocTypeSel_randBG=[VoiceTypeSel(IndVocOnly(randperm(length(IndVocOnly)))); VoiceTypeSel(IndBackground)];
    RR=0;
    VocRandBG = zeros(1,length(VoiceTypeSel));
    NVPC = zeros(NstimTypeCM,1);
    for vtR=1:NstimTypeCM
        stR=StimTypeCM(vtR);
        selectorR=strcmp(VocTypeSel_randBG, stR);
        selectorIndR=find(selectorR);
        VocRandBG(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
        RR=RR+length(selectorIndR);
        NVPC(vtR) = length(selectorIndR);
    end
    confusion_matrix_vocalizationsRandBG = SF_ConfMat(:, VocRandBG);
    confusion_matrix_vocalizationsRandBG = confusion_matrix_vocalizationsRandBG./sum(sum(confusion_matrix_vocalizationsRandBG));
    CM_IV_RandBG{bb} = confusion_matrix_vocalizationsRandBG;
    List_VocRandBG{bb} = VocRandBG;
    Nb_VocPerCatBG{bb} = NVPC;
    
    VocTypeSel_rand=VoiceTypeSel(randperm(length(VoiceTypeSel)));
    RR=0;
    VocRand = zeros(1,length(VoiceTypeSel));
    NVPC = zeros(NstimTypeCM,1);
    for vtR=1:NstimTypeCM
        stR=StimTypeCM(vtR);
        selectorR=strcmp(VocTypeSel_rand, stR);
        selectorIndR=find(selectorR);
        VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
        RR=RR+length(selectorIndR);
        NVPC(vtR) = length(selectorIndR);
    end
    confusion_matrix_vocalizationsRand = SF_ConfMat(:, VocRand);
    confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand));
    CM_IV_Rand{bb} = confusion_matrix_vocalizationsRand;
    List_VocRand{bb} = VocRand;
    Nb_VocPerCat{bb} = NVPC;
end
RandMat.CM_IV_RandBG = CM_IV_RandBG;
RandMat.List_VocRandBG = List_VocRandBG;
RandMat.Nb_VocPerCatBG = Nb_VocPerCatBG;
RandMat.CM_IV_Rand = CM_IV_Rand;
RandMat.List_VocRand = List_VocRand;
RandMat.Nb_VocPerCat = Nb_VocPerCat;
RandMat.subject = MAT.subject;
RandMat.originalfile = MatfilePath;
RandMat.VoiceTypeSel=VoiceTypeSel; %added on 08/15/2014
[Path, Matfile] = fileparts(MatfilePath);

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',MAT.subject,['RandPVoi_' Matfile(9:end) '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['RandPVoi_' Matfile(9:end) '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',MAT.subject,['RandPVoi_' Matfile(9:end) '.mat']);
end

save(calfilename, '-struct', 'RandMat');
fprintf(1,'done making calculus on %s\nData save under %s\n', MatfilePath, calfilename);
end