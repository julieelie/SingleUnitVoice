function []=Cal_MIperarea_ConfVoice(Voicematfile)
FigFlag=0;

%% Cell aray containing the different classifications of the stims studied
StimTypeCM=cell(1,4);
VoiceTypeSel=StimTypeCM;
StimTypeCM{1}=cell(11,1);
StimTypeCM{2}={'A' 'C' 'B'};
StimTypeCM{3}={'F' 'M' 'G'};
StimTypeCM{4}={'AF' 'AM' 'CF' 'CM' 'BG'};

%% Load file
MAT=load(Voicematfile);

%% Calculate values of MI per area for the observed individual
% vocalization matrices
MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
if length(MAXWinCT)>1 %for highly random matrices, values of MI can be negative and equal whatever the window size then there is no max
    MAXWinCT=MAXWinCT(1);
end
sprintf('Calculate values of MI per area in observed Matrix\n')
Mat = MAT.confusionMatrix{MAXWinCT}; % we are choosing the window size that gives the best value of MI conf in the CT matrix
VocTypeConfMat = MAT.confusionMatrixCT{MAXWinCT};

if FigFlag>0
    figure(1)
    subplot(1,2,1)
    imagesc(Mat)
    colorbar;
    xlabel('Predicted vocalization');
    ylabel('Actual vocalization');
    title(sprintf('Single vocalizations\n'));
    subplot(1,2,2)
    imagesc(VocTypeConfMat)
    colorbar;
    xlabel('Predicted vocalization');
    ylabel('Actual vocalization');
    title(sprintf('Voc Type \n'));
end

% construct the cell array of indices for each name category
VoiceTypeSel{1} = MAT.VoiceTypeSel;
StimTypeCM{1}=unique(VoiceTypeSel{1});
IBG = find(strcmp(StimTypeCM{1}, 'BG'));
StimTypeCM{1} = [StimTypeCM{1}(1:(IBG-1)); StimTypeCM{1}((IBG+1):end); StimTypeCM{1}(IBG)];
NstimTypeCM=length(StimTypeCM{1});
cat1 = cell(NstimTypeCM,1);
for cc = 1:NstimTypeCM
    cat1{cc}=find(strcmp(StimTypeCM{1}(cc), VoiceTypeSel{1}));
end

% construct the cell array of indices for age category
VoiceTypeSel{2}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{2}{vv}=VoiceTypeSel{1}{vv}(1);
end
NstimTypeCM=length(StimTypeCM{2});
cat2 = cell(NstimTypeCM,1);
for cc = 1:NstimTypeCM
    cat2{cc}=find(strcmp(StimTypeCM{2}(cc), VoiceTypeSel{2}));
end

% construct the cell array of indices for sex category
VoiceTypeSel{3}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{3}{vv}=VoiceTypeSel{1}{vv}(2);
end
NstimTypeCM=length(StimTypeCM{3});
cat3 = cell(NstimTypeCM,1);
for cc = 1:NstimTypeCM
    cat3{cc}=find(strcmp(StimTypeCM{3}(cc), VoiceTypeSel{3}));
end

% construct the cell array of indices for agesex category
VoiceTypeSel{4}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{4}{vv}=VoiceTypeSel{1}{vv}(1:2);
end
NstimTypeCM=length(StimTypeCM{4});
cat4 = cell(NstimTypeCM,1);
for cc = 1:NstimTypeCM
    cat4{cc}=find(strcmp(StimTypeCM{4}(cc), VoiceTypeSel{4}));
end

CAT={cat1 cat2 cat3 cat4};
%Calculate the values of MI in the different areas of the IV
%confusion matrix
[ mi_tot, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea_voice(Mat, CAT, FigFlag);

%% Store values in a structure
MI_Calculations.MI_tot=mi_tot;
MI_Calculations.MI_diag_uni_cat=mi_diag_uni_cat;
MI_Calculations.MI_real_error_uni=mi_real_error_uni;
MI_Calculations.CAT=CAT;
MI_Calculations.VoiceTypeSel=VoiceTypeSel;
MI_Calculations.StimTypeCM=StimTypeCM;

%% Add that new structure to the file
save(Voicematfile, 'MI_Calculations', '-append');
clear MAT MI_Calculations CAT mi_tot mi_diag_uni_cat mi_real_error_uni VoiceTypeSel StimTypeCM Mat
end
