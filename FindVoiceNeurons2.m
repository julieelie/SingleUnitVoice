
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')
NbCell_estimate = 1500;
FigFlag =0;

%% Cell aray containing the different classifications of the stims studied
StimTypeCM_predict=cell(1,4);
StimTypeCM_predict{1}=cell(15,1);%there must be max 15 different names per confusion matrix
StimTypeCM_predict{2}={'A' 'C' 'B'};
StimTypeCM_predict{3}={'F' 'M' 'G'};
StimTypeCM_predict{4}={'AF' 'AM' 'CF' 'CM' 'BG'};


%% Retrieve spike shape info
%cd /Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/SpikeShape
cd /auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode/SpikeShape
Ldir=pwd;
SpikeShape=load(strcat(Ldir,'/spikeShapeResults.mat'),'allNames', 'indG5', 'spikeType');

%% Retrieve Histology info
cd /auto/fhome/julie/Documents/Histo
% cd /Volumes/FREDERIC/HistoJulie/Stacks
AnaLdir=pwd;
AnatTxt=dir(fullfile(AnaLdir, 'List_h5files*Histo.txt'));
Anat_root = cell(NbCell_estimate, 5);
nbTxt = length(AnatTxt);
Lanat = 0;
for tt=1:nbTxt
    textfile = AnatTxt(tt).name;
    fid = fopen(fullfile(AnaLdir, textfile));
    dataAnat = textscan(fid, '%s\t%s\t%f\t%f\t%s');
    anatnb = size(dataAnat{1}, 1);
    for ll = 1: anatnb
        Lanat=Lanat+1;
        for col=1:5
            Anat_root{Lanat,col}=dataAnat{col}(ll);
        end
    end
end
Anat_root = Anat_root(1:Lanat,:);


%% Setup some output variables
% cell containing the path of each unit
List_matfilepath = cell(NbCell_estimate, 1);

% cell containing the anatomy info of each unit
List_anat = cell(NbCell_estimate, 5);

% vectors containing optimal windows for each unit
optWin=zeros(NbCell_estimate,1); % based on the sound file matrix
optWin_BGwo=zeros(NbCell_estimate,1); % based on the sound file matrix without BG
optWin_CT= zeros(NbCell_estimate,1); % based on the calls type matrix

% vectors containing mutual information for each unit
MI_confB=zeros(NbCell_estimate,9);
MI_confBBGwo=zeros(NbCell_estimate,9);
MI_confBCT=zeros(NbCell_estimate,9);
MI_confBCTBGwo=zeros(NbCell_estimate,9);

% vector containing the mutual information of the best matrix for each unit
MIB_OptCT = zeros(NbCell_estimate,1);
MIBCT_OptCT = zeros(NbCell_estimate,1);
MIB_Opt = zeros(NbCell_estimate,1);
MIBCT_Opt = zeros(NbCell_estimate,1);

% vector containing the size of the confusion matrices for each unit
ConfSize = zeros(NbCell_estimate,2);
ConfSizeCT = zeros(NbCell_estimate,2);

% Vectors containing the mutual information per area for all units
MI_tot=zeros(NbCell_estimate,1);
MI_diag_uni_cat=zeros(NbCell_estimate,4);% the four column correspond to the four categories (names, age, sex, ageXsex)
MI_real_error_uni=zeros(NbCell_estimate,4);

AvMI_tot_Rand = zeros(NbCell_estimate,2); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
SDMI_tot_Rand = zeros(NbCell_estimate,2);
AVMI_diag_uni_cat_Rand=cell(2,1);
AVMI_diag_uni_cat_Rand{1}=nan(NbCell_estimate,length(StimTypeCM_predict));
AVMI_diag_uni_cat_Rand{2}=AVMI_diag_uni_cat_Rand{1};
SDMI_diag_uni_cat_Rand=AVMI_diag_uni_cat_Rand;
MI_uni_diag_cat_Rand = cell(NbCell_estimate,1);
MI_uni_diag_cat_RandBG = cell(NbCell_estimate,1);


PCC_cat = cell(length(StimTypeCM_predict),1);
for iss=1:length(PCC_cat);
    PCC_cat{iss} = cell(NbCell_estimate,1);
end
IS3=PCC_cat;
IS2 = zeros(NbCell_estimate,length(StimTypeCM_predict));% the four column correspond to the four categorisatios (names, age, sex, ageXsex)
II2 = zeros(NbCell_estimate,length(StimTypeCM_predict));% the four column correspond to the four categorisatios (names, age, sex, ageXsex)

% Vector of the mean spiking rate for each unit calculated from WholeVoc files on first extracts
MeanSR = zeros(NbCell_estimate,1);

% Vector of the spike shape and sort quality
Spike_shape = zeros(NbCell_estimate,1);

% Vector of the code for the subject
SUBJ = zeros(NbCell_estimate,2);
SUBJECTS = {'YelBlu6903F' 'BlaBro09xxF' 'LblBlu2028M' 'GreBlu9508M' 'WhiBlu5396M' 'WhiWhi4522M'};

%Vector of the code for the anatomical zones
ZONES = zeros(NbCell_estimate,2);
ZONES_List = {'CMM', 'CML', 'L1', 'L2A', 'L2B', 'L3', 'L', 'NCM', 'HP'};
            

% Matrix containing the invariance value for each cell and each category
H_invariance = cell(length(StimTypeCM_predict),1);
H_invariance{1}=nan(NbCell_estimate,length(StimTypeCM_predict{1}));%there must be max 15 different names per confusion matrix
H_invariance{2}=nan(NbCell_estimate,length(StimTypeCM_predict{2}));%2 ages C and A + B for background stims
H_invariance{3}=nan(NbCell_estimate,length(StimTypeCM_predict{3}));%2 sexes F and M + B for background stims
H_invariance{4}=nan(NbCell_estimate,length(StimTypeCM_predict{4}));%2 ages C and A x 2 sexes (F and M) + B for background stims
H_invariancemax = H_invariance;
II_percat = H_invariance;

% cell of Vectors containing the number of vocalizations per category, the total number of vocalizations and the proba of each category for all units
NVocPerCat = cell(length(StimTypeCM_predict),1);
VocCatNames = NVocPerCat;
for cl=1:length(StimTypeCM_predict)
    NVocPerCat{cl}=nan(NbCell_estimate,length(StimTypeCM_predict{cl}));
    VocCatNames{cl}=cell(NbCell_estimate,1);
end
NVoc = zeros(NbCell_estimate,1);
ProbaVocPerCat = NVocPerCat;

%cell of Vectors containing the mean spike rate per category per unit
SpikeRate_Cat = cell(length(StimTypeCM_predict),1);
for cl=1:length(StimTypeCM_predict)
    SpikeRate_Cat{cl}=nan(NbCell_estimate,length(StimTypeCM_predict{cl})-1); %We don't need the spike rate during background that's why: StimTypeCM-1
end

%cell of Matrices containing the number of trial per vocalization per unit
NbTrial_voc = NVocPerCat;

%Cell arrays containing the mean spike rate and the TDT name for each vocalization for each
%unit
SpikeRate_Voc = cell(NbCell_estimate,1);
TDT_names = cell(NbCell_estimate,1);

% Cell arrays containing the real StimTypeCM and VoiceTypeSel of each unit
StimTypeCM_allUnits = cell(NbCell_estimate,1);
VoiceTypeSel_allUnits = cell(NbCell_estimate,1);

%% Start to loop through individuals and files and extract information
cd /auto/k6/julie/matfile
input_dir=pwd;
Subjects = dir(input_dir);
ii=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11 && ~strcmp(Indiv, 'WhiBlu5396M')% we are not keeping data from WhiBlu5396M here because there is only one individual per sex and age category
        Idir=fullfile(input_dir, Indiv);
        
        sprintf('Harvesting data of %s\n', Indiv)
        % retrieve Confusion matrices files
        matfiles=dir(fullfile(Idir,'ConfVoi*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        Indices=find(SS_Ind);
        LM=length(Indices);
        
        %retrieve Random confusion matrices files
        Randmatfiles=dir(fullfile(Idir, 'RandPVoi*.mat'));
        rm=length(Randmatfiles);
        SS_Ind=zeros(rm,1);
        for ff=1:rm
            if ~isempty(strfind(Randmatfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        IndicesRM=find(SS_Ind);
        RM=length(IndicesRM);
        
        % retrieve WholeVoc files to calculate Mean rate
        WVmatfiles = dir(fullfile(Idir, 'WholeVoc*.mat'));
        wm=length(WVmatfiles);
        SS_Ind=zeros(wm,1);
        for ff=1:wm
            if ~isempty(strfind(WVmatfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        IndicesWV=find(SS_Ind);
        WV=length(IndicesWV);
        
         % retrieve FirstVoc files to extract the number
        % of trials for all stims within same category
        FVmatfiles = dir(fullfile(Idir, 'FirstVoc*.mat'));
        fm=length(FVmatfiles);
        FSS_Ind=zeros(fm,1);
        for ff=1:fm
            if ~isempty(strfind(FVmatfiles(ff).name, 'ss'))
                FSS_Ind(ff)=1;
            end
        end
        IndicesFV=find(FSS_Ind);
        FV=length(IndicesFV);
    
    
        
        % check that we have the same number of files for all of them
        if LM~=RM
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of randfiles (%d)\n', LM, RM)
            break
        end
        if LM~=WV
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of WholeVocfiles (%d)\n', LM, WV)
            break
        end
        if LM~=FV
            sprintf('WARNING: the nb of observed matrices (%d) is different from the nb of FirstVocfiles (%d)\n', LM, FV)
            break
        end
        
  %% loop through files and calculate all variables
        for hh=1:LM
            Matfile = matfiles(Indices(hh)).name;
            sprintf('Loading %s\n', Matfile)
            MatfilePath=fullfile(Idir, Matfile);
            MAT=load(MatfilePath);
            
            ii=ii+1;
            %% keep track of the site ID
            List_matfilepath{ii}=MatfilePath;
            
            %% Code the subject sex and identity
            sprintf('Subject ID\n')
            if strcmp(MAT.subject, 'BlaBro09xxF')
                SUBJ(ii,:)=[1 0];
            elseif strcmp(MAT.subject, 'LblBlu2028M')
                SUBJ(ii,:)=[2 1];
            elseif strcmp(MAT.subject, 'GreBlu9508M')
                SUBJ(ii,:)=[3 1];
            elseif strcmp(MAT.subject, 'WhiBlu5396M')
                SUBJ(ii,:)=[4 1];
            elseif strcmp(MAT.subject, 'WhiWhi4522M')
                SUBJ(ii,:)=[5 1];
            end
            
            %% Construct anatomy vectors
            sprintf('Anatomy\n')
            for ff = 1:size(Anat_root, 1)
                TempID = cell2mat(Anat_root{ff, 1});
                if strcmp(Matfile(9:end-8), TempID(1:end-3))
                    List_anat(ii, :) = Anat_root(ff, :);
                end
            end

            if strcmp(List_anat{ii,2}, 'CMM')
                ZONES(ii,1)=1;
            elseif strcmp(List_anat{ii,2}, 'CML')
                ZONES(ii,1)=2;
            elseif strcmp(List_anat{ii,2}, 'L1')
                ZONES(ii,1)=3;
            elseif strcmp(List_anat{ii,2}, 'L2A')
                ZONES(ii,1)=4;
            elseif strcmp(List_anat{ii,2}, 'L2B')
                ZONES(ii,1)=5;
            elseif strcmp(List_anat{ii,2}, 'L3')
                ZONES(ii,1)=6;
            elseif strcmp(List_anat{ii,2}, 'L')
                ZONES(ii,1)=7;
            elseif strcmp(List_anat{ii,2}, 'NCM')
                ZONES(ii,1)=8;
            elseif strcmp(List_anat{ii,2}, 'HP')
                ZONES(ii,1)=9;
            end
            if strcmp(List_anat{ii,5}, 'L')%lef or right hemisphere
                ZONES(ii,2)=1;
            end
            
            
            %% Collect info about the observed matrices
            sprintf('Observed Confusion Matrices\n')
            Winsize=MAT.winSize;
            MAXWin=find(MAT.mi_confusionB==max(MAT.mi_confusionB));
            MAXWinCT=find(MAT.mi_confusionCT==max(MAT.mi_confusionCT));
            MAXWin_BGwo=find(MAT.mi_confusionB_BGwo==max(MAT.mi_confusionB_BGwo));
            if length(MAXWinCT)>1 %for highly random matrices, values of MI can be negative and equal whatever the window size then there is no max
                MAXWinCT=MAXWinCT(1);
                MAXWin=MAXWin(1);
                MAXWin_BGwo=MAXWin_BGwo(1);
            end
            optWin(ii)=Winsize(MAXWin);
            optWin_BGwo(ii)=Winsize(MAXWin_BGwo);
            optWin_CT(ii)=Winsize(MAXWinCT);
            MI_confB(ii,:)=MAT.mi_confusionB;
            MI_confBBGwo(ii,:)=MAT.mi_confusionB_BGwo;
            MI_confBCT(ii,:)=MAT.mi_confusionCT;
            MI_confBCTBGwo(ii,:)=MAT.mi_confusionCT_BGwo;
            ConfSize(ii,:) = size(MAT.confusionMatrix{1});
            ConfSizeCT(ii,:) = size(MAT.confusionMatrixCT{1});
            MIB_OptCT(ii) = MI_confB(ii,find(Winsize==optWin_CT(ii)));
            MIB_Opt(ii) = MI_confB(ii,find(Winsize==optWin(ii)));
            MIBCT_OptCT(ii) = MI_confBCT(ii,find(Winsize==optWin_CT(ii)));
            MIBCT_Opt(ii) = MI_confBCT(ii,find(Winsize==optWin(ii)));
            
            %Store the values of MI in the different areas of the IV
            %confusion matrix
            MI_tot(ii)=MAT.MI_Calculations.MI_tot;
            MI_diag_uni_cat(ii,:)=MAT.MI_Calculations.MI_diag_uni_cat;
            MI_real_error_uni(ii,:)=MAT.MI_Calculations.MI_real_error_uni;
            CAT=MAT.MI_Calculations.CAT;
            
            %% Load some needed variables (confusion matrix, stim names...) for following calculations
            Mat = MAT.confusionMatrix{find(Winsize==optWin_CT(ii))};
            VoiceTypeSel=MAT.MI_Calculations.VoiceTypeSel;
            StimTypeCM=MAT.MI_Calculations.StimTypeCM;
            VoiceTypeSel_allUnits{ii}=VoiceTypeSel;
            StimTypeCM_allUnits{ii} = StimTypeCM;
           
            
            %% Calculate the indices of invariance based on entropy for each cell and each category
            sprintf('Index of Invariance per category\n')
            for cl=1:length(CAT);
                cat = CAT{cl};
                NstimTypeCM=length(StimTypeCM{cl});
                for cc = 1:NstimTypeCM
                    Nb_stim = length(cat{cc});
                    Psum=0;
                    for rr=1:Nb_stim
                        rro = cat{cc}(rr);
                        for co = 1:Nb_stim
                            cco = cat{cc}(co);
                            Psum=Psum + Mat(rro,cco);
                        end
                    end
                    Mat_for_entropy = reshape(Mat(cat{cc}, cat{cc}),Nb_stim.*Nb_stim,1);
                    Mat_for_entropy = Mat_for_entropy./Psum;%some values are going to be NAN but we don't worry too much since these are non discriminant units
                    MZeroIndices=find(Mat_for_entropy==0);
                    Mat_for_entropy(MZeroIndices)=1;%to make sure entropy is zero when p=0 because 1.*log2(1)=0
                    H_invariance{cl}(ii,cc) = sum(-Mat_for_entropy.*log2(Mat_for_entropy));
                    H_invariancemax{cl}(ii,cc) = -log2(1./(Nb_stim.*Nb_stim));
                end
                II_percat{cl}(ii,:) = H_invariance{cl}(ii,:)./H_invariancemax{cl}(ii,:);
            end
            %% Calculate the alternative Index of Invariance II2
            % For all the category zones calculate observed entropy then max entropy if the joint probability was smoothed within each row 
            sprintf('Index of Invariance\n')
            zero_ind = Mat == 0;
            prob_matrix_for_entropy = Mat;
            prob_matrix_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
            for cl=1:length(CAT);
                cat=CAT{cl};
                Hobs_Inv = 0;
                Hmax_Inv = 0;
                nstim = size(prob_matrix_for_entropy, 1);
                Ncat = length(cat);
                for aa = 1:nstim
                    ProbaCat = 0;
                    for nc = 1:Ncat
                            icat = cat{nc};
                            if ~isempty(intersect(aa, icat))
                                break
                            end
                    end
                    nstim_local = length(icat);
                    for bb = 1:nstim_local;
                        cc=icat(bb);
                        Hobs_Inv = Hobs_Inv + prob_matrix_for_entropy(aa,cc) * log2(prob_matrix_for_entropy(aa,cc));
                        ProbaCat = ProbaCat + prob_matrix_for_entropy(aa,cc);
                    end
                    Hmax_Inv = Hmax_Inv + ProbaCat*log2(ProbaCat/nstim_local);
                end
                II2(ii,cl) = Hobs_Inv/Hmax_Inv;
            end

            %% Calculate an index of selectivity for the different call categories for that unit
            sprintf('Indices of Selectivity\n')
            for cl=1:length(CAT)
                % First calculate the CT Matrix
                NstimTypeCM= length(StimTypeCM{cl});
                confusion_matrix_CallType = zeros(NstimTypeCM, NstimTypeCM);
                for vtR=1:NstimTypeCM
                    stR=StimTypeCM{cl}(vtR);
                    selectorR=strcmp(VoiceTypeSel{cl}, stR);
                    selectorIndR=find(selectorR);
                    for vtC = 1:NstimTypeCM
                        stC=StimTypeCM{cl}(vtC);
                        selectorC=strcmp(VoiceTypeSel{cl}, stC);
                        selectorIndC=find(selectorC);
                        confusion_matrix_CallType(vtR,vtC)=sum(sum(Mat(selectorIndR, selectorIndC)));
                    end
                end

                %Convert joint probabilities to conditional probailities
                ProbaCat = sum(confusion_matrix_CallType, 2);
                repProbaCat = repmat(ProbaCat, 1, size(confusion_matrix_CallType,2));
                confusion_matrix_CallType_cond = confusion_matrix_CallType ./ repProbaCat;

                % isolate the diagonal vectors containing the percentage of correct
                % classification of each call category
                Vect_class_cond = diag(confusion_matrix_CallType_cond,0);
                PCC_cat{cl}{ii} = Vect_class_cond; %here we still have the background

                % calculating IS index of selectivity based on the entropy of the
                % diagonal of the CT matrix compare to the entropy of a uni diag
                % CT matrix without taking the background into account
                Vect_class = diag(confusion_matrix_CallType,0);
                Vect_class2 = Vect_class./sum(Vect_class);
                Hobs2 = 0;
                Hmax2 = -log2(NstimTypeCM-1);
                zero_ind = Vect_class == 0;
                Vect_class_for_entropy2 = Vect_class2;
                Vect_class_for_entropy2(zero_ind) = 1; 
                for nt = 1:(NstimTypeCM-1)
                    Hobs2 = Hobs2 + Vect_class_for_entropy2(nt)*log2(Vect_class_for_entropy2(nt));
                end
                IS2(ii,cl) = 1- Hobs2/Hmax2;% Normalized probabilities so to calculate more accurately the entropy in the diagonal because entropy calculation is complete when sum(proba)=1


                % calculating IS3 index of selectivity3 on the matrix diagonal
                % without taking the background
                IS3_temp=zeros((NstimTypeCM-1),1);
                for vt=1:(NstimTypeCM-1)
                    Vect_class_temp = Vect_class_cond(1:(end-1));
                    Vect_class_temp(vt)=[];
                    IS3_temp(vt)=log2(Vect_class_cond(vt)/mean(Vect_class_temp));
                end
                IS3{cl}{ii}=IS3_temp;
            end


 %% Upload the random matrices and calculate MI per area
            sprintf('Collect MI per area for random matrices\n')
            for rr=1:RM
                RandMatfile = Randmatfiles(IndicesRM(rr)).name;
                if strcmp(Matfile(9:end), RandMatfile(10:end))
                    RMAT = load(fullfile(Idir, RandMatfile));
                    break
                end
            end
            
            
            AvMI_tot_Rand(ii,:) =RMAT. MI_Calculations.AvMI_tot_Rand;
            SDMI_tot_Rand(ii,:) = RMAT.MI_Calculations.SDMI_tot_Rand;
            fprintf('the MI_confusion of the IV Matrix is %f for that cell\nThe average MI of random matrices is %f+/-%f\n', MI_tot(ii), AvMI_tot_Rand(ii,2),SDMI_tot_Rand(ii,2))
            for rt=1:2
                AVMI_diag_uni_cat_Rand{rt}(ii,:)=mean(RMAT. MI_Calculations.MI_uni_diag_cat_RandBG,1);
                SDMI_diag_uni_cat_Rand{rt}(ii,:)=std(RMAT. MI_Calculations.MI_uni_diag_cat_RandBG,1);
            end
            %MI_uni_diag_cat_Rand{ii} = RMAT.MI_Calculations.MI_uni_diag_cat_Rand;
            MI_uni_diag_cat_RandBG{ii} = RMAT.MI_Calculations.MI_uni_diag_cat_RandBG;
            
            
%% find the mean spike rate of the unit
            sprintf('Mean spike rate\n')
            for ww=1:WV
                WVMatfile = WVmatfiles(IndicesWV(ww)).name;
                if strcmp(Matfile(9:end), WVMatfile(10:end))
                    WMAT = load(fullfile(Idir, WVMatfile));
                    break
                end
            end
            FirstWholeVoc = find(WMAT.Voc_orders==1);
            MeanSR_sections=cell2mat(WMAT.MeanRate(FirstWholeVoc));
            MeanSR(ii)=nanmean(MeanSR_sections);
            
%% find the mean spike rate for each first vocalization
            sprintf('Individual Voc Mean spike rate\n')
            SpikeRate_Voc{ii} = MeanSR_sections;
            TDT_names{ii} = WMAT.TDT_wavfiles(FirstWholeVoc);
            
%% find the mean spike rate of the unit per call category, the number of voc per category and the proba of a given category
            sprintf('Category Mean spike rate\n')
            WavfilesFirstSect = WMAT.Original_wavfiles(FirstWholeVoc);
            VocTypeFirstSect=cell(length(StimTypeCM),1);
            for vv=1:length(FirstWholeVoc)
                [path,filename,ext]=fileparts(WavfilesFirstSect{vv});
                VocTypeFirstSect{1}{vv} = strcat(filename(13),filename(12),filename(14),'_',filename(1:10));
                VocTypeFirstSect{2}{vv} = filename(13);
                VocTypeFirstSect{3}{vv} = filename(12);
                VocTypeFirstSect{4}{vv} = filename(12:13);
            end 
            
            NVoc(ii) = size(MAT.VoiceTypeSel,1);
            for cl=1:length(StimTypeCM)
                NstimTypeCM=length(StimTypeCM{cl});
                Icat_all = cell(NstimTypeCM-1, 1);
                for cc=1:(NstimTypeCM-1)%There is no value of mean spike rate for backgroung here
                    cat = StimTypeCM{cl}(cc);
                    Icat_all{cc} = find(strcmp(VocTypeFirstSect{cl},cat));
                    SpikeRate_Cat{cl}(ii,cc) = nanmean(MeanSR_sections(Icat_all{cc}));
                end
                VocCatNames{cl}{ii} = StimTypeCM{cl};
                for cc=1:NstimTypeCM
                    cat = StimTypeCM{cl}(cc);
                    NVocPerCat{cl}(ii,cc) = sum(strcmp(VoiceTypeSel{cl}, cat));
                    ProbaVocPerCat{cl}(ii,cc) = NVocPerCat{cl}(ii,cc)./NVoc(ii);
                end
            
            end
            
%% Find the sum of trials for all vocalizations per category
            sprintf('Nb Trials for all vocalizations within each category\n')
            for fw=1:FV
                FVMatfile = FVmatfiles(IndicesFV(fw)).name;
                if strcmp(Matfile(9:end), FVMatfile(10:end))
                    FMAT = load(fullfile(Idir, FVMatfile));
                    break
                end
            end
            VocTypeFirstVoc=cell(length(StimTypeCM),1);
            for vv=1:length(FMAT.Original_wavfiles)
                    [path,filename,ext]=fileparts(FMAT.Original_wavfiles{vv});
                    VocTypeFirstVoc{1}{vv} = strcat(filename(13),filename(12),filename(14),'_',filename(1:10));
                    VocTypeFirstVoc{2}{vv} = filename(13);
                    VocTypeFirstVoc{3}{vv} = filename(12);
                    VocTypeFirstVoc{4}{vv} = filename(12:13); 
            end
            for cl=1:length(StimTypeCM)
                NstimTypeCM = length(StimTypeCM{cl});
                % retrieve values for non background sounds
                for cc=1:(NstimTypeCM-1)
                    cat = StimTypeCM{cl}(cc);
                    Icat_local = find(strcmp(VocTypeFirstVoc{cl},cat));
                    Nsections=length(Icat_local);
                    for section = 1:Nsections
                        NbTrial_voc{cl}(ii,cc) = NbTrial_voc{cl}(ii,cc) + length(FMAT.Trials{Icat_local(section)});
                    end
                end
                %retrieve values for background sounds
                BGInd = find(strcmp(MAT.VoiceTypeSel, 'BG'));
                TDTnames = MAT.TDTwav(BGInd);
                for bg=1:length(TDTnames)
                    TDTname = TDTnames(bg);
                    ITDTname = find(strcmp(FMAT.TDT_wavfiles, TDTname));
                    NbTrial_voc{cl}(ii,NstimTypeCM) = NbTrial_voc{cl}(ii,NstimTypeCM) + length(FMAT.Trials_BG{ITDTname});
                end
            end            
%% find the spike shape of the unit
            sprintf('Spike shape\n')
            uu=1;
            while strcmp(Matfile(9:end-4), SpikeShape.allNames{uu}(31:end-3))==0
                uu=uu+1;
            end
            if ~isempty(find(SpikeShape.indG5==uu))
                sshape=SpikeShape.spikeType(find(SpikeShape.indG5==uu));
                if strcmp(sshape, 'unclassified')
                    Spike_shape(ii)=1;
                elseif strcmp(sshape, 'large')
                    Spike_shape(ii)=2;
                elseif strcmp(sshape, 'narrow')
                    Spike_shape(ii)=3;
                elseif strcmp(sshape, 'wide')
                    Spike_shape(ii)=4;
                end
            end
        end
    end
end

%% Compile results into structures
sprintf('Compile results\n')
List_matfilepath=List_matfilepath(1:ii);
List_anat=List_anat(1:ii,:);
SUBJ = SUBJ(1:ii,:);
ZONES = ZONES(1:ii,:);
Spike_shape = Spike_shape(1:ii);
StimTypeCM_allUnits = StimTypeCM_allUnits(1:ii);
VoiceTypeCell_allUnits = VoiceTypeCell_allUnits(1:ii);


Conf.optWin=optWin(1:ii);
Conf.optWin_CT=optWin_CT(1:ii);
Conf.optWin_BGwo=optWin_BGwo(1:ii);

Conf.MI_confB=MI_confB(1:ii,:);
Conf.MI_confBBGwo=MI_confBBGwo(1:ii,:);
Conf.MI_confBCT=MI_confBCT(1:ii,:);
Conf.MI_confBCTBGwo=MI_confBCTBGwo(1:ii,:);
Conf.ConfSize = ConfSize(1:ii, :);
Conf.ConfSizeCT = ConfSizeCT(1:ii, :);
Conf.MIB_OptCT = MIB_OptCT(1:ii);
Conf.MIB_Opt = MIB_Opt(1:ii);
Conf.MIBCT_OptCT = MIBCT_OptCT(1:ii);
Conf.MIBCT_Opt = MIBCT_Opt(1:ii);

MI_perArea.MI_tot=MI_tot(1:ii);
MI_perArea.MI_diag_uni_cat=MI_diag_uni_cat(1:ii,:);
MI_perArea.MI_real_error_uni=MI_real_error_uni(1:ii,:);

MI_perArea.AvMI_tot_Rand = AvMI_tot_Rand(1:ii,:); % first column is for rand matrices, second column for rand matrices with Background sections not shuffled
MI_perArea.SDMI_tot_Rand = SDMI_tot_Rand(1:ii, :);
MI_perArea.AVMI_diag_uni_cat_RandBG=AVMI_diag_uni_cat_Rand{2}(1:ii,:);
MI_perArea.SDMI_diag_uni_cat_RandBG=SDMI_diag_uni_cat_Rand{2}(1:ii,:);
MI_perArea.MI_uni_diag_cat_Rand = MI_uni_diag_cat_Rand(1:ii,1);
MI_perArea.MI_uni_diag_cat_RandBG = MI_uni_diag_cat_RandBG(1:ii,1);


Invariance.Globalvalue = II2(1:ii,:);
for cl=1:length(StimTypeCM)
    H_invariance{cl}=H_invariance{cl}(1:ii,:);
    H_invariancemax{cl}=H_invariancemax{cl}(1:ii,:);
    II_percat{cl}=II_percat{cl}(1:ii,:);
    PCC_cat{cl}=PCC_cat{cl}(1:ii);
    IS3{cl}=IS3{cl}(1:ii);
    NVocPerCat{cl} = NVocPerCat{cl}(1:ii,:);
    ProbaVocPerCat{cl} = ProbaVocPerCat{cl}(1:ii,:);
    VocCatNames{cl} = VocCatNames{cl}(1:ii);
    NbTrial_voc{cl} = NbTrial_voc{cl}(1:ii,:);
    SpikeRate_Cat{cl} = SpikeRate_Cat{cl}(1:ii,:);
end
IS2 = IS2(1:ii,:);
Invariance.HPerCat = H_invariance;
Invariance.HPerCatmax = H_invariancemax;
Invariance.IIPerCat = II_percat;
Selectivity.PCC_cat = PCC_cat;
Selectivity.DiagEntropyNorm = IS2;
Selectivity.DiagLRI = IS3;

Discrimination.NVocPerCat = NVocPerCat;
Discrimination.ProbaVocPerCat = ProbaVocPerCat;
Discrimination.NVoc = NVoc(1:ii);
Discrimination.VocCatNames = VocCatNames;
Discrimination.NbTrial_voc = NbTrial_voc;
MeanSR.MeanSR=MeanSR(1:ii);
MeanSR.SpikeRate_Cat = SpikeRate_Cat;
MeanSR.SpikeRate_Voc = SpikeRate_Voc(1:ii);
MeanSR.TDT_names = TDT_names(1:ii);

VoiceIndex.Observed = MI_perArea.MI_diag_uni_cat(:,1) ./ MI_perArea.MI_tot;
VoiceIndex.PvalueBGP = normpdf((MI_perArea.MI_diag_uni_cat(:,1) - MI_perArea.AVMI_diag_uni_cat_RandBG(:,1)) ./ MI_perArea.SDMI_diag_uni_cat_RandBG(:,1));

AgeIndex.Observed = MI_perArea.MI_diag_uni_cat(:,2) ./ MI_perArea.MI_tot;
AgeIndex.PvalueBGP = normpdf((MI_perArea.MI_diag_uni_cat(:,2) - MI_perArea.AVMI_diag_uni_cat_RandBG(:,2)) ./ MI_perArea.SDMI_diag_uni_cat_RandBG(:,2));

SexIndex.Observed = MI_perArea.MI_diag_uni_cat(:,3) ./ MI_perArea.MI_tot;
SexIndex.PvalueBGP = normpdf((MI_perArea.MI_diag_uni_cat(:,3) - MI_perArea.AVMI_diag_uni_cat_RandBG(:,3)) ./ MI_perArea.SDMI_diag_uni_cat_RandBG(:,3));

AgeSexIndex.Observed = MI_perArea.MI_diag_uni_cat(:,4) ./ MI_perArea.MI_tot;
AgeSexIndex.PvalueBGP = normpdf((MI_perArea.MI_diag_uni_cat(:,4) - MI_perArea.AVMI_diag_uni_cat_RandBG(:,4)) ./ MI_perArea.SDMI_diag_uni_cat_RandBG(:,4));
%% Save data
save(fullfile(resultsDirectory, 'VoiceAnalysis.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'VoiceAnalysis.mat'))



