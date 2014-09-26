function []=Cal_MIperarea_RandVoice(Randmatfile)

%% define some output variables
StimTypeCM=cell(1,4);
VoiceTypeSel=StimTypeCM;

AvMI_tot_Rand = nan(1,2);
SDMI_tot_Rand = nan(1,2);
AVMI_diag_uni_cat_Rand=cell(2,1);
AVMI_diag_uni_cat_Rand{1}=nan(1,length(StimTypeCM));
AVMI_diag_uni_cat_Rand{2}=AVMI_diag_uni_cat_Rand{1};
SDMI_diag_uni_cat_Rand=AVMI_diag_uni_cat_Rand;


%% Load file
RMAT=load(Randmatfile);

%% Cell aray containing the different classifications of the stims studied
StimTypeCM{1}=cell(11,1);
StimTypeCM{2}={'A' 'C' 'B'};
StimTypeCM{3}={'F' 'M' 'G'};
StimTypeCM{4}={'AF' 'AM' 'CF' 'CM' 'BG'};

% construct the cell array of name categories
VoiceTypeSel{1} = RMAT.VoiceTypeSel;
StimTypeCM{1}=unique(VoiceTypeSel{1});
IBG = find(strcmp(StimTypeCM{1}, 'BG'));
StimTypeCM{1} = [StimTypeCM{1}(1:(IBG-1)); StimTypeCM{1}((IBG+1):end); StimTypeCM{1}(IBG)];

% construct the cell array of age categories
VoiceTypeSel{2}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{2}{vv}=VoiceTypeSel{1}{vv}(1);
end

% construct the cell array of sex categories
VoiceTypeSel{3}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{3}{vv}=VoiceTypeSel{1}{vv}(2);
end

% construct the cell array of agesex categories
VoiceTypeSel{4}=cell(size(VoiceTypeSel{1}));
for vv=1:length(VoiceTypeSel{1})
    VoiceTypeSel{4}{vv}=VoiceTypeSel{1}{vv}(1:2);
end


rt = 2;%I set it to 2 on 08/11 2014 because we are not using values of random matrices without BG fixed
            while rt<=2
                if rt==1
                    LRM=length(RMAT.CM_IV_Rand);%number of random matrices
                elseif rt==2
                    LRM=length(RMAT.CM_IV_RandBG);%number of random matrices BG fixed
                end
                MI_tot_rand=zeros(LRM,1);
                MI_diag_uni_cat_rand=zeros(LRM,length(StimTypeCM));
                for rm =1:LRM
                    if rt==1
                        Rmat = RMAT.CM_IV_Rand{rm};
                    elseif rt==2
                        Rmat = RMAT.CM_IV_RandBG{rm};
                    end

                    % construct the cell array of indices for each random category
                    CAT_rand=cell(length(StimTypeCM),1);
                    for cl=1:length(StimTypeCM)
                        if rt==1 % random matrices
                            fprintf(1,'The code is not written here (lne515)!!!!');
                            return
                        elseif rt==2 % random BG Matrices
                            NstimTypeCM = length(StimTypeCM{cl});
                            Nb_VocPerCat = zeros(NstimTypeCM,1);
                            for vtR=1:NstimTypeCM
                                stR=StimTypeCM{cl}(vtR);
                                selectorR=strcmp(VoiceTypeSel{cl}, stR);
                                Nb_VocPerCat(vtR) = sum(selectorR);
                            end
                        end
                        NstimTypeRM = length(Nb_VocPerCat);
                        cat_rand = cell(NstimTypeRM,1);
                        nni = 0;
                        for cc = 1:NstimTypeRM
                                cat_rand{cc} = (nni+1):(Nb_VocPerCat(cc)+nni);
                                nni = nni+Nb_VocPerCat(cc);
                        end
                        CAT_rand{cl} = cat_rand;
                    end
                

                    [ mi_tot, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea_voice(Rmat, CAT_rand);
                    MI_tot_rand(rm)=mi_tot;
                    MI_diag_uni_cat_rand(rm,:)=mi_diag_uni_cat;    
                end
                AvMI_tot_Rand(rt) = mean(MI_tot_rand);
                SDMI_tot_Rand(rt) = std(MI_tot_rand);
                AVMI_diag_uni_cat_Rand{rt}(1,:)=mean(MI_diag_uni_cat_rand,1);
                SDMI_diag_uni_cat_Rand{rt}(1,:)=std(MI_diag_uni_cat_rand,1);
                if rt==1
                    MI_Calculations.MI_uni_diag_cat_Rand = MI_diag_uni_cat_rand;
                elseif rt==2
                    MI_Calculations.MI_uni_diag_cat_RandBG = MI_diag_uni_cat_rand;
                end
                rt = rt +1;
            end


%% Store values in a structure
MI_Calculations.AvMI_tot_Rand=AvMI_tot_Rand;
MI_Calculations.SDMI_tot_Rand=SDMI_tot_Rand;
MI_Calculations.CAT_rand=CAT_rand;
MI_Calculations.VoiceTypeSel=VoiceTypeSel;
MI_Calculations.StimTypeCM=StimTypeCM;

%% Add that new structure to the file
save(Randmatfile, 'MI_Calculations', '-append');
clear RMAT MI_Calculations CAT_rand mi_tot mi_diag_uni_cat mi_real_error_uni VoiceTypeSel StimTypeCM RMat
end
