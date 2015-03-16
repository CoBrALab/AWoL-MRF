%clearvars -except loaded_data; close all; clc; 
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using MRF model with simple Spanning Tree (no iterations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nameing:
%I = intensity image
%manL = manual labels
%maj_L = majority labels
%PD_p = probability distribution of HC labels
%fuse_L2 = segmentations produced by AWoL MRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters to set
% Input
input_file='/data/chamal/projects/nikhil/baby_data/folds/fold2/BABY_DATASET_Part1.mat';
input_struct='BABY_out';

% Output file prefix (L or R is added depending on HC_side)
output_file_pre='/data/chamal/projects/nikhil/baby_data/folds/fold2/Complete_Perf_BABY_fold2_test1';

foldID='fold2';
fold_csv_file='/data/chamal/projects/nikhil/baby_data/folds/fold2/fold2_data_BABY.csv';
% Mixing proportion / seed nbr count
deg = 10;
% patch size
patch = 5;
% Confidence thresholds
T_L1_orig = .6;
T_L0_orig = .8;


% MRF_beta (homogeneity for MRF)
mrf_beta = -.2;

% ATL types (2: even, 1:odd)
ATL_type = 0;
if ATL_type == 1
    output_file = strcat(output_file_pre,'_odd');
    maxATL = 5;
elseif ATL_type == 2
    output_file = strcat(output_file_pre,'_even');
    maxATL = 4;
else % comp datasets
    output_file = strcat(output_file_pre,'_comp');
    maxATL = 9;
end

% Max number of atlas and templates (indices not the actual count)
maxTMPL= 5;

% Hippocampus side (2: left, 1:right)
HC_side = 1;
if HC_side == 1
    output_file = strcat(output_file,'_Right_',date);
    fold_csv_file = strcat(fold_csv_file,'_Right_',date);
else
    output_file = strcat(output_file,'_Left_',date);
    fold_csv_file = strcat(fold_csv_file,'_Left_',date);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_file

%loaded_data = load(input_file, input_struct);
StructField = fieldnames(loaded_data);
input_data = loaded_data.(StructField{1});

maxSubjects = length(input_data);
Perf = zeros(maxSubjects,6);

%label_stats[grey, +, -, unreachable by Seeds(resolved by maj_vote)]
label_stats= zeros(maxSubjects,6);
nf=1;
bad_import = zeros(maxATL,maxTMPL,maxSubjects);
for ATL=1:1:maxATL
    for TMPL = 1:1:maxTMPL
        for sub=1:maxSubjects
            MaxL = ATL*(4*TMPL-1); %TMPL:3,7,9,11,15,19
            Dice = zeros(3,1);
            
            fp = zeros(2,1);
            fn = zeros(2,1);
            
            if HC_side == 1
                I = input_data(sub).img_R;
                manL = input_data(sub).manL_R;
                %Maj vote
                maj_L = input_data(sub).majL_R(:,:,:,ATL,TMPL);
                PDp = input_data(sub).PDp_R(:,:,:,ATL,TMPL);
                %PDn = MaxL - PDp; %input_data(sub).PDn(:,:,:,cadL);
            else
                I = input_data(sub).img_L;
                manL = input_data(sub).manL_L;
                %Maj vote
                maj_L = input_data(sub).majL_L(:,:,:,ATL,TMPL);
                PDp = input_data(sub).PDp_L(:,:,:,ATL,TMPL);
                %PDn = MaxL - PDp; %input_data(sub).PDn(:,:,:,cadL);
            end
            
            
            [Jaccard,D,rfp,rfn] = getDiceScore(manL,maj_L);
            Dice(1) = D;
            fp(1) = rfp;
            fn(1) = rfn;
            
            if(max(PDp(:))<= MaxL)
                
                Pr_update=ones(size(PDp,1),size(PDp,2),size(PDp,3)); %count of number times a voxel is assigned updated probabiltiy (used in EM)
                Pr = PDp./MaxL;
                
                I_size = size(I);
                
                fuse_label_map = [1 0];
                
                mix_ratio_met=1;
                mr_it = 0;
                T_L0 = T_L0_orig;
                T_L1 = T_L1_orig;
                
                global_mix_T = 100;
                
                %Make sure # of uncertain labels are less than know HC labels
                %This rarely happens but can be controlled by gloabl_mix_T
                while (mix_ratio_met > 0 && mr_it < 5)
                    
                    %assign partitioned labels: 0:Background, 1:HC, 0.5:
                    %uncertain
                    Base_L0 = PDp;
                    Base_L0(Base_L0 >=round(MaxL*T_L1)) = MaxL;
                    Base_L0(Base_L0 <= round(MaxL*(1-T_L0))) = -MaxL;
                    Base_L0(Base_L0 > -MaxL & Base_L0 < MaxL) = 0;
                    
                    Base_L0 = Base_L0./MaxL;
                    Base_L0 = (Base_L0 + 1)./2;
                    
                    LabelSizes0(ATL, TMPL, sub) = sum(Base_L0(:)==0);
                    LabelSizes1(ATL, TMPL, sub) = sum(Base_L0(:)==1);
                    LabelSizes05(ATL, TMPL, sub) = sum(Base_L0(:)==0.5);
                    round([ LabelSizes1(ATL, TMPL, sub), LabelSizes05(ATL, TMPL, sub),  ...
                        100*LabelSizes05(ATL, TMPL, sub)/ LabelSizes1(ATL, TMPL, sub)])
                    
                    % check if image has way too many uncertain labels
                    if  100*LabelSizes05(ATL, TMPL, sub)/ LabelSizes1(ATL, TMPL, sub) < global_mix_T
                        mix_ratio_met = 0;
                    else %Reduce the threshold to keep uncetain labels under limit (rarely the case)
                        if mod(mr_it,2) == 0
                            T_L0 = T_L0 - 0.05;
                        else
                            T_L0 = T_L0 - 0.05;
                        end
                    end
                    mr_it = mr_it + 1;
                end
                %Collect nbr and gradient information
                ps = 1;
                L=Base_L0;
                for i=1:I_size(1)
                    for j=1:I_size(2)
                        for k=1:I_size(3)
                            if i>ps && j>ps && k>ps && i< I_size(1)-ps+1 && j < I_size(2)-ps+1 && k < I_size(3)-ps+1
                                %grad = [I(i,j,k)-I(i-1,j,k), I(i+1,j,k)-I(i,j,k), I(i,j,k)-I(i,j-1,k),I(i,j+1,k)-I(i,j,k), I(i,j,k)-I(i,j,k-1),I(i,j,k+1)-I(i,j,k)];
                                DirGrad = [I(i-1,j,k)-I(i+1,j,k), I(i,j-1,k)-I(i,j+1,k), I(i,j,k-1)-I(i,j,k+1)];
                                Info_struct_array(i,j,k).DG = DirGrad;
                                grad = I(i-ps:i+ps,j-ps:j+ps,k-ps:k+ps)- I(i,j,k);
                                Info_struct_array(i,j,k).IG = abs(grad).^1;
                            else
                                Info_struct_array(i,j,k).IG = 1000*ones(3,3,3);
                                Info_struct_array(i,j,k).DG = zeros(1,3);
                            end
                        end
                    end
                end
                % Avoid hitting the edges..
                L(1:2,:,:) = maj_L(1:2,:,:);
                L(:,1:2,:) = maj_L(:,1:2,:);
                L(:,:,1:2) = maj_L(:,:,1:2);
                L(size(L,1)-1:size(L,1),:,:) = maj_L(size(L,1)-1:size(L,1),:,:);
                L(:, size(L,2)-1:size(L,2),:) = maj_L(:, size(L,2)-1:size(L,2),:);
                L(:, :, size(L,3)-1:size(L,3)) = maj_L(:, :, size(L,3)-1:size(L,3));
                
                NBR = getNBRInfo(L,deg);
                AP_p = NBR.AP_p(:);
                AP_n = NBR.AP_n(:);
                if max(NBR.AP_p(:)) > 1
                    Seeds = NBR.Seeds;
                    seedsize(sub)=size(Seeds,1)
                    
                    %Start the walk on tree
                    label_up_p=zeros(I_size(1),I_size(2),I_size(3));
                    label_up_n=zeros(I_size(1),I_size(2),I_size(3));
                    %fuse_L = maj_L;
                    fuse_L2 = maj_L;
                    
                    for s=1:min(size(Seeds,1),500)
                        new_L = L;
                        stop = 0;
                        p(1) = Seeds(s,1);
                        p(2) = Seeds(s,2);
                        p(3) = Seeds(s,3);
                        
                        prev_candidate1 = p;
                        
                        % Generate spanning tree for local patch
                        MST_Output = getMSTSequencePlus(L,Info_struct_array,I,p,patch);
                        
                        walk_seq = MST_Output.STP;
                        walk_patch_stats = MST_Output.patch_stats;
                        walk_seq_array(s).seq = walk_seq;
                        if length(walk_seq) > 0 && sum(walk_patch_stats) ~= 0
                            t=1;
                            Wsize(s,:) = size(walk_seq);
                            
                            %If we want to randomize walk seq
                            %[Ws,Wi] = sort(rand(Wsize(1),1));
                            while t <= Wsize(s,1) && stop~=1
                                label_wt = zeros(2,1);
                                label_Vwt = zeros(2,1);
                                labelV = zeros(2,3);
                                
                                %     Entire patch
                                % May want to scale by std(I(:))
                                localI = I(p(1)-ps:p(1)+ps,p(2)-ps:p(2)+ps,p(3)-ps:p(3)+ps);
                                i1=1;
                                i0=1;
                                beta=1;
                                
                                %MRF notes
                                % Not summing over entire patch because except
                                % of the node under considearation all other
                                % singleton potnetial remain the same while
                                % performing iterative updates. Also doubleton
                                % cliques consist of only the 6 nbrs of p
                                
                                %MRF singleton component
                                mrf_single(1) = log(sqrt(2*pi)*walk_patch_stats(2)) + ((I(p(1),p(2),p(3)) - walk_patch_stats(1))^2)./(2*(walk_patch_stats(2)^2));
                                mrf_single(2) = log(sqrt(2*pi)*walk_patch_stats(4)) + ((I(p(1),p(2),p(3)) - walk_patch_stats(3))^2)./(2*(walk_patch_stats(4)^2));
                                
                                %MRF doubleton component
                                doubleton_tmp = [fuse_L2(p(1)-1,p(2),p(3)), fuse_L2(p(1)+1,p(2),p(3)),...
                                    fuse_L2(p(1),p(2)-1,p(3)), fuse_L2(p(1),p(2)+1,p(3)),...
                                    fuse_L2(p(1),p(2),p(3)-1), fuse_L2(p(1),p(2),p(3)+1)];
                                
                                mrf_double(1) = mrf_beta*(sum(doubleton_tmp==1) - sum(doubleton_tmp==0));
                                mrf_double(2) = mrf_beta*(sum(doubleton_tmp==0) - sum(doubleton_tmp==1));
                                
                                cliques = [mrf_single;mrf_double];
                                %MRF total clique energy
                                mrf_energy = [mrf_single(1)+mrf_double(1), mrf_single(2)+mrf_double(2)];
                                [probMRFWt probMRF] = min(mrf_energy);
                                new_L(p(1),p(2),p(3)) = mod(probMRF,2);
                                
                                if new_L(p(1),p(2),p(3)) == 1
                                    label_up_p(p(1),p(2),p(3)) = label_up_p(p(1),p(2),p(3)) + 1;
                                else
                                    label_up_n(p(1),p(2),p(3)) = label_up_n(p(1),p(2),p(3)) + 1;
                                end
                                [fused_L_value,fused_L_ind] = max([label_up_p(p(1),p(2),p(3)),label_up_n(p(1),p(2),p(3))]);
                                fuse_L2(p(1),p(2),p(3)) =  fuse_label_map(fused_L_ind);
                                vts = [manL(p(1),p(2),p(3)), maj_L(p(1),p(2),p(3)),fuse_L2(p(1),p(2),p(3))];
                                t=t+1;
                                if t <= Wsize(s,1)
                                    Next_candidate = walk_seq(t,:);
                                    p = Next_candidate;
                                end
                            end
                        end
                    end
                    
                    [Jaccard,D,rfp,rfn] = getDiceScore(manL,fuse_L2);
                    Dice(2) = D;
                    fp(2) = rfp;
                    fn(2) = rfn;
                    
                else
                    fuse_L2=maj_L;
                    Dice(2) = Dice(1);
                    fp(2) = fp(1);
                    fn(2) = fn(1);
                end
            else
                % exception handling for bad imports
                bad_import(ATL, TMPL, sub) = 1;
                fuse_L2=maj_L;
                Dice(2) = Dice(1);
                fp(2) = fp(1);
                fn(2) = fn(1);
            end
            % Back to nifti conversion: to visualize / compare dice values etc.
            %nii_out = make_nii(fuse_L);
            %save_nii(nii_out, subjects(sub));
            %fuseEM = EM_fun(Pr(:),fuse_L(:));
            Dice
            AWoL_Array(ATL,TMPL,sub).fuse_L2 = fuse_L2;
            manVol(ATL, TMPL, sub) = sum(manL(:)==1);
            majVol(ATL, TMPL, sub) = sum(maj_L(:)==1);
            AwolVol(ATL, TMPL, sub) = sum(fuse_L2(:)==1);
            
            csv_str = {foldID, num2str(2*ATL),num2str(4*TMPL-1),input_data(sub).subject,num2str(HC_side),num2str(Dice(1)),num2str(Dice(2))};
            csv_save(fold_csv_file,csv_str);
            Perf(sub,:) = [Dice(1),Dice(2),fp(1),fp(2),fn(1),fn(2)];
            majDice(ATL, TMPL, sub) = Dice(1);
            AwolDice(ATL, TMPL, sub) = Dice(2);
            
            clear manL;
            clear maj_L;
            clear L_array;
            clear fuse_L;
            clear fuse_L2;
            clear L0;
            clear Base_L0;
            clear I;
            clear L;
            clear new_L;
            
            fprintf('%f ATL number ',ATL);
            fprintf(' %f percent subjects complete\n', 100*sub/maxSubjects)
        end
        cadPerfMean(ATL,TMPL,:) = mean(Perf);
        cadPerfMax(ATL,TMPL,:) = max(Perf);
        cadPerfStd(ATL,TMPL,:) = std(Perf);
        
        fprintf('%f ATL number ',ATL);
        fprintf(' %f TMPL number complete\n', TMPL);
        
    end
    fprintf('%f ATL number complete\n',ATL);
end
save(output_file,'AWoL_Array','manVol','majVol','AwolVol','majDice','AwolDice','-v7.3');