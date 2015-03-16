function Output = getMSTSequencePlus(L,I,Int,p,patch)
%use I.NL (neighbor labels info) and move in 1-nbr-patches (step size of 3)
%using only +x coordinate nbr labels for now
%patch size = (2*patch_length + 1)^2
sizeL = size(L);
truncXMax = 0;
truncXMin = 0;
truncYMax = 0;
truncYMin = 0;
truncZMax = 0;
truncZMin = 0;
% row == Y axis col == X Axis

%add another dim here
if p(1)+patch > sizeL(1)
    truncXMax = p(1) + patch - sizeL(1);
end
if p(1) - patch < 0
    truncXMin = - p(1) + patch;
end
if p(2)+patch > sizeL(2)
    truncYMax = p(2) + patch - sizeL(2);
end
if p(2) - patch < 0
    truncYMin = - p(2) + patch;
end
if p(3)+patch > sizeL(3)
    truncZMax = p(3) + patch - sizeL(3);
end
if p(3) - patch < 0
    truncZMin = - p(3) + patch;
end

patch = patch - max([truncXMax, truncXMin, truncYMax, truncYMin, truncZMax, truncZMin]) - 2;
if patch < 0
    STP = 0;
else
    patch_length = 2*patch  + 1;
    
    %find candidnates
    candidates(1,:) = p;
    n=1;
    n1=1;
    n0=1;
    for i=1:patch_length
        for j=1:patch_length
            for k=1:patch_length
                shiftX = p(1) - patch - 1 + i;
                shiftY = p(2) - patch - 1 + j;
                shiftZ = p(3) - patch - 1 + k;
                
                %NL(2) is the +X nbr
                %                 % Need to see why this happens...
                %                 if shiftX <= sizeL(1) && shiftY<=sizeL(2) && shiftZ <= sizeL(3)
                if L(shiftX,shiftY, shiftZ) == 0.5 && norm(p-[shiftX,shiftY,shiftZ]) > 0
                    n = n+1;
                    candidates(n,:) = [shiftX,shiftY,shiftZ];
                elseif L(shiftX,shiftY, shiftZ) == 1 && norm(p-[shiftX,shiftY,shiftZ]) > 0
                    L1_Points(n1) = Int(shiftX,shiftY,shiftZ);
                    n1=n1+1;
                elseif L(shiftX,shiftY, shiftZ) == 0 && norm(p-[shiftX,shiftY,shiftZ]) > 0
                    L0_Points(n0) = Int(shiftX,shiftY,shiftZ);
                    n0=n0+1;
                end
                %                end
            end
        end
    end
    candidates;
    %assign weights
    if n>1
        dist_wt = 100; %wt to avoid jumps over nbrs
        weightM = zeros(n,n);
        for wx = 1:n
            for wy = wx:n
                if wx == wy
                    weightM(wx,wy) = 0;
                elseif norm(candidates(wx,:,:)-candidates(wy,:,:)) > 1
                    weightM(wx,wy) = dist_wt*norm(candidates(wx,:,:)-candidates(wy,:,:));
                    weightM(wy,wx) = weightM(wx,wy);
                else
                    currentP = candidates(wx,:,:);
                    nbrP_shift = candidates(wy,:,:) - candidates(wx,:,:);
                    nbrP_shiftX = nbrP_shift(1);
                    nbrP_shiftY = nbrP_shift(2);
                    nbrP_shiftZ = nbrP_shift(3);
                    weightM(wx,wy) = abs(I(currentP(1),currentP(2),currentP(3)).IG(nbrP_shiftX+2,nbrP_shiftY+2,nbrP_shiftZ+2));
                    weightM(wy,wx) = weightM(wx,wy);
                end
            end
        end
        weightM;
        MST = prim(weightM);
        M = [MST(:,1), MST(:,2)];
        n=1;
        seq(n) = M(1,1);
        STP(n,:) = candidates(seq(n),:);
        sizeM = size(M);
        for j = 1:sizeM(1)
            for k = 1:2
                if sum(seq(1:n) == M(j,k)) == 0
                    n = n+1;
                    seq(n) = M(j,k);
                    STP(n,:) = candidates(seq(n),:);
                end
            end
        end
    elseif n == 1
        STP = p;
    else
        STP = 0;
    end
end
Output.STP = STP;
if n0 == 1 || n1 == 1
    Output.patch_stats = 0;
else
    Output.patch_stats = [mean(L1_Points),std(L1_Points),mean(L0_Points),std(L0_Points)];
end