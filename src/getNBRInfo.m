function NBR = getNBRInfo(L,N)
%L:Label Matrix, N=degree of the vertex (number of neighbors)
n=1;
totalNBR=26;
NBR.AP_p=[];
NBR.AP_n=[];
NBR.Seeds=[];
%Lm = L(2:L_size(1)-1,2:L_size(2)-1,2:L_size(3)-1);
[greyPx, greyPy, greyPz] =ind2sub(size(L),find(L==0.5));
greyP = [greyPx greyPy greyPz];
sz = size(L); %# size of the cube along each dimension
seed_indices=[];
for m=1:length(greyPx)
    coords = greyP(m,:);
    %# generate increments to get the 26 neighbors around a 3D point
    [X, Y, Z] = ndgrid([-1 0 1], [-1 0 1], [-1 0 1]);
    nb = [X(:) Y(:) Z(:)];
    nb(ismember(nb,[0 0 0],'rows'),:) = [];  %# remove the row [0 0 0]
    
    %# for each 3D point, compute its neighbors
    szMod = repmat(sz, [size(nb,1) 1]);
 
    cc = bsxfun(@plus, nb, coords); %# find 26 neighbors of coords(i,:)
    for g=1:length(cc)
       nbrL(g) =  L(cc(g,1),cc(g,2),cc(g,3));
    end
    NBR.AP_p(greyP(m,1),greyP(m,2),greyP(m,3)) = sum(nbrL==1);
    NBR.AP_n(greyP(m,1),greyP(m,2),greyP(m,3)) = sum(nbrL==0);
    
   %if NBR.AP_p(greyP(m,1),greyP(m,2),greyP(m,3)) + NBR.AP_n(greyP(m,1),greyP(m,2),greyP(m,3)) >= N
   if NBR.AP_p(greyP(m,1),greyP(m,2),greyP(m,3)) >= N
       NBR.Seeds(n,:) = greyP(m,:);
       seed_indices(n) = m;
       n = n+1;
   end   
end

