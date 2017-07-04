function [xanswerMin,xanswerMax] = Find3BestWorst(FvalStoreAll,xanswerStore, varargin)
% Find the  best Fvals
% no worst params calculation done yet.

if nargin>2
    nanswers = varargin{1};
else
    nanswers = 3; 
end

nDim = length(size(FvalStoreAll));
if nargin>3
    if varargin{2}
    nIndModels = size(FvalStoreAll, nDim);
    end
    ind_dim = nDim;
else
    nIndModels = 1;
    ind_dim = nDim+1;
end



xanswerMin = zeros(nanswers,size(xanswerStore,2), nIndModels );
xanswerMax = zeros(nanswers,size(xanswerStore,2), nIndModels );

% Search for the maximum and minimum values of F
for jj=1:nDim
    indx{jj}=1:size(FvalStoreAll,jj);
end
for model=1:nIndModels
    indx{ind_dim} = model;
    FvalStore = FvalStoreAll( indx{:} );
    tmp = sort(FvalStore(:),'ascend');
    for jj=1:nanswers
    FMin1 = tmp(jj);
    xanswerMin(jj,:, model) = xanswerStore(FvalStore == FMin1,:, model);
    end
% FMin2 = tmp(2);
% locMin2 = find(FvalStore == FMin2);
% xanswerMin(2,:) = xanswerStore(locMin2,:);
% FMin3 = tmp(3);
% locMin3 = find(FvalStore == FMin3);
% xanswerMin(3,:) = xanswerStore(locMin3,:);
% FMax1 = tmp(length(FvalStore(:)));
% locMax1 = find(FvalStore == FMax1);
% xanswerMax(1,:) = xanswerStore(locMax1(1),:);
% FMax2 = tmp(length(FvalStore(:))-1);
% locMax2 = find(FvalStore == FMax2);
% xanswerMax(2,:) = xanswerStore(locMax2(1),:);
% FMax3 = tmp(length(FvalStore(:))-2);
% locMax3 = find(FvalStore == FMax3);
% xanswerMax(3,:) = xanswerStore(locMax3(1),:);
end

xanswerMin = squeeze(xanswerMin);
xanswerMax = squeeze(xanswerMax);
end
