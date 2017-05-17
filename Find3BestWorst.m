function [xanswerMin,xanswerMax] = Find3BestWorst(FvalStore,xanswerStore)
% Find the 3 best and 3 worst Fvals

xanswerMin = zeros(3,size(xanswerStore,2));
xanswerMax = zeros(3,size(xanswerStore,2));

% Search for the maximum and minimum values of F
tmp = sort(FvalStore(:),'ascend');
FMin1 = tmp(1);
locMin1 = find(FvalStore == FMin1);
xanswerMin(1,:) = xanswerStore(locMin1,:); %#ok<*FNDSB>
FMin2 = tmp(2);
locMin2 = find(FvalStore == FMin2);
xanswerMin(2,:) = xanswerStore(locMin2,:);
FMin3 = tmp(3);
locMin3 = find(FvalStore == FMin3);
xanswerMin(3,:) = xanswerStore(locMin3,:);
FMax1 = tmp(length(FvalStore(:)));
locMax1 = find(FvalStore == FMax1);
xanswerMax(1,:) = xanswerStore(locMax1(1),:);
FMax2 = tmp(length(FvalStore(:))-1);
locMax2 = find(FvalStore == FMax2);
xanswerMax(2,:) = xanswerStore(locMax2(1),:);
FMax3 = tmp(length(FvalStore(:))-2);
locMax3 = find(FvalStore == FMax3);
xanswerMax(3,:) = xanswerStore(locMax3(1),:);


end