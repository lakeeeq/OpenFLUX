function Nout = bSplineMat_lite(knotSeq,tData,order)
%%%this lite script takes a single time point (not a vector)
noKnots = numel(knotSeq);
if tData == knotSeq(end)
    tData = tData - 1e-8;
end

%general requirements
ti = knotSeq';
dataTdiff = tData - ti; %left or term 1 form 

%contruct Ni,0
term1 = dataTdiff(1:end-1);
term2 = -dataTdiff(2:end);
Ni_0 = term1>=0 & term2>0;
Nin = +Ni_0;
% N_store{end+1} = Nin;
for curOrder = 1:order-1
    knotDiff = knotSeq(curOrder+1:end)-knotSeq(1:end-curOrder);
    knotDiff(knotDiff==0) = inf;
    knotDiff = 1./knotDiff;
    term1_denom = diag(knotDiff(1:end-1));
    term2_denom = diag(knotDiff(2:end));
    term1 = dataTdiff(1:end-1-curOrder);
    term2 = -dataTdiff(2+curOrder:end);
    Nin = term1_denom*term1.*Nin(1:end-1) + term2_denom*term2.*Nin(2:end);
end
Nout = Nin;