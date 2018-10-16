function corr_vect = cVectGen2(iDist, ele, MIDlength)
%author Lake-Ee Quek, AIBN
%function corr_vect = cVectGen2(iDist, ele, MIDlength)
%generate correction vector for a given natural isotope distribution vector and number of atomic elements
%iDist = isotopomer distribution row vector of the atom (from from m0, m1,m2,...)
%ele = number of atoms
%MIDlength = MID vector length (i.e, mass increments, include m0)
% clearvars
% clc
% iDist = [0.9 0.06 0.03 0.01];
% ele = 7;
% MIDlength = 5;
%doesnt add up to 1 unless MIDlength = ele+(iDit-1)*2

corr_vect = zeros(MIDlength,1);

noIsotope = length(iDist);
if ele < MIDlength
    sampleTable = zeros(1,ele);
else
    sampleTable = zeros(1,MIDlength);
end

for i = 1:noIsotope-1
    noSamples = floor((MIDlength-1)/i);
    if  noSamples > ele
        noSamples = ele;
    end
    sampleTable = [sampleTable ones(1,noSamples)*i];
end

if ele>MIDlength
    combi = combnk(sampleTable,MIDlength);
    combi = [zeros(size(combi,1),ele-MIDlength) combi];
else
    combi = combnk(sampleTable,ele);
end

combi_sum = sum(combi,2);

for i = 1:MIDlength
    hit = combi_sum == i-1;
    combi_hit_unique = unique(combi(hit,:),'rows');
    
    probabilitySum = 0;
    for j = 1:size(combi_hit_unique,1)
        %         totalPerms = unique(perms(combi_hit_unique(j,:)),'rows'); %for actual comparison
        
        %calculate analytically
        [isotope_unique,~,index_B] = unique(combi_hit_unique(j,:));
        countOcc = zeros(1,numel(isotope_unique));
        eleNow = ele;
        noPerms = 1;
        probMultiply = 1;
        for k = 1:numel(isotope_unique)
            noOcc = sum(index_B==k);
            countOcc(k) = noOcc;
            probMultiply = probMultiply*iDist(isotope_unique(k)+1)^noOcc;
            if eleNow>noOcc
                noPerms = noPerms*nchoosek(eleNow,noOcc);
            end
            eleNow = eleNow - noOcc;
        end
        %         disp([size(totalPerms,1) noPerms]);%validate occurance
        probabilitySum = probabilitySum + noPerms*probMultiply;
    end
    corr_vect(i) = probabilitySum;
end
end