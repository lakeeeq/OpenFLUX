function CM = stepCV(CV,columns)
%author Lake-Ee Quek, AIBN
%to construction correction matrix (CM) from correction vector (CV)
%function CM = stepCV(CV,columns)
%specify CV: correction vector
%specify columns: number of columns for correction matrix (CM)

%convert CV to column vector if CV is a row vector
if size(CV,2) > size(CV,1)
    CV = CV';
end

CM=[];
for i = 1:columns;
    CM(i:i+length(CV)-1,i) = CV;
end
