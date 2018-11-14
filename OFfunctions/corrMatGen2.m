function CM = corrMatGen2(output_row, input_size, molform)
%author Lake-Ee Quek, AIBN
%function CM = corrMatGen(output_size, input_size, molform)
%to generate correction matrix from a given molecular formula
%example: outputMDV = CM * inputMDV
%output_size: length of output MDV (truncated length)
%input_size: lenght of input MDV (carbon backbone + 1)
%molform: molecular formula without carbon backbone
output_size=max(output_row);
if nargin == 2 || isempty(molform)==1
    CM = eye(output_size,input_size);
    return
end

try
    natDist;
catch
    C_dist=[];
    H_dist=[];
    N_dist=[];
    O_dist=[];
    Si_dist=[];
    S_dist=[];
    P_dist=[];
end
%natural isotope distribution ref: van winden (2002), BnB
if isempty(C_dist)
    C_dist=[0.9893 0.0107];
end
if isempty(H_dist)
    H_dist=[0.999885 0.000115];
end
if isempty(N_dist)
    N_dist=[0.99632 0.00368];
    
end
if isempty(O_dist)
    O_dist=[0.99757 0.00038 0.00205];
end
if isempty(Si_dist)
    Si_dist=[0.922297 0.046832 0.030872];
    
end
if isempty(S_dist)
    S_dist=[0.9493 0.0076 0.0429 0.0002];
end
if isempty(P_dist)
    P_dist=[1 0];
end

pos = isletter(molform);
atom={};
atom{1} ='';
for i = 1:length(pos)
    if pos(i)==1
        atom{end} = strcat(atom{end},molform(i));
    elseif isempty(atom{end}) == 0
        atom{end+1} = '';
    end
end
atom(end)=[];


coeff={};
coeff{1}='';
for i = 1:length(pos)
    if pos(i)==0
        coeff{end} = strcat(coeff{end},molform(i));
    elseif isempty(coeff{end}) == 0
        coeff{end+1} = '';
    end
end
if length(coeff) ~= length(atom)
    %     fprintf('error occured in %s\n',m);
    fprintf('number of atom label (%1.0f labels) did not match coefficients (%1.0f coefficients)\n', length(atom), length(coeff));
    fprintf('press Ctrl C to terminate');
    pause(inf)
end
for i = 1:length(atom)
    atom{i} = strcat(atom{i}, '_dist,', coeff{i});
end
CM=1;
for i = 1:length(atom)
    CM = cauchy(CM,eval(['cVectGen2(' atom{i} ',output_size)'])');
end
CM = stepCV(CM,input_size);
CM = CM(1:output_size,:);
CM = CM(output_row,:);


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
        for i = 1:columns
            CM(i:i+length(CV)-1,i) = CV;
        end
        
    end
end