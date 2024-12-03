function [data, headers, everything] = FNreadcsv (fname,rowsANDcols)
% READING CSV FILEs on all systems
% fname is name of the CSV file (must be comma seperated, not point-vergule).
% rowsANDcols used for indicate if row and column names exist
% default is [1,1] meaning both row and column names exist
% a zero indicates nothing e.g. [0,0] means no row names or column names,
% [0,1] means just column names and [1,0] means just row names.


% set default
if nargin==1
    rows=1 ; cols=1;  
else
    rows=rowsANDcols(1); 
    cols=rowsANDcols(2); 
end
%fname='pricesDAILY.csv'; rows=1 ; cols=1;  % testing ...

% open file
fid = fopen(fname,'rt');
tmp = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
tmp = tmp{1};

%# split and concatenate the rest
result = regexp(tmp,',','split');
result = cat(1,result{:});

% convert missing values to NaN
result1=result;
result(strcmp(result,''))={NaN};
result1(strcmp(result1,''))={'-999.999'};
%result1(strcmp(result1,''))={'NaN'};  % results in bugs related to str2num performance !

% get row and column names
data1=str2num(char(result1((cols+1):end,(rows+1):end)));
data1(data1==-999.999)=NaN;

%data1=str2double(char(result((rows+1):end,(cols+1):end)));
if isempty(data1); 
    %error('ERROR: check data is numeric and headers specified correctly (non-numeric conversion error)'),return ;
    warning('WARNING: non-numeric data OR incorrect row/col specification')
    data = result(cols+1:end,rows+1:end);
    dataCELL=data;
      
    nums=cellfun(@str2double,dataCELL);
    dataCELL(~isnan(nums))=num2cell(nums(~isnan(nums)));
    dataCELL(strcmp(dataCELL,''))={NaN};
    
else
    data=reshape(data1,size(result)-[cols , rows]);
    dataCELL=num2cell(data);
end;


if rows==1
    rownames=result((cols+1):end,1);
    everything= [rownames,dataCELL] ;
else
    rownames={};
    everything=dataCELL;
end

if cols==1; 
    colnames=result(1,1:end);  
    everything=[ [colnames] ; everything ] ;
else
    colnames={};
end;

headers.rownames=rownames;
headers.colnames=colnames;

end