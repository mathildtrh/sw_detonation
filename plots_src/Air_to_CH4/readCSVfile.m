function [Header, data] = readCSVfile(filename)

fid = fopen(filename, 'r');
Header = {fscanf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s', [1 10]) };


data = fscanf(fid, '%e,%e,%e,%e,%e,%e,%e,%e,%e,%e', [1 10]);

while ~feof(fid)
   data = [data; fscanf(fid, '%e,%e,%e,%e,%e,%e,%e,%e,%e,%e', [1 10]) ];
end

fclose(fid);


