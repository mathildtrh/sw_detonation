function [Header, data] = hfile(filename)

fid = fopen(filename, 'r');

nHeaderRows = 2;
nColumns = 15;

nHeaderElements = nHeaderRows*nColumns;

Header = {fscanf(fid, '%s', [1 1]) };

for i = 2:nHeaderElements
   Header = [Header fscanf(fid, '%s', [1 1]) ]; 
    i = i+1;
end


data = fscanf(fid, '%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e', [1 15]);

while ~feof(fid)
   data = [data; fscanf(fid, '%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e', [1 15]) ];
end

fclose(fid);


