function ign = ignition(temp_array)
% This function detects an ignition in the temperature evolution of a
% particle, by determining whether temperature becomes higher than the
% initial temperature

% if ign = true then ignition occurs
maximum = max(temp_array);
ign = (maximum > temp_array(1)+100);
end