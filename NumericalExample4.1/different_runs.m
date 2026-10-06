clear; clc; close all

r_start = 10;
r_step  = 10;
r_max   = 100;

r = r_start;

while r <= r_max

    fprintf('Running script with r = %d\n', r);
    scriptFile = 'dOpInf_Results.m';

    % Read original script
    txt = fileread(scriptFile);

    % Replace the line containing "r = ..."
    txt_mod = regexprep(txt, ...
        '(?m)^\s*r\s*=\s*\d+\s*;', ...
        sprintf('r = %d;', r));

    % Write modified script
    fid = fopen(scriptFile, 'w');
    fwrite(fid, txt_mod);
    fclose(fid);

    % Save current r for next iteration
    save('r_vals.mat','r','r_step','r_max');

    % Run the script (clear all inside is fine)
    run(scriptFile);
    
    load('r_vals.mat');
    % Increment r
    r = r + r_step;
end
