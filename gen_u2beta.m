% Script to generate values beta = beta(u).
% These values used in CAP to interpolate values of beta.
%

n = 1e3;    % NOTE this number needs to match NBETA in cap.h

u = linspace(0, 3*pi/4, n)';
beta = u2beta(u);

fid = fopen('generated_u2beta.c','w');

fprintf(fid, '// data generated using matlab script gen_u2beta.m\n\n');
fprintf(fid, 'GENU2BETA generated_u2beta[NBETA] = \n{\n');
for i = 1 : (n-1)
    fprintf(fid, '\t{%17.15f, %17.15f},\n', u(i), beta(i));
end
fprintf(fid, '\t{%17.15f, %17.15f}\n', u(i+1), beta(i+1));
fprintf(fid, '};\n');

fclose(fid);

