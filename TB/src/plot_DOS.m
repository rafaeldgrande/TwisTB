% plot_DOS.m  -- called as a script from main.m when task == 3
% Computes the DOS via Gaussian smearing over the uniform k-mesh and
% saves a .dat file plus a .png figure.
%
% Workspace variables used (set by main.m before this script is called):
%   tb_bands   (neigs x knum_tot)  eigenvalues in eV
%   noccs_k    (knum_tot x 1)      number of occupied bands per k-point
%   neigs, knum_tot                dimensions
%   outfname, dirname              output file base name and directory

fprintf('--> Computing DOS ... ')

sigma = 0.02;    % Gaussian broadening (eV)
ne    = 2000;    % energy grid points

emin = min(tb_bands(:)) - 0.5;
emax = max(tb_bands(:)) + 0.5;
egrid = linspace(emin, emax, ne);

dos = zeros(1, ne);
pref = 1.0 / (sigma * sqrt(2*pi));
for ib = 1 : neigs
    for ik = 1 : knum_tot
        dos = dos + pref * exp(-0.5 * ((egrid - tb_bands(ib,ik)) / sigma).^2);
    end
end
dos = dos / knum_tot;

% Fermi level: VBM from the last occupied band at each k, take maximum
n_occ = round(median(noccs_k));
vbm   = max(tb_bands(n_occ, :));
egrid_shifted = egrid - vbm;

fprintf('done\n')

% Save .dat
fprintf('--> Saving DOS in %s_DOS.dat ... ', outfname)
fname_dat = fullfile(dirname, [outfname, '_DOS.dat']);
fid = fopen(fname_dat, 'w');
fprintf(fid, '# E-Ef(eV)   DOS(states/eV/cell)\n');
fprintf(fid, '%.6f  %.8f\n', [egrid_shifted; dos]);
fclose(fid);
fprintf('done\n')

% Plot
fprintf('--> Plotting DOS ... ')
fig = figure('visible','off');
plot(egrid_shifted, dos, 'r-', 'LineWidth', 1.5);
xline(0, 'k--', 'LineWidth', 0.8);
xlabel('E - E_{VBM} (eV)');
ylabel('DOS (states/eV/cell)');
title(outfname, 'Interpreter', 'none');
xlim([-4, 4]);
grid on;
fname_png = fullfile(dirname, [outfname, '_DOS.png']);
saveas(fig, fname_png);
close(fig);
fprintf('done\n')
fprintf('--> Saved: %s\n', fname_dat)
fprintf('--> Saved: %s\n', fname_png)
