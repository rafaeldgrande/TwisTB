%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to find intralayer and interlayer
% neighbors for each orbital
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% orbitals - orbitals class
% norbs - Number of orbitals
% cell - lattice vectors of moire cell (ang)
% ncell1,ncell2 - number of repeating cells in
%                 a1 and a2 directions, respectively
% aXM -  minimum distance between a metal atom centre
%        projected chalcogen atom centre onto metal plane
% nspin - npsins = 2 with SOC, nspin = 1 otherwise
%
% Outputs:
% orbitals - orbitals class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale
% dz2-pz neighbor list added by Kemal Atalar
% Parallelised (parfor) by R. Del Grande
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orbitals] = opt_find_neigh_3rdnn(orbitals,mcell,ncell1,ncell2,...
     bondlength,aXX2,nspin,type1,type2,type3,type4,type5,type6,...
     ia1,ia2,Rot,multilayer,nlayer,interlayer_int,Interpd,flipped)

% Make this an input later, X-M interlayer coupling cutoff
aXM = 6;

strain = 0.1;
tol = strain * bondlength;

% ── Pre-extract orbital properties into plain arrays for parfor broadcast ──
% orbitals is a MATLAB object array; broadcasting it to parfor workers is
% expensive and incompatible with some class implementations.  We read all
% fields we need once here and pass plain numeric arrays into the loops.
tot_norbs_all = size(orbitals, 2);
orb_centre_all  = zeros(tot_norbs_all, 3);
orb_Pcentre_all = zeros(tot_norbs_all, 3);
orb_M2sym_all   = zeros(tot_norbs_all, 1);
orb_M1sym_all   = zeros(tot_norbs_all, 1);
orb_rel_idx_all = zeros(tot_norbs_all, 1);
for iorb = 1:tot_norbs_all
    orb_centre_all(iorb, :)  = orbitals(iorb).Centre;
    pc = orbitals(iorb).Pcentre;
    if ~isempty(pc)
        orb_Pcentre_all(iorb, :) = pc;
    end
    orb_M2sym_all(iorb)      = orbitals(iorb).M2sym;
    orb_M1sym_all(iorb)      = orbitals(iorb).M1sym;
    orb_rel_idx_all(iorb)    = orbitals(iorb).Rel_index;
end

% ── Main neighbor loop ────────────────────────────────────────────────────
for spin = -nspin + 1 : 2 : nspin - 1
    for layer = 1 : nlayer
        if(flipped(layer)==0)
            a1 = ia1*Rot{layer}';
            a2 = ia2*Rot{layer}';
        elseif(flipped(layer)==1)
            a1 = ia1*Rot{layer}';
            a2 = [ia2(1),-ia2(2),ia2(3)]*Rot{layer}';
        end

        bl        = bondlength(layer);
        tol_l     = tol(layer);
        sqrt3_bl  = sqrt(3)*bl;
        sqrt3_tol = sqrt(3)*tol_l;
        two_bl    = 2.0*bl;
        two_tol   = 2.0*tol_l;

        for isym = [-1,1]
            t_isym = tic; fprintf('    layer %i/%i, isym %+i ... ', layer, nlayer, isym)
            orb1   = findobj(orbitals,'Layer',layer,'Spin',spin,'M1sym',isym);
            norbs  = size(orb1,1);
            coor   = [reshape([orb1(:).Centre],[3,norbs])'];
            hindex = [[orb1(:).Ham_index]'];

            supercell   = zeros(ncell1^2*norbs, 2);
            super_hidx  = zeros(ncell1^2*norbs, 3);
            ir = 0;
            for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                    ir = ir + 1;
                    supercell((ir-1)*norbs+1:ir*norbs, 1:2) = ...
                        coor(:,1:2) + n1*mcell(1,1:2) + n2*mcell(2,1:2);
                    super_hidx((ir-1)*norbs+1:ir*norbs, :) = ...
                        [hindex(:), n1*ones(norbs,1), n2*ones(norbs,1)];
                end
            end
            dist = sqrt( bsxfun(@plus, sum(coor(:,1:2).^2,2), sum(supercell.^2,2)') ...
                         - 2*(coor(:,1:2)*supercell') );

            % Pre-allocate cell arrays for parfor outputs
            nn1_list = cell(norbs,1); nn1_cen = cell(norbs,1); nn1_hop = cell(norbs,1);
            nn2_list = cell(norbs,1); nn2_cen = cell(norbs,1); nn2_hop = cell(norbs,1);
            nn3_list = cell(norbs,1); nn3_cen = cell(norbs,1); nn3_hop = cell(norbs,1);

            % Broadcast scalars for parfor
            isym_bc  = isym;
            layer_bc = layer;

            parfor ii = 1:norbs
                iorb     = hindex(ii);
                my_cen   = orb_centre_all(iorb, :);
                ir_val   = orb_rel_idx_all(iorb);
                my_M2sym = orb_M2sym_all(iorb);

                % ── NN1 ──────────────────────────────────────────────────
                dl = find(abs(dist(ii,:) - bl) <= tol_l);
                j1  = super_hidx(dl, 1);
                n1v = super_hidx(dl, 2);  n2v = super_hidx(dl, 3);
                c1  = orb_centre_all(j1, :) + n1v.*mcell(1,:) + n2v.*mcell(2,:);
                if length(j1) > 9
                    error('d/p orbital: too many NN1 neighbours (%d)', length(j1))
                end
                nn1_list{ii} = j1';  % row vector: build_H_3rdnn uses size(list,2)
                nn1_cen{ii}  = c1;
                nn1_hop{ii}  = assign_hoppings(a1, a2, c1 - my_cen, 1, layer_bc, ...
                                   my_M2sym * orb_M2sym_all(j1)', ...
                                   ir_val, orb_rel_idx_all(j1)', ...
                                   type1,type2,type3,type4,type5,type6);

                % ── NN2 ──────────────────────────────────────────────────
                dl = find(abs(dist(ii,:) - sqrt3_bl) <= sqrt3_tol);
                j2  = super_hidx(dl, 1);
                n1v = super_hidx(dl, 2);  n2v = super_hidx(dl, 3);
                c2  = orb_centre_all(j2, :) + n1v.*mcell(1,:) + n2v.*mcell(2,:);
                if length(j2) > 18
                    error('d/p orbital: too many NN2 neighbours (%d)', length(j2))
                end
                nn2_list{ii} = j2';
                nn2_cen{ii}  = c2;
                nn2_hop{ii}  = assign_hoppings(a1, a2, c2 - my_cen, 2, layer_bc, ...
                                   my_M2sym * orb_M2sym_all(j2)', ...
                                   ir_val, orb_rel_idx_all(j2)', ...
                                   type1,type2,type3,type4,type5,type6);

                % ── NN3 ──────────────────────────────────────────────────
                dl = find(abs(dist(ii,:) - two_bl) <= two_tol);
                j3  = super_hidx(dl, 1);
                n1v = super_hidx(dl, 2);  n2v = super_hidx(dl, 3);
                c3  = orb_centre_all(j3, :) + n1v.*mcell(1,:) + n2v.*mcell(2,:);
                if length(j3) > 9
                    error('d/p orbital: too many NN3 neighbours (%d)', length(j3))
                end
                nn3_list{ii} = j3';
                nn3_cen{ii}  = c3;
                if isym_bc > 0
                    % Filter neighbours by M1sym (replaces findobj in original)
                    valid3 = (orb_M1sym_all(j3) == isym_bc);
                    jr3    = orb_rel_idx_all(j3(valid3))';
                    nn3_hop{ii} = assign_hoppings(a1, a2, c3 - my_cen, 3, layer_bc, ...
                                      my_M2sym * orb_M2sym_all(j3)', ...
                                      ir_val, jr3, ...
                                      type1,type2,type3,type4,type5,type6);
                end
            end  % parfor ii

            % ── Assign results back to orbitals (serial) ─────────────────
            for ii = 1:norbs
                iorb = hindex(ii);
                orbitals(iorb).Intra_nn1_list     = nn1_list{ii};
                orbitals(iorb).Intra_nn1_centre   = nn1_cen{ii};
                orbitals(iorb).Intra_nn1_hoppings = nn1_hop{ii};
                orbitals(iorb).Intra_nn2_list     = nn2_list{ii};
                orbitals(iorb).Intra_nn2_centre   = nn2_cen{ii};
                orbitals(iorb).Intra_nn2_hoppings = nn2_hop{ii};
                orbitals(iorb).Intra_nn3_list     = nn3_list{ii};
                orbitals(iorb).Intra_nn3_centre   = nn3_cen{ii};
                if ~isempty(nn3_hop{ii})
                    orbitals(iorb).Intra_nn3_hoppings = nn3_hop{ii};
                end
            end
            fprintf('done [%.1fs]\n', toc(t_isym))

        end  % for isym
    end  % for layer

    clear dist coor orb1 supercell norbs super_hidx

    % ── Interlayer neighbors ──────────────────────────────────────────────
    if(multilayer && interlayer_int)
        for layer = 1 : nlayer - 1

            % ── X-X (p-p) interlayer coupling ────────────────────────────
            t_inter = tic; fprintf('    interlayer X-X (layer %i-%i) ... ', layer, layer+1)
            orb1   = findobj(orbitals,'Layer',layer,  'Spin',spin,'l',1);
            orb2   = findobj(orbitals,'Layer',layer+1,'Spin',spin,'l',1);
            norbs1 = size(orb1,1);
            norbs2 = size(orb2,1);
            coor1  = [reshape([orb1(:).Pcentre],[3,norbs1])'];
            coor2  = [reshape([orb2(:).Pcentre],[3,norbs2])'];
            hindex1 = [[orb1(:).Ham_index]'];
            hindex2 = [[orb2(:).Ham_index]'];

            ir = 0;
            supercell2  = zeros(ncell1^2*norbs2, 3);
            super_hidx2 = zeros(ncell1^2*norbs2, 3);
            for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                    ir = ir + 1;
                    supercell2((ir-1)*norbs2+1:ir*norbs2, :) = ...
                        coor2 + n1*mcell(1,:) + n2*mcell(2,:);
                    super_hidx2((ir-1)*norbs2+1:ir*norbs2, :) = ...
                        [hindex2(:), n1*ones(norbs2,1), n2*ones(norbs2,1)];
                end
            end
            dist2 = sqrt( bsxfun(@plus, sum(coor1.^2,2), sum(supercell2.^2,2)') ...
                          - 2*(coor1*supercell2') );

            inter_nn2_list = cell(norbs1,1);
            inter_nn2_cen  = cell(norbs1,1);

            parfor ii = 1:norbs1
                iorb = hindex1(ii);
                dl   = find(dist2(ii,:) <= aXX2);
                n1v  = super_hidx2(dl, 2);  n2v = super_hidx2(dl, 3);
                inter_nn2_list{ii} = super_hidx2(dl, 1)';
                inter_nn2_cen{ii}  = orb_Pcentre_all(super_hidx2(dl,1), :) ...
                                     - orb_Pcentre_all(iorb, :) ...
                                     + n1v.*mcell(1,:) + n2v.*mcell(2,:);
            end

            for ii = 1:norbs1
                iorb = hindex1(ii);
                orbitals(iorb).Inter_nn2_list   = inter_nn2_list{ii};
                orbitals(iorb).Inter_nn2_centre = inter_nn2_cen{ii};
            end
            fprintf('done [%.1fs]\n', toc(t_inter))

            clear orb1 orb2

            % ── X-M (d-p) interlayer coupling ─────────────────────────────
            t_inter = tic; fprintf('    interlayer X-M (layer %i-%i) ... ', layer, layer+1)
            orb1  = findobj(orbitals,'Layer',layer,  'Spin',spin, ...
                            {'Rel_index',9,'-or','Rel_index',6, ...
                             '-or','Rel_index',7,'-or','Rel_index',8});
            orb2  = findobj(orbitals,'Layer',layer+1,'Spin',spin, ...
                            {'Rel_index',3,'-or','Rel_index',6, ...
                             '-or','Rel_index',7,'-or','Rel_index',8});
            diorb = [orb1(:).Ham_index];
            djorb = [orb2(:).Ham_index];
            clear orb1 orb2

            s1 = size(diorb,2);
            s2 = size(djorb,2);

            inter_mx_list = cell(s1,1);
            inter_mx_cen  = cell(s1,1);

            parfor ii = 1:s1
                iorb     = diorb(ii);
                iorb_rel = orb_rel_idx_all(iorb);
                in4 = 0;
                in4_list   = zeros(1,1);
                in4_centre = zeros(1,3);
                for jj = 1:s2
                    jorb     = djorb(jj);
                    jorb_rel = orb_rel_idx_all(jorb);
                    if((iorb_rel == 9 && jorb_rel ~= 3) || ...
                       (iorb_rel ~= 9 && jorb_rel == 3))
                        for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                            for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                                if(iorb_rel == 9)
                                    shift = n1*mcell(1,:) + n2*mcell(2,:) ...
                                            + orb_centre_all(jorb,:);
                                    d = shift - orb_Pcentre_all(iorb,:);
                                else
                                    shift = n1*mcell(1,:) + n2*mcell(2,:) ...
                                            + orb_Pcentre_all(jorb,:);
                                    d = shift - orb_centre_all(iorb,:);
                                end
                                if(norm(d) < aXM)
                                    in4 = in4 + 1;
                                    in4_list(in4)      = jorb;
                                    in4_centre(in4,:)  = d;
                                end
                            end
                        end
                    end
                end
                inter_mx_list{ii} = in4_list(1:in4);
                inter_mx_cen{ii}  = in4_centre(1:in4,:);
            end

            for ii = 1:s1
                iorb = diorb(ii);
                orbitals(iorb).Inter_mx_n_list   = inter_mx_list{ii};
                orbitals(iorb).Inter_mx_n_centre = inter_mx_cen{ii};
            end
            fprintf('done [%.1fs]\n', toc(t_inter))

        end  % for layer (interlayer)
    end

end  % for spin

clear in_list in_centre inn_list inn_centre innn_list innn_centre ...
      in2_list in2_centre in3_list in3_centre dis1 dis2
clear liorb ljorb pdiorb pdjorb dpiorb dpjorb diorb djorb
end
