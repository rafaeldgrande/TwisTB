%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Serial version of opt_build_and_diag  %%
%% (no Parallel Computing Toolbox needed) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Hmat,tb_bands,tb_vecs,gradH_x] = opt_build_and_diag_H_3rdnn(task,Hmat,...
              orbitals,norbs,tot_norbs,multilayer,nspin,nlayer,neigs,...
              all_kpts,knum_tot,interlayer_int,lsoc,T,onsite,compute_eigvecs,...
              pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
              pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,flipped,...
              read_ham,save_workspace,reduced_workspace,num_workers,gradient,lambda,...
              pc,dirname,parallel,write_ham,onsite_moire,ef_strength,g1)

if(~read_ham)
   diagH = zeros(tot_norbs);
   temp_norbs = 0;
   for layer = 1 : nlayer
      for iorb = 1:norbs(layer)
          orb1 = orbitals(iorb + temp_norbs);
          if( orb1.Rel_index == 7 || orb1.Rel_index == 8)
             dr =-orb1.dr*g1;
          else
             dr = 0.0;
          end
          diagH(orb1.Ham_index,orb1.Ham_index) = diagH(orb1.Ham_index,orb1.Ham_index) + onsite(orb1.Rel_index,layer) + dr;
          if(ef_strength ~= 0.0)
             if(orb1.l == 2)
                z_coor = dot(orb1.Centre,[0 0 1]);
             elseif(orb1.l == 1)
                z_coor = dot(orb1.Pcentre,[0 0 1]);
             end
             diagH(orb1.Ham_index,orb1.Ham_index) = diagH(orb1.Ham_index,orb1.Ham_index) + ef_strength*z_coor;
          end
      end
      temp_norbs = sum(norbs(1:layer));
   end
   clear orb1;

   Mxz = generate_mask(orbitals,tot_norbs,nlayer,flipped);
   if(reduced_workspace)
      Mxz = sparse(Mxz);
   end

   fprintf('--> Starting serial diagonalisation of TB Hamiltonian ... \n')

   % Serial loop over k-points
   tb_bands_global = zeros(neigs, knum_tot);
   if(compute_eigvecs)
      tb_vecs_global = zeros(tot_norbs, neigs, knum_tot);
   end
   if(write_ham)
      Hmat = cell(1, knum_tot);
   end
   if(gradient)
      gradH_x_out = cell(1, knum_tot);
   end

   for ik = 1 : knum_tot
      if(reduced_workspace)
         H = sparse(diagH);
         gH = sparse(zeros(tot_norbs,tot_norbs));
      else
         H = diagH;
         gH = zeros(tot_norbs);
      end
      k = all_kpts(ik,:);
      if(gradient)
         [H,gH] = build_H_3rdnn(H,orbitals,Mxz,k, ...
                                 nlayer,multilayer,nspin,interlayer_int,lsoc,T, ...
                                 pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
                                 pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,gradient,gH);
      else
         H = build_H_3rdnn(H,orbitals,Mxz,k, ...
                                 nlayer,multilayer,nspin,interlayer_int,lsoc,T, ...
                                 pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
                                 pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,gradient);
      end

      if(write_ham)
         Hmat{ik} = H;
      end
      if(gradient)
         gradH_x_out{ik} = gH;
      end

      if(reduced_workspace)
         H = H + lambda*sparse(eye(size(H)));
         if (compute_eigvecs)
            [tb_vecs_global(:,:,ik),D] = eigs(H,neigs,'SM');
            [tb_bands_global(:,ik),ind] = sort(diag(real(D)));
            tb_vecs_global(:,:,ik) = tb_vecs_global(:,ind,ik);
         else
            D = eigs(H,neigs,'SM');
            [tb_bands_global(:,ik),ind] = sort(real(D));
         end
      else
         if (compute_eigvecs)
            [V,D] = eig(H,'vector');
            [tb_bands_global(:,ik),ind] = sort(real(D));
            tb_vecs_global(:,:,ik) = V(:,ind);
         else
            D = eig(H,'vector');
            [tb_bands_global(:,ik),ind] = sort(real(D));
         end
      end
      fprintf('k-vector number %i completed \n',ik);
   end

   fprintf('done\n\n')

   if(~write_ham)
      Hmat = 0.0;
   end
   if(gradient)
      gradH_x = gradH_x_out;
      clear gradH_x_out;
   else
      gradH_x = 0.0;
   end

   if(save_workspace)
      fprintf('--> Saving Hamiltonian matrix in H.mat ... ')
      save(fullfile(dirname,'H'),'Hmat');
      fprintf('done\n\n')
   end
   if(save_workspace && task==5)
      fprintf('--> Saving Hamiltonian gradient matrix in gradH.mat ... ')
      save(fullfile(dirname,'gradH'),'gradH_x');
      fprintf('done\n\n')
   end
   clear diagH H gH

else
   fprintf('--> Starting serial diagonalisation of TB Hamiltonian on read Hamiltonian ... ')

   tb_bands_global = zeros(neigs, knum_tot);
   if(compute_eigvecs)
      tb_vecs_global = zeros(tot_norbs, neigs, knum_tot);
   end

   gradH_x = 0.0;

   for ik = 1 : knum_tot
      H = cell2mat(Hmat(ik));
      if(reduced_workspace)
         H = H + lambda*sparse(eye(size(H)));
         if (compute_eigvecs)
            [tb_vecs_global(:,:,ik),D] = eigs(H,neigs,'SM');
            [tb_bands_global(:,ik),ind] = sort(diag(real(D)));
            tb_vecs_global(:,:,ik) = tb_vecs_global(:,ind,ik);
         else
            D = eigs(H,neigs,'SM');
            [tb_bands_global(:,ik),ind] = sort(real(D));
         end
      else
         if (compute_eigvecs)
            [V,D] = eig(H,'nobalance','vector');
            [tb_bands_global(:,ik),ind] = sort(real(D));
            tb_vecs_global(:,:,ik) = V(:,ind);
         else
            D = eig(H,'nobalance','vector');
            [tb_bands_global(:,ik),ind] = sort(real(D));
         end
      end
      fprintf('k-vector number %i completed \n',ik)
   end
   fprintf('done\n\n')
end

clear D;
if(~reduced_workspace && compute_eigvecs)
    clear V;
end
if(ismatrix(T))
   clear T;
end
if(lsoc && iscell(Hsoc))
   clear Hsoc
end

% Return results (no gather needed for serial)
tb_bands = real(tb_bands_global);
clear tb_bands_global;
if(compute_eigvecs)
   tb_vecs = tb_vecs_global;
   clear tb_vecs_global;
else
   tb_vecs = 0.0;
end

end
