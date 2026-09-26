function MPC_objfun = MPC_objfun(x, Phi_extend, Gamma_extend, r_p, Q, R, delta_U, Phiref, REF)
    Y =  Phi_extend*x + Gamma_extend*delta_U + Phiref*REF;    
 %persistent done
 %  if isempty(done)
 %      disp('Phi_extend=');
 %      disp(Phi_extend);
 %      disp('Phiref');
 %      disp(Phiref);
 %      disp('REF');
 %      disp(REF);
 %      done = true;
 %  end
    MPC_objfun = 1/2*(r_p-Y)'*Q*(r_p-Y) + 1/2*delta_U'*R*delta_U;
end

