function audit=inelast_accept_step(audit,abs_residual,rel_residual,state,V,gammadot,theta,sigma_e,retries,step_time,postFolder)
% Audit the original pressure tolerance before any accepted-state projection.
pressure_ok=isreal([abs_residual,rel_residual]) && all(isfinite([abs_residual,rel_residual])) && ...
    abs_residual<=audit.pressure_tolerance && rel_residual<=audit.pressure_tolerance;
state_ok=isreal(state) && all(isfinite(state)) && V>0 && all(gammadot>0) && all(theta>0) && all(sigma_e>0);
if ~pressure_ok || ~state_ok
    audit.nonconverged_steps=audit.nonconverged_steps+double(~pressure_ok);
    audit.invalid_state_steps=audit.invalid_state_steps+double(~state_ok);
    audit.failure_time=step_time;
    save(fullfile(postFolder,'inelast_acceptance_failure.mat'),'audit','abs_residual','rel_residual','retries','step_time','state','-v7');
    error('inelast:InvalidAcceptedStep','Refusing invalid accepted step at %.16g s (pressure residuals %.6g, %.6g; tolerance %.6g; pressure_ok=%d, state_ok=%d).',step_time,abs_residual,rel_residual,audit.pressure_tolerance,pressure_ok,state_ok);
end
audit.accepted_step_count=audit.accepted_step_count+1;
audit.maximum_absolute_normalized_pressure_residual=max(audit.maximum_absolute_normalized_pressure_residual,abs_residual);
audit.maximum_relative_pressure_residual=max(audit.maximum_relative_pressure_residual,rel_residual);
audit.maximum_pressure_retries=max(audit.maximum_pressure_retries,retries);
audit.minimum_sigma_e=min(audit.minimum_sigma_e,min(sigma_e));
audit.minimum_gammadot=min(audit.minimum_gammadot,min(gammadot));
audit.minimum_theta=min(audit.minimum_theta,min(theta));
end
