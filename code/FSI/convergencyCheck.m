function converged = convergencyCheck(iteration, analysisParams)

displacement_converged = false;
pressure_converged = false;

%% Convergency on displacement
if iteration.iterations_error <= analysisParams.simtol
	displacement_converged = true;
end

%% Convergency on pressure
if iteration.iterations_error_cPrel <= analysisParams.simtol
    pressure_converged = true;
end

converged = displacement_converged && pressure_converged;
end