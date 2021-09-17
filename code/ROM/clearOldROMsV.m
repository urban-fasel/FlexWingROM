function clearOldROMsV(paramFSI,V)
home = pwd;
sure = input(sprintf('sure to delete the existing Snapshot matrices and ROMs at V = %dm/s ?',V));
if sure==0
    error('Stop');
end
try
    cd(strcat(home,filesep,'data'))
    try
        rmdir(strcat(paramFSI.wingParams.airfoil,filesep,'ROM',filesep,...
            'V',num2str(V),filesep,'Snapshot'),'s')
    catch
    end
    cd(home)
catch
end
end