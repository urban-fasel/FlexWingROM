function clearOldROMsV(paramFSI,V,k)
home = pwd;
if k == 1
    sure = input(sprintf('sure to delete the existing Snapshot matrices and ROMs at V = %dm/s ?',V));
else
    sure = 1;
end
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