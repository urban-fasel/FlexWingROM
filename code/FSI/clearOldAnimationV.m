function clearOldAnimationV(paramFSI)
home = pwd;
sure = input(['sure to delete the existing ' paramFSI.animationName ' FSI animation data ?']);
if sure==0
    error('Stop');
end
try
    cd(strcat(home,filesep,'data'))
    try
        rmdir(strcat(paramFSI.wingParams.airfoil,filesep,'Animation',filesep,paramFSI.animationName),'s')
    catch
    end
    cd(home)
catch
end
end