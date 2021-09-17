function paramFSI = saveParamFSI(wingDesign,simParam,wingModelStructure,wingModelAero,storeAllData)

paramFSI.wingParams.b       = wingDesign.span;
paramFSI.wingParams.chord   = wingDesign.chord;
paramFSI.wingParams.S_mw    = wingDesign.span*wingDesign.chord;
paramFSI.wingParams.airfoil = wingDesign.airfoil;
paramFSI.cog                = wingModelStructure.wP.cog;

paramFSI.nModes                                 = simParam.nmodes + 2;
paramFSI.analysisParams.num_airfoil_nodes_panel = simParam.num_airfoil_nodes_panel;
paramFSI.analysisParams.n_seg_PM                = simParam.n_seg_PM;
paramFSI.analysisParams.rho                     = simParam.rho;
paramFSI.analysisParams.simtol                 = simParam.simtol;

paramFSI.K_MODES    = wingModelStructure.K_MODES;
paramFSI.M_MODES    = wingModelStructure.M_MODES;
paramFSI.C_MODES    = wingModelStructure.C_MODES;

paramFSI.y_segLimit       = wingModelAero.y_segLimit      ;
paramFSI.y_segMidpoint    = wingModelAero.y_segMidpoint   ;
paramFSI.xMesh_loop       = wingModelAero.xMesh_loop      ;
paramFSI.yMesh_loop       = wingModelAero.yMesh_loop      ;
paramFSI.zMesh_loop       = wingModelAero.zMesh_loop      ;
paramFSI.cLdistr_yPos     = wingModelAero.cLdistr_yPos    ;

paramFSI.M1               = wingModelAero.M1              ;
paramFSI.N1               = wingModelAero.N1              ;
paramFSI.TPSU_X           = wingModelAero.TPSU_X          ;
paramFSI.TPSU_Y           = wingModelAero.TPSU_Y          ;
paramFSI.TPSU_Z           = wingModelAero.TPSU_Z          ;
paramFSI.TPSD_X           = wingModelAero.TPSD_X          ;
paramFSI.TPSD_Y           = wingModelAero.TPSD_Y          ;
paramFSI.TPSD_Z           = wingModelAero.TPSD_Z          ;
paramFSI.IDWuX            = wingModelAero.IDWuX           ;
paramFSI.IDWuY            = wingModelAero.IDWuY           ;
paramFSI.IDWuZ            = wingModelAero.IDWuZ           ;
paramFSI.IDWdX            = wingModelAero.IDWdX           ;
paramFSI.IDWdY            = wingModelAero.IDWdY           ;
paramFSI.IDWdZ            = wingModelAero.IDWdZ           ;
paramFSI.viscPre          = wingModelAero.viscPre;
paramFSI.K_alpha_ind_ELLT = wingModelAero.K_alpha_ind_ELLT;

paramFSI.fiACTMax         = wingModelAero.fiACTMax;
paramFSI.fiACTMin         = wingModelAero.fiACTMin;
paramFSI.fiACT_LAMax      = wingModelAero.fiACT_LAMax;
paramFSI.fiACT_LAMin      = wingModelAero.fiACT_LAMin;

%% save parameters for FSI and ROM
save(['data/parsim_FSI_', wingDesign.airfoil, '.mat'],'paramFSI');


%% save all parameters
if storeAllData
    save(['data/par_all_FSI_', wingDesign.airfoil, '.mat'],'paramFSI','wingDesign','simParam','wingModelStructure','wingModelAero');
end
