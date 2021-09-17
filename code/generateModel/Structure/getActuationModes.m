function [x_ASET_L, x_ASET_R, ACTparam] = getActuationModes(wingParams, wingDesign, numRibs, FE_Grid, K_ASET, DOF_RSET, forceACTinp)

% Calculate actuation mode

nDOF = 6*size(FE_Grid.ID, 1);

forces_ACTl = zeros(nDOF/6, 6);
forces_ACTr = zeros(nDOF/6, 6);

ActFront_ID = [];
ActRear_ID = [];
ActFront_IDR = [];
ActRear_IDR = [];

% left wing
ActFront_Pos1 = zeros(numRibs-2,3);
ActRear_Pos1 = zeros(numRibs-2,3);
nV = zeros(numRibs-2,3);

for iAct = numRibs-2:-2:numRibs-2*wingDesign.nRibsC
    ActFront_ID(iAct) = find(FE_Grid.X(:) == wingParams.actuation.actuationSparAP(iAct,1) & FE_Grid.Y(:) == wingParams.actuation.actuationSparAP(iAct,2) & FE_Grid.Z(:) == wingParams.actuation.actuationSparAP(iAct,3));
    ActRear_ID(iAct) = find(FE_Grid.X(:) == wingParams.actuation.actuationsdvAP(iAct,1) & FE_Grid.Y(:) == wingParams.actuation.actuationsdvAP(iAct,2) & FE_Grid.Z(:) == wingParams.actuation.actuationsdvAP(iAct,3));
    ActFront_Pos1(iAct,:) = [FE_Grid.X(ActFront_ID(iAct)), FE_Grid.Y(ActFront_ID(iAct)), FE_Grid.Z(ActFront_ID(iAct))];
    ActRear_Pos1(iAct,:) = [FE_Grid.X(ActRear_ID(iAct)), FE_Grid.Y(ActRear_ID(iAct)), FE_Grid.Z(ActRear_ID(iAct))];
    nV(iAct,:) = (ActFront_Pos1(iAct,:)-ActRear_Pos1(iAct,:))./norm(ActFront_Pos1(iAct,:)-ActRear_Pos1(iAct,:));
    forces_ACTl(ActFront_ID(iAct),1:3) = -nV(iAct,:)*forceACTinp.L;
    forces_ACTl(ActRear_ID(iAct),1:3) = nV(iAct,:)*forceACTinp.L;
end

forces_ACTl_transp = forces_ACTl';
forces_ACT_L = forces_ACTl_transp(:);
forces_ACT_L(DOF_RSET) = [];
x_ASET_L = K_ASET\forces_ACT_L; %Find actuation mode Left


% right wing
ActFront_Pos1R = zeros(numRibs-2,3);
ActRear_Pos1R = zeros(numRibs-2,3);
nVR = zeros(numRibs-2,3);
for iAct = numRibs-2:-2:numRibs-2*wingDesign.nRibsC
    ActFront_IDR(iAct) = find(FE_Grid.X(:) == wingParams.actuation.actuationSparAP(iAct,1) & FE_Grid.Y(:) == wingParams.actuation.actuationSparAP(iAct,2) & FE_Grid.Z(:) == -wingParams.actuation.actuationSparAP(iAct,3));
    ActRear_IDR(iAct) = find(FE_Grid.X(:) == wingParams.actuation.actuationsdvAP(iAct,1) & FE_Grid.Y(:) == wingParams.actuation.actuationsdvAP(iAct,2) & FE_Grid.Z(:) == -wingParams.actuation.actuationsdvAP(iAct,3));
    ActFront_Pos1R(iAct,:) = [FE_Grid.X(ActFront_IDR(iAct)), FE_Grid.Y(ActFront_IDR(iAct)), FE_Grid.Z(ActFront_IDR(iAct))];
    ActRear_Pos1R(iAct,:) = [FE_Grid.X(ActRear_IDR(iAct)), FE_Grid.Y(ActRear_IDR(iAct)), FE_Grid.Z(ActRear_IDR(iAct))];
    nVR(iAct,:) = (ActFront_Pos1R(iAct,:)-ActRear_Pos1R(iAct,:))./norm(ActFront_Pos1R(iAct,:)-ActRear_Pos1R(iAct,:));
    forces_ACTr(ActFront_IDR(iAct),1:3) = -nVR(iAct,:)*forceACTinp.R;
    forces_ACTr(ActRear_IDR(iAct),1:3) = nVR(iAct,:)*forceACTinp.R;
end

forces_ACTr_transp = forces_ACTr';
forces_ACT_R = forces_ACTr_transp(:);
forces_ACT_R(DOF_RSET) = [];
x_ASET_R = K_ASET\forces_ACT_R; %Find actuation mode Right


% store actuation IDs in ACTparam struct
ACTparam.ActFront_ID = ActFront_ID;
ACTparam.ActRear_ID = ActRear_ID;
ACTparam.ActFront_IDR = ActFront_IDR;
ACTparam.ActRear_IDR = ActRear_IDR;
ACTparam.nV = nV;
ACTparam.nVR = nVR;

end