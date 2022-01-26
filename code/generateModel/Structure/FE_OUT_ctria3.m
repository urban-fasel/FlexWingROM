function meshOut = FE_OUT_ctria3(mirrorWing, EID, PID, G1, G2, G3)

%%
%
%  By Giulio Molinari
%
%  Modified by Urban Fasel
%
%
% EID		Element identification number. (0 < Integer < 100,000,000)
% PID		Property identification number of a PSHELL, PCOMP, PCOMPG or PLPLANE
%			entry. (Integer > 0; Default = EID)
% Gi		Grid point identification numbers of connection points. (Integers > 0, all unique)


% Check for consistency in vector sizes
nel = numel(EID);
if (numel(PID) > 1)
	if (numel(PID) ~= nel)
		error('FE_OUT_ctria3:VectorSizes', 'Inconsistent number of elements in PID with respect to element IDs.');
	end
else
	PID = repmat(PID,nel,1);
end


% meshOut array
meshOut = [];
for ii = 1:length(EID)
    meshOut = [meshOut; [EID(ii), PID(ii), G1(ii), G2(ii), G3(ii)]];
end


% mirror wing
dID = 1000000;
EID = EID + dID;
G1 = G1 + dID;
G2 = G2 + dID;
G3 = G3 + dID;

for i = 1:nel
    if ~isempty(mirrorWing)
        if find(G1(i)==mirrorWing(:,1)+dID)
            G1(i) = G1(i)-dID;
        end
        if find(G2(i)==mirrorWing(:,1)+dID)
            G2(i) = G2(i)-dID;
        end
        if find(G3(i)==mirrorWing(:,1)+dID)
            G3(i) = G3(i)-dID;
        end
    end
end

for ii = 1:length(EID)
    meshOut = [meshOut; [EID(ii), PID(ii), G1(ii), G2(ii), G3(ii)]];
end