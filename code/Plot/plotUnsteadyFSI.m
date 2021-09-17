function plotUnsteadyFSI(simOut,simUnsteadyParam)

time = simUnsteadyParam.dT:simUnsteadyParam.dT:simUnsteadyParam.time;

figure
plot(time,simOut.cL)
xlabel('time, s')
ylabel('cL')

figure
plot(time,simOut.cD)
xlabel('time, s')
ylabel('cD')

figure
plot(time,simOut.cRoll)
xlabel('time, s')
ylabel('cRoll')

figure
plot(time,simOut.bendingModeAmplitude)
xlabel('time, s')
ylabel('bending mode amplitude')

end