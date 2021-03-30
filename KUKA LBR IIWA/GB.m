% Pavel Trutman
% pavel.trutman@cvut.cz

addpath('../');

% store your desired pose in the variable Mh
MhFeasible =   [    0.4405    0.6064    0.6620 -281.8821;
                    0.7712   -0.6331    0.0666  483.7939;
                    0.4595    0.4812   -0.7465  276.0602;
                         0         0         0    1.0000];
MhInfeasible = [    0.4528   -0.1592    0.8773  544.6913;
                    0.0416    0.9866    0.1576 -606.4869;
                   -0.8906   -0.0348    0.4534   39.5770;
                         0         0         0    1.0000];
MhFailed =     [    0.4293   -0.9028    0.0254  386.1023;
                    0.8547    0.4151    0.3118  -85.1766;
                   -0.2920   -0.1122    0.9498    4.5084;
                         0         0         0    1.0000];
Mh = MhFeasible;

% store your theta_hat in the variable qHat
qHat = zeros(7, 1);

% set your weights of the objective function
w = ones(7, 1);
w = w./sum(w);

% set power of the objective function n and the relaxation order
n = 4;
relaxOrder = 2;

% path to the Maple binary
maplePath = 'maple';

% load kinematic parameters of the manipulator
load('manipulator.mat');

e4 = [1 2 3 5 6 7];

% prepare for Maple
MhFile = [tempname, '.mat'];
eqFile = tempname;
save(MhFile, 'Mh');
system(['touch ', eqFile]);

% run Maple
mapleCMD = [maplePath, ' -c ''MhFile := \`', MhFile, '\`'' -c ''eqFile := \`', eqFile, '\`'' -c ''read \`symReduction.mpl\`'' -c ''quit'''];
[mapleStatus, output] = system(mapleCMD, '-echo');
assert(mapleStatus == 0, output);

% format the equations
system(['sed -i -e ''1 s/\([cs]\)\([1-7]\)/\1(\2)/g'' ', eqFile]);
system(['sed -i -e ''3 s/\([cs]\)\([1-7]\)/\1Sol(\2, j)/g'' ', eqFile]);

% unknowns
mpol('c', 7);
mpol('s', 7);
ce4 = c(e4);
se4 = s(e4);
qHate4 = qHat(e4);
we4 = w(e4);

% parse the equations
eqStr = splitlines(fileread(eqFile));
eval(eqStr{1});
I = cg';

% clean up
delete(MhFile);
delete(eqFile);

% add cos(th_i)^2 + sin(th_i)^2 = 1 equations
I = [I; ce4.^2 + se4.^2 - 1];

% add joint limit inequalities
Iineq = [-(ce4 + 1).*tan(M.thetaLimLow(e4)./2) + se4; (ce4 + 1).*tan(M.thetaLimHigh(e4)./2) - se4];

% add variable constraints
Iineq = [Iineq; 1 - ce4.^2; 1 - se4.^2];

% objective function
f = sum(2*we4.*(-ce4.*cos(qHate4) -se4.*sin(qHate4) + 1 + 1).^n);

mset('yalmip', true);
mset(sdpsettings('solver', 'mosek'));

mset('testol', 1e-1);
%mset('ranktol', 1e-1);

% solve in GloptiPoly
P = msdp(min(f), I == 0, Iineq >= 0, relaxOrder);
[status, obj] = msol(P);
if status == -1
  fprintf('\nThis is an infeasible pose.\n');
elseif status == 0
  fprintf('\nProblem for this pose can not be solved.\n');
elseif status == 1
  fprintf('\nPose successfully solved.\n');
end

if status == 1
  % recover solution
  numSol = size(double(ce4), 3);
  sSol = zeros(7, numSol);
  cSol = zeros(7, numSol);
  eval(eqStr{2});
  sSol(e4, :) = squeeze(double(se4));
  cSol(e4, :) = squeeze(double(ce4));
  cSol(4, :) = repmat(cg0, [1 numSol]);
  for j = 1:numSol
    eval(eqStr{3});
    sSol(4, j) = cg1;
  end
  angles = atan2(sSol, cSol);
  
  % show errors
  for j = 1:numSol
    MhAngles = double(DHFKT(M, angles(:, j)));

    % compute errors
    errorC = norm(MhAngles(1:3, 4) - Mh(1:3, 4));
    errorR = acos(1/2*(trace(Mh(1:3, 1:3)\MhAngles(1:3, 1:3))-1));

    fprintf(['Position error: ', num2str(errorC), ' mm\n']);
    fprintf(['Rotation error: ', num2str(rad2deg(errorR)), ' deg\n']);
  end
else
  angles = zeros(7, 0);
end