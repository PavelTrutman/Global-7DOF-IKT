% Pavel Trutman
% pavel.trutman@cvut.cz

addpath('../');

% store your desired pose in the variable Mh
MhFeasible.M =       [[41/105, -88/105, 8/21, -11807/100]; [8/21, 11/21, 16/21, 389/100]; [-88/105, -16/105, 11/21, -24931/100]; [0, 0, 0, 1]];
MhFeasible.maple =   'MhR := Matrix(3, 4, [[41/105,-88/105,8/21,-11807/100],[8/21,11/21,16/21,389/100],[-88/105,-16/105,11/21,-24931/100]]);';
MhInfeasible.M =     [[43/93, 20/93, -80/93, -3841/100]; [-20/93, -85/93, -32/93, 791/20]; [-80/93, 32/93, -35/93, 13219/50]; [0, 0, 0, 1]];
MhInfeasible.maple = 'MhR := Matrix(3, 4, [[43/93,20/93,-80/93,-3841/100],[-20/93,-85/93,-32/93,791/20],[-80/93,32/93,-35/93,13219/50]]);';
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
% choose one of the random manipulators
manipulatorId = 1;
load(['manipulator.', num2str(manipulatorId), '.mat']);

% prepare for Maple
eqFile = tempname;
system(['touch ', eqFile]);

% run Maple
mapleCMD = [maplePath, ' -c ''', Mh.maple, ''' -c ''eqFile := \`', eqFile, '\`'' -c ''manId := ', num2str(manipulatorId), ''' -c ''read \`symReduction.mpl\`'' -c ''quit'''];
[mapleStatus, output] = system(mapleCMD, '-echo');
assert(mapleStatus == 0, output);

% format the equations
system(['sed -i -e ''1 s/\([cs]\)\([1-7]\)/\1(\2)/g'' ', eqFile]);

% unknowns
mpol('c', 7);
mpol('s', 7);

% parse the equations
eqStr = splitlines(fileread(eqFile));
eval(eqStr{1});
I = cg';

% clean up
delete(eqFile);

% add cos(th_i)^2 + sin(th_i)^2 = 1 equations
I = [I; c.^2 + s.^2 - 1];

% add joint limit inequalities
Iineq = [-(c + 1).*tan(M.thetaLimLow./2) + s; (c + 1).*tan(M.thetaLimHigh./2) - s; c - max(cos(M.thetaLimLow), cos(M.thetaLimHigh))];

% add variable constraints
Iineq = [Iineq; 1 - c.^2; 1 - s.^2];

% objective function
f = sum(2*w.*(-c.*cos(qHat) -s.*sin(qHat) + 1 + 1).^n);

mset('yalmip', true);
mset(sdpsettings('solver', 'mosek'));

mset('testol', 1e-2);
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
  sSol = squeeze(double(s));
  cSol = squeeze(double(c));
  angles = atan2(sSol, cSol);

  % show errors
  for j = 1:size(cSol, 3)
    MhAngles = double(DHFKT(M, angles(:, j)));

    % compute errors
    errorC = norm(MhAngles(1:3, 4) - Mh.M(1:3, 4));
    errorR = acos(1/2*(trace(Mh.M(1:3, 1:3)\MhAngles(1:3, 1:3))-1));

    fprintf(['Position error: ', num2str(errorC), ' mm\n']);
    fprintf(['Rotation error: ', num2str(rad2deg(errorR)), ' deg\n']);
  end
else
  angles = zeros(7, 0);
end
