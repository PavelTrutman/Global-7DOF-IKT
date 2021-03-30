% Pavel Trutman
% pavel.trutman@cvut.cz

addpath('../');

% store your desired pose in the variable Mh
MhFeasible =   [[41/105, -88/105, 8/21, -11807/100]; [8/21, 11/21, 16/21, 389/100]; [-88/105, -16/105, 11/21, -24931/100]; [0, 0, 0, 1]];
MhInfeasible = [[43/93, 20/93, -80/93, -3841/100]; [-20/93, -85/93, -32/93, 791/20]; [-80/93, 32/93, -35/93, 13219/50]; [0, 0, 0, 1]];
MhFailed =     [[69/95, -6/19, 58/95, -1089/5]; [42/95, 17/19, -6/95, 15977/100]; [-10/19, 6/19, 15/19, -4607/50]; [0, 0, 0, 1]];
Mh = MhFeasible;

% store your theta_hat in the variable qHat
qHat = zeros(7, 1);

% set your weights of the objective function
w = ones(7, 1);
w = w./sum(w);

% set power of the objective function n and the relaxation order
n = 4;
relaxOrder = 2;

% load kinematic parameters of the manipulator
% choose one of the random manipulators
manipulatorId = 1;
load(['manipulator.', num2str(manipulatorId), '.mat']);

% unknowns
mpol('c', 7);
mpol('s', 7);

% polynomials
p = DHT_GP(M, c, s, 3)*DHT_GP(M, c, s, 4)*DHT_GP(M, c, s, 5) - iDHT_GP(M, c, s, 2)*iDHT_GP(M, c, s, 1)*Mh*iDHT_GP(M, c, s, 7)*iDHT_GP(M, c, s, 6);
p = p(1:3, :);
p = transpose(p); % to match equation numbering with Maple
I = p(:);

% add cos(th_i)^2 + sin(th_i)^2 = 1 equations
I = [I; c.^2 + s.^2 - 1];

% add joint limit inequalities
Iineq = [-(c + 1).*tan(M.thetaLimLow./2) + s; (c + 1).*tan(M.thetaLimHigh./2) - s];

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
    errorC = norm(MhAngles(1:3, 4) - Mh(1:3, 4));
    errorR = acos(1/2*(trace(Mh(1:3, 1:3)\MhAngles(1:3, 1:3))-1));

    fprintf(['Position error: ', num2str(errorC), ' mm\n']);
    fprintf(['Rotation error: ', num2str(rad2deg(errorR)), ' deg\n']);
  end
else
  angles = zeros(7, 0);
end