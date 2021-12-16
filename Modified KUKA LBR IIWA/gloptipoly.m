% Pavel Trutman
% pavel.trutman@cvut.cz

addpath('../');

% store your desired pose in the variable Mh
MhFeasible =   [    0.3427    0.3200    0.8832  239.6992;
                   -0.4003    0.9003   -0.1709  355.4527;
                   -0.8499   -0.2950    0.4367  321.9285;
                         0         0         0    1.0000];
MhInfeasible = [   -0.5853    0.8026    0.1151  214.4608;
                    0.0837   -0.0814    0.9932 -119.0223;
                    0.8065    0.5909   -0.0195  342.5528;
                         0         0         0    1.0000];
MhFailed =     [    0.4556    0.3670   -0.8110 -284.5424;
                   -0.2883   -0.8011   -0.5245 -199.1473;
                   -0.8422    0.4728   -0.2592  809.9252;
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

% load kinematic parameters of the manipulator
load('manipulator.mat');

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
    errorC = norm(MhAngles(1:3, 4) - Mh(1:3, 4));
    errorR = acos(1/2*(trace(Mh(1:3, 1:3)\MhAngles(1:3, 1:3))-1));

    fprintf(['Position error: ', num2str(errorC), ' mm\n']);
    fprintf(['Rotation error: ', num2str(rad2deg(errorR)), ' deg\n']);
  end
else
  angles = zeros(7, 0);
end
