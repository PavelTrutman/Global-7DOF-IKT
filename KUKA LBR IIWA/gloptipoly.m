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
Iineq = [-(c + 1).*tan(M.thetaLimLow./2) + s; (c + 1).*tan(M.thetaLimHigh./2) - s; c - min(cos(M.thetaLimLow), cos(M.thetaLimHigh))];

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
