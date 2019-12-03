% Pavel Trutman
% pavel.trutman@cvut.cz


% store your desired pose in the variable Mh
%Mh = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

% path to the Maple binary
maplePath = 'maple';

load('manipulator.mat');

e4 = [1 2 3 5 6 7];

% prepare for Maple
MhFile = [tempname, '.mat'];
eqFile = tempname;
save(MhFile, 'Mh');
system(['touch ', eqFile]);

% run Maple
mapleCMD = [maplePath, ' -c ''MhFile := \`', MhFile, '\`'' -c ''eqFile := \`', eqFile, '\`'' -c ''read \`KUKA_LBR.mpl\`'' -c ''quit'''];
[status, output] = system(mapleCMD, '-echo');
%[status, output] = system(mapleCMD);

% format the equations
system(['sed -i -e ''1 s/\([cs]\)\([1-7]\)/\1(\2)/g'' ', eqFile]);
system(['sed -i -e ''3 s/\([cs]\)\([1-7]\)/\1Sol(\2, j)/g'' ', eqFile]);

% unknowns
mpol('c', 7);
mpol('s', 7);
ce4 = c(e4);
se4 = s(e4);

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
f = sum((ce4 - 1).^2) - 7;

mset('yalmip', true);
mset(sdpsettings('verbose', 0, 'solver', 'mosek'));
mset('verbose', 0);

mset('testol', 1e-1);
%mset('ranktol', 1e-1);

for relaxOrder = 2
  P = msdp(min(f), I == 0, Iineq >= 0, relaxOrder);
  [status, obj] = msol(P);
  if status == 1
    break
  end
end
fprintf(['status: ', num2str(status), '\n']);

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
  if ~dispatch
    for j = 1:numSol
      MhAngles = double(DHDKT(M, angles(:, j)));

      % compute errors
      errorC = norm(MhAngles(1:3, 4) - Mh(1:3, 4));
      errorR = acos(1/2*(trace(Mh(1:3, 1:3)\MhAngles(1:3, 1:3))-1));

      fprintf(['Position error: ', num2str(errorC), '\n']);
      fprintf(['Rotation error: ', num2str(rad2deg(errorR)), ' deg\n']);
    end
  end
else
  angles = zeros(7, 0);
end
