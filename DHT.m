% Pavel Trutman
% pavel.trutman@cvut.cz

function M = DHT(M, th, i)
  [M1, M2] = DH(M, th, i);
  M = M1*M2;
end