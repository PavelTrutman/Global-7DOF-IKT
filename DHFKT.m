% Pavel Trutman
% pavel.trutman@cvut.cz

function M = DHFKT(M, th)
  M = DHT(M, th, 1)*DHT(M, th, 2)*DHT(M, th, 3)*DHT(M, th, 4)*DHT(M, th, 5)*DHT(M, th, 6)*DHT(M, th, 7);
end