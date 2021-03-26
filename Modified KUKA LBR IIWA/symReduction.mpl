with(ListTools):
with(LinearAlgebra):
with(PolynomialTools):
with(combinat, choose):
with(Groebner):
with(MatrixPolynomialAlgebra):
with(CodeGeneration):
interface(rtablesize=10):
interface(warnlevel=0):
Digits:=30:
RatResolution:=10000:

kernelopts(numcpus=1):

read(`../DH.mpl`);

# mechanism description
MechParams := {
a1=0,     a2=100,  a3=0,     a4=0,    a5=0,     a6=0,    a7=0,
d1=340,   d2=0,    d3=400,   d4=0,    d5=400,   d6=0,    d7=126,
A1=-Pi/2, A2=Pi/2, A3=-Pi/2, A4=Pi/2, A5=-Pi/2, A6=Pi/2, A7=0};

Mechanism := {op(MechParams),
              lambda1=CS2CSR(subs(MechParams,A1),RatResolution)[1],mu1=CS2CSR(subs(MechParams,A1),RatResolution)[2],
              lambda2=CS2CSR(subs(MechParams,A2),RatResolution)[1],mu2=CS2CSR(subs(MechParams,A2),RatResolution)[2],
              lambda3=CS2CSR(subs(MechParams,A3),RatResolution)[1],mu3=CS2CSR(subs(MechParams,A3),RatResolution)[2],
              lambda4=CS2CSR(subs(MechParams,A4),RatResolution)[1],mu4=CS2CSR(subs(MechParams,A4),RatResolution)[2],
              lambda5=CS2CSR(subs(MechParams,A5),RatResolution)[1],mu5=CS2CSR(subs(MechParams,A5),RatResolution)[2],
              lambda6=CS2CSR(subs(MechParams,A6),RatResolution)[1],mu6=CS2CSR(subs(MechParams,A6),RatResolution)[2],
              lambda7=CS2CSR(subs(MechParams,A7),RatResolution)[1],mu7=CS2CSR(subs(MechParams,A7),RatResolution)[2]};

# import end effector pose
Mh := Matrix(4,4,[[lx,mx,nx,rx],[ly,my,ny,ry],[lz,mz,nz,rz],[0,0,0,1]]):
MhV := ImportMatrix(MhFile)[2][1..3, 1..4];
MhR := E2ER(MhV,RatResolution);
MhRV := map(x->x[1]=x[2],convert(<convert(Mh[1..3,1..4],Vector)|convert(MhR,Vector)>,listlist));

# equations
lrE:=M2L((dhT(3).dhT(4).dhT(5)- dhInv(dhT(2)).dhInv(dhT(1)).<<MhR>,<0|0|0|1>>.dhInv(dhT(7)).dhInv(dhT(6)))[1..3,1..4]):
eq:=simplify(subs([op(Mechanism),op(MhRV)],M2L(lrE)));

eqX:=numer~(normal~(([op(eq),s1^2+c1^2-1,s2^2+c2^2-1,s3^2+c3^2-1,s4^2+c4^2-1,s5^2+c5^2-1,s6^2+c6^2-1,s7^2+c7^2-1]))):

infolevel[GroebnerBasis] := 0:
B := Basis(eqX, tdeg(op(indets(eqX))));

if not evalb(B = [1]) then

  idxDegTwo := ListTools[SearchAll](2, map(degree, B)):
  eqQ := B[[idxDegTwo]]:

  BNew := Basis(eqQ, tdeg(op(indets(eqQ))));

  if evalb(B = BNew) then

    for i from 1 to nops(eqQ) do
      eqQ[i]:= evalf(eqQ[i]/LeadingCoefficient(eqQ[i], tdeg(op(indets(eqQ))))):
    end do:
    eqMatlab := Matlab(eqQ, output=string):
    fd := fopen(eqFile, WRITE, TEXT):
    FileTools[Text][WriteString](fd, eqMatlab):
    fclose(fd):
  end if;
end if;