% Decoupled Load Flow Method

s2 = 0.5 + 1i;
s3 = -1.5 + 0.3i;

a = [9.999 -5.006 -4.990;
    -5.006 13.014 -8.008
    -4.990 -8.008 13.026;];

b= [-24.999 9.997 15.002;
    9.997 -25.989 15.992;
    15.002 15.992 -30.988;];

ybus = complex(a,b);
v1 = 1;
del1 = 0;
v2 = 1;
del2 = 0;
v3 = 1;
del3 = 0;


x = [0;
    0;
    1;];

delta =[0;
    0;];

Pact = [-1.5;
    0.5;];

Qact = 1;

J1 = zeros(2,2);

for i = 1:3
    fprintf("del2 = %d\n", (del2));
    fprintf("del3 = %d\n", (del3));
    fprintf("v2 = %d\n\n", (v2));

    p2 = (abs(ybus(2,1)) * abs(v1) * abs(v2)) * cos(del2-del1- angle(ybus(2,1))) + (abs(ybus(2,2)) * abs(v2) * abs(v2)) * cos(angle(ybus(2,2))) + (abs(ybus(2,3)) * abs(v2) * abs(v3)) * cos(del2-del3-angle(ybus(2,3))); 
    q2= -(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))-abs(v2)*abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2))))-(abs(v2)*abs(v3)*abs(ybus(2,3)))*sin((angle(ybus(2,3))+del3-del2));
    p3=(abs(v3)*abs(v1)*abs(ybus(3,1))*cos((angle(ybus(3,1)))+del1-del3))+abs(v3)*abs(v3)*abs(ybus(3,3))*cos((angle(ybus(3,3))))+(abs(v2)*abs(v3)*abs(ybus(3,2)))*cos((angle(ybus(3,2))+del2-del3));

    J1(1,1)=(abs(v3)*abs(v1)*abs(ybus(3,1))*sin((angle(ybus(3,1)))+del1-del3))+(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
    J1(1,2)=-(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));

    J1(2,1)=(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((-angle(ybus(2,3)))-del3+del2));
    J1(2,2)=(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))+(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
    
    J4 =(abs(v1)*abs(ybus(2,1))*sin((-angle(ybus(2,1)))-del1+del2))+2*(abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2)))))+ (abs(ybus(3,2))*abs(v3)*sin(-angle(ybus(2,3))+del2-del3));
    Pcal = [p3;
        p2;];

    delP = Pact- Pcal;

    del_delta = inv(J1)*delP;
    delta = delta + del_delta;
    del3 = delta(1,1);
    del2 = delta(2,1);
   
    Qcal = q2;
    delQ = Qact- Qcal;
    delv = 1/J4 * delQ;

    v2 = v2 + delv;
    v1 = v2;
    fprintf("Iter  %d\n",i);
    disp("Pcal");
    disp(Pcal);
    disp("delP");
    disp(delP);
    disp("delta");
    disp(delta);
    fprintf("del2 = %d\n", rad2deg(del2));
    fprintf("del3 = %d\n", rad2deg(del3));
    fprintf("v2 = %d\n\n", (v2));

end
