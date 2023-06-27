%NEWTON RAPHSON METHOD
clc;

a = [0 0.02 0.01;
    0.02 0 0.0125;
    0.01 0.0125 0];

b = [0 0.04 0.03;
    0.04 0 0.025;
    0.03 0.025 0;];

LineData = complex(a,b); 

ybus = calculateYbus(LineData , 3);

disp('ybus=');
disp(ybus);

BusTypes = [1 2 3;]; 
Base_mva = 100;


v1=1.05+0i;
v2=1+0i;
v3=1+0i;

del3=angle(v3);
del1=angle(v1);
del2=angle(v2);

%% phasor values of the voltage at buses 2 and 3

for i=1:3
    p2 = (abs(ybus(2,1)) * abs(v1) * abs(v2)) * cos(del2-del1- angle(ybus(2,1))) + (abs(ybus(2,2)) * abs(v2) * abs(v2)) * cos(angle(ybus(2,2))) + (abs(ybus(2,3)) * abs(v2) * abs(v3)) * cos(del2-del3-angle(ybus(2,3))); 
    q2=-(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))-abs(v2)*abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2))))-(abs(v2)*abs(v3)*abs(ybus(2,3)))*sin((angle(ybus(2,3))+del3-del2));
    p3=(abs(v3)*abs(v1)*abs(ybus(3,1))*cos((angle(ybus(3,1)))+del1-del3))+abs(v3)*abs(v3)*abs(ybus(3,3))*cos((angle(ybus(3,3))))+(abs(v2)*abs(v3)*abs(ybus(3,2)))*cos((angle(ybus(3,2))+del2-del3));
    q3 = (abs(ybus(3,1)) * abs(v1) * abs(v3)) * sin(del3-del1- angle(ybus(3,1))) + (abs(ybus(3,2)) * abs(v2) * abs(v3)) * sin(del3-del2-angle(ybus(3,2))) - (abs(ybus(3,3)) * abs(v3)*abs(v3)) * sin(angle(ybus(3,3))); 
    
    J(1,1)=(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))+(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
    J(1,2)=(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((-angle(ybus(2,3)))-del3+del2));
    J(1,3)=(abs(v1)*abs(ybus(2,1))*cos(del2-del1-angle(ybus(2,1))))+2*(abs(v2)*abs(ybus(2,2))*cos((angle(ybus(2,2)))))+(abs(v3)*abs(ybus(2,3))*cos(del2-del3-angle(ybus(2,3))));
    J(1,4)= (abs(v2)*abs(ybus(2,3))*cos((-angle(ybus(2,3)))-del3+del2));
    
    J(2,1)=-(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
    J(2,2)=(abs(v3)*abs(v1)*abs(ybus(3,1))*sin((angle(ybus(3,1)))+del1-del3))+(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
    J(2,3)=(abs(v3)*abs(ybus(3,2))*cos((angle(ybus(3,2)))+del2-del3));
    J(2,4) =(abs(v1)*abs(ybus(3,1))*cos((-angle(ybus(3,1))-del1+del3)))+(abs(v2)*abs(ybus(3,2))*cos((-angle(ybus(3,2)))-del2+del3))+ 2*(abs(ybus(3,3))* abs(v3) * cos(angle(ybus(3,3))));
    
    J(3,1)=(abs(v2)*abs(v1)*abs(ybus(2,1))*cos((-angle(ybus(2,1)))-del1+del2))+(abs(v2)*abs(v3)*abs(ybus(2,3))*cos((-angle(ybus(2,3)))+del2-del3));
    J(3,2)=-(abs(v2)*abs(v3)*abs(ybus(2,3))*cos((-angle(ybus(2,3)))+del2-del3));
    J(3,3)=(abs(v1)*abs(ybus(2,1))*sin((-angle(ybus(2,1)))-del1+del2))+2*(abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2)))))+ (abs(ybus(3,2))*abs(v3)*sin(-angle(ybus(2,3))));
    J(3,4) = abs(ybus(2,3))*abs(v2)*sin(del2-del3-angle(ybus(2,3)));

    J(4,1)= -(abs(ybus(3,2))*abs(v2)*abs(v3)*cos(del3-del2-angle(ybus(3,2))));
    J(4,2)=(abs(v3)*abs(v1)*abs(ybus(3,1))*cos((-angle(ybus(3,1)))-del1+del3))+(abs(v2)*abs(v3)*abs(ybus(3,2))*cos((-angle(ybus(3,2)))-del2+del3));
    J(4,3)= abs(ybus(3,2))*abs(v3)*sin(del3-del2-angle(ybus(3,2)));
    J(4,4)=(abs(v1)*abs(ybus(3,1))*sin(-(angle(ybus(3,1)))))+(abs(v2)*abs(ybus(3,2))*sin(-(angle(ybus(3,2)))))+2*(abs(v3)*abs(ybus(3,3))*sin((angle(ybus(3,3)))));

    
    A=[del2;del3;abs(v2);abs(v3)];
    M = [-2.556+0.5;
        -1.386+0.5;
        -1.102+1;
        -0.45+1.5;];
    X = [0;
        0;
        1;
        1;];
    delx = inv(J)*M;
    delx(1,1) = rad2deg(delx(1,1));
    delx(2,1) = rad2deg(delx(2,1));
    X = X + delx;
    del2 = X(1,1);
    del3 = X(2,1);
    v2 = X(3,1);
    v3 = X(4,1);
    fprintf("Iter %d\n\n",i);
    disp("Jacobian Matrix");
    disp(J);
    fprintf("\n");
    fprintf("V2 = %f\nv3 = %f\n\n\n",v2,v3);
end


%% Function to calculate Ybus
function Ybus = calculateYbus(LineData, numBuses)
    Ybus = zeros(numBuses, numBuses);
    for i = 1:numBuses
        for j = 1:numBuses
            if i ~= j
                Ybus(i,j) = -1 / LineData(i,j);
            else
                for k = 1:numBuses
                    if i ~= k
                        Ybus(i,j) = Ybus(i,j) + 1 / LineData(i,k);
                    end
                end
            end
        end
    end
end
