%Hfunc.m
%Calculating H function
function H = Hfunc(J)
%     H = sqrt(det(J * J')); %Redundant manipulator
    H = sqrt(det(J' * J)); %Underactuated manipulator
end

