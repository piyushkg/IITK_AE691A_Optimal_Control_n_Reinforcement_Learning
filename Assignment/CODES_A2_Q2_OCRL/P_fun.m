function P_dot = P_fun(t, P, A, B, Q, R)
P = reshape(P, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
P_dot = - A.'*P - P*A + P*B* inv(R)*B.'*P - Q; %Determine derivative
P_dot = P_dot(:); %Convert from "n"-by-"n" to "n^2"-by-1




