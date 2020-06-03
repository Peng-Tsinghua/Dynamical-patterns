function [value,isterminal,direction] = nbusEventsFcn(t,y,~,~,~,~,K,~,~)
n=length(K);
% The value that we want to be zero: maxV=1e6, max(abs(the))=500
value = [max(y(n+1:2*n))-1e6;max(abs(y(1:n)))-500]; 
isterminal = [1;1];  % Halt integration 
direction = [0;0];   % The zero can be approached from either direction
end