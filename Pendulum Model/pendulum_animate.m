speed = 10;

x=[];z=[];

figure;
for i=1:speed:length(u)
   
   ui = u(:,i);
   x(:,i) = ui(1:7:end);
   z(:,i) = ui(3:7:end);

end

t2 = 0:1e-2:t(end);
x2 = interp1(t,x',t2)';
z2 = interp1(t,z',t2)';

for i=1:length(t2)
   clf;
   plot(x2(:,i),z2(:,i),':ok')
   xlim([-1.25 1.25]);
   ylim([-1.5 1]);
   grid on
   drawnow
   
   %disp(t2(i))
end