function y = sawtooth(x,c,r)
% y = sawtooth(x,c,r)
% sawtooth periodic function with contraction duration c and relaxation duration r
% restrictions: c and r should be > 0 and < pi
% example:
% x = -10:0.01:10;
% c = 1;
% r = 2;
% y = sawtooth(x,c,r);
% plot(x,y)

y = zeros(size(x));
for i = 1:length(x)
    y(i) = sawtoothPeriod(mod(x(i)+pi,2*pi)-pi,c,r);
end

y = (((y+1)/2).^2).*2-1; % smooth valley
    
end

function y = sawtoothPeriod(x,c,r)
    if x < -c
        y = 1;
    elseif x < 0
        y = -2/c*x-1;
    elseif x < r
        y = 2/r*x-1;
    else
        y = 1;
    end
end