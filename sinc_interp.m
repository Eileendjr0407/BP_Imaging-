function [out,x_out] = sinc_interp(in,x_in,x_new,N,win)
%
% sinc_interp
%
% Sinc-based (band-limited)interpolation.
%
% INPUTS
% in = data sequence to be interpolated
% x_in = vector of sample locations corresponding to data samples of 'in'.
% Must be uniformly spaced at some interval dx_in. Must be same
% length as 'in'.
% x_new = vector of desired sample locations. Must be uniformly spaced at
% some interval dx_out.
% N = order of interpolating sinc, in units of max(dx_in,dx_out).
% Must be odd.
% win = 1 if Hamming window applied to interpolation kernel, otherwise no
% window. (win=1 is recommended.)
%
% OUTPUTS
% out = interpolated data sequence corresponding to sample locations in
% x_out.
% x_out = sample locations of output vector. This will be a subset of the
% locations in x_new; relative span of x_new and x_in, and filter
% end effects, may limit x_out to not include some of the values in
% x_new.
%
% Mark A. Richards
% February 2007
if (mod(N,2) ~= 1)
disp(' ')
disp(' ** Error: sinc_interp : filter order not odd.')
disp([' ** Filter order input = ',int2str(N)])
disp(' ')
return
end
Nhalf = (N-1)/2;
d_in = x_in(2) - x_in(1);
% d_out = x_new(2) - x_new(1);
% Now figure out interpolating filter impulse response in continuous time.
% This will be a sinc function, possibly windowed. Bandwidth of the sinc
% LPF frequency response is based on the larger sampling interval of the
% two grids. Specifically, the unwindowed impulse response is h(x) =
% sin(pi*x/del)/(pi*x) and dt = max(dt1,dt2). To add to this, we specify
% how many sample increments we will go out on the tails, where 1 increment
% is dt; and then we also apply a hamming window of the same length. The
% Hamming formula in continuous time is w(t) = 0.54 + 0.46*cos(pi*t/dt).
% del = max(d_in,d_out); 
del=d_in;
% find the values within x_new that can be successfully interpolated from the
% values available in x_in, i.e. where end effects won't kill us.
% x_in
% x_new
% Nhalf
% x_new(1)-Nhalf*del
% x_in(1)
% x_new(end)+Nhalf*del
% x_in(end)
index = find( (x_new-Nhalf*del >= x_in(1) ) & ...
(x_new+Nhalf*del <= x_in(end)) );   %�ҵ�Ҫ��ֵ��λ��
if (isempty(index))
disp(' ')
disp(' ** Error: sinc_interp : Requested output samples cannot be interpolated')
disp(' ')
return
end
x_out = x_new(index);   % index ��ʾ����ֵ��λ��
out = zeros(size(x_out));
% step through the output samples one at a time, interpolating a value for
% each one from the input samples.
for k = 1:length(x_out)
% first find the span of the interpolating filter on the x axis
x_current = x_out(k);
x_low = x_current - Nhalf*del;
x_high = x_current + Nhalf*del;
% compute the *relative* position of each input sample within this span
% compared to the current output sample location; these will be the
% values at which the interpolating kernel filter response will be
% needed.
index_rel = find( (x_in >= x_low ) & ...
(x_in <= x_high) );
x_rel = x_in(index_rel) - x_current;
% Now compute and apply the sinc weights. First fix any spots where
% x_rel = 0; these will cause the sinc function to be undefined. Then
% add in the window, if used, and apply to the data to compute the
% output point.
trouble = find(x_rel==0);
if (~isempty(trouble))
    x_rel(trouble) = x_rel(trouble)+eps; % this will prevent division by zero
end
h = (sin(pi*x_rel/del)/pi./x_rel);
w = ones(size(h));
if (win == 1)
    w = 0.54 + 0.46*cos(pi*x_rel/del/(Nhalf+1));
    h = h.*w;
end
h = h/sum(h);
out(k) = sum( in(index_rel).*h);
end
end 