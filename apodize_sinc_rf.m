function result = apodize_sinc_rf(window,TBWF,FA,dt)
%flip angle in radians
    gamma = 2*pi*42.577*10^6;
    w = hann(window); 
    x = linspace(-TBWF,TBWF,window);
    rfPulse = sinc(x);

    summation = 0;
    % window and find area under
    for g = 1:window
        rfPulse(g) = rfPulse(g).*w(g);
        summation = summation+rfPulse(g)*(dt);
    end

    scalefactor = FA/(gamma*summation);
    rfPulse = scalefactor*rfPulse;
    excitation_checker_value = sum(rfPulse)*gamma*dt;
    result = rfPulse;
end
