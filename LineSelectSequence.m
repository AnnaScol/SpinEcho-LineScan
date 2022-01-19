function kSpace_Slice = LineSelectSequence(line_number,xSteps, ySteps,gradAmp, FOV,durPE,timestepsREADOUTPrephaseEvent)

gamma = 2*pi*42.577*10^6;

kSpace_Slice  = zeros(ySteps,total_num_lines_x);
gradPEmax = pi * (ySteps-1)/(gamma*FOV) * 1/durPE;
gradAmp(2,timestepsREADOUTPrephaseEvent) =  (1-(line_number-1)/((ySteps-1)/2)) * gradPEmax; %Y gradients in Tesla per meter

tic
for trIndex=1:total_num_lines_x
    disp(trIndex);  

    k=1:xSteps;  
    j=1:ySteps;

        for i=1:zSteps            
            r = reshape(pos(:,k,j,i),[3,xSteps*ySteps]);
            mT = zeros(xSteps*ySteps,1);
            mZ = ones(xSteps*ySteps,1);
            
            dB0 = (gradAmp(:,1)'*r)'; 
            [mT,mZ] =  bloch(dt, dB0,rfPulse(1), reshape(T1(k,j),[xSteps*ySteps,1]),...
                                                 reshape(T2(k,j),[xSteps*ySteps,1]), mT, mZ);  
                
            for t=2:nTimeSteps 
                dB0 = (gradAmp(:,t)'*r)';
                [mT,mZ] =  bloch(dt, dB0,rfPulse(t), reshape(T1(k,j),[xSteps*ySteps,1]),...
                                                     reshape(T2(k,j),[xSteps*ySteps,1]), mT, mZ);
                
                if(adc(1,t)>0)
                    PDt = PD(k,j);
                    kSpace_Slice(round(adc(t)),trIndex) = kSpace_Slice(round(adc(t)),trIndex)+sum(mT.*reshape(PDt(k,j),[xSteps*ySteps,1]));
                end
                
            end  %end of time loop                     
        end %end z steps
end %end of TR loop

end