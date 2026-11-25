function [durations, signs] = pulseToDurations(direction, angle,RabiFreq,RabiX0)
    %RabiFreq = 12;
    direction = direction / norm(direction);
    
    if isequal(direction, [1, 0, 0])
        angles = [angle, 0, 0];  % Corrected syntax
    end
    if isequal(direction, [0, 1, 0])
        angles = [0 ,angle, 0];
    else
        R = vrrotvec2mat([direction, angle]);
        angles = rotm2eul_XYX(R);
    end
    
    % Assuming angleToDuration function is defined elsewhere
    anglesTime = angleToDuration(angles,RabiFreq, RabiX0);
    durations = abs(anglesTime);
    signs = sign(anglesTime);
end


%         function [f_combined_pulse, total_duration] = createPulse( direction, angle, phi_0)     % not used currently
%             if ~exist('phi_0', 'var')
%                 phi_0 = 0;
%             end
%             
%             % Normalize the direction vector
%             direction = direction / norm(direction);
%             % Calculate the corresponding Bloch sphere angles (theta, phi)
%             R = vrrotvec2mat([direction, angle]);
%             [alpha, beta, gamma] = rotm2eul(R, 'XYX');
%             
%             alphaTime = RoundDurationByFrequencies(angleToDuration(alpha), AWG.SampleRate);
%             betaTime = RoundDurationByFrequencies(angleToDuration(beta), AWG.SampleRate);
%             gammaTime = RoundDurationByFrequencies(angleToDuration(gamma), AWG.SampleRate);
%             
%             total_duration = alphaTime + betaTime + gammaTime;
%             dt = 1 / sampleRate;
%             t = 0 : dt : total_duration - dt; 
%             
%             % Define pulse functions
%             w = 2*pi * frequency;
%             f_pulse = cell(3, 1);
%             f_pulse{1} = @(t) (0 <= t & t < alphaTime) .* cos(w*t - phi_0);
%             f_pulse{2} = @(t) (alphaTime <= t & t < (alphaTime + betaTime)) .* cos(w*t + pi/2 - phi_0);
%             f_pulse{3} = @(t) ((alphaTime + betaTime) <= t & t <  (alphaTime+betaTime+gammaTime)) .* cos(w*t - phi_0);
%             
%             % Create combined pulse function
%             f_combined_pulse = sum(cellfun(@(f) f(t), f_pulse));
%         end


function time = angleToDuration(angle,RabiFreq,RabiX0)
        %RabiFreq=12; RabiX0=0;
    time = angle / (2*pi*RabiFreq) + RabiX0;
end

function angle = durationToAngle( time)
    RabiFreq=12;
    RabiX0 =0;
    angle = 2*pi*RabiFreq * (time + RabiX0);
end


function [eul] = rotm2eul_XYX( R)
    
    psi = mod(atan2(-R(2,1),R(3,1)),pi);
    theta = mod(acos(R(1,1)),2*pi);
    phi = mod(atan2(R(1,2), R(1,3)),pi);
    
    % Pack the Euler angles into a 1-by-3 row vector
    eul = [psi, theta, phi];
end


