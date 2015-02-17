function mustart = AvoidValleyJump(x,gamma,mustart)
% % valley jumping avoidance routine
    if x>3
        valley = diff(gamma(:,x-3:x)'); 
    
        for v=1:size(valley,2)
            if sign(valley(1,v)) ~= sign(valley(2,v)) && ...
                    sign(valley(2,v)) ~= sign(valley(3,v)) % if oscillating estimations
                mustart(v) = mustart(v)*.5; % lower the step confidence
            else
%                 mustart(v) = mustart_init(v);
            end
        end
    end
end