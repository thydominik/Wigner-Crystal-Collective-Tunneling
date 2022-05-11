function [R] = rotation_matrix(delta_phi, e, f, N)
    R = zeros(3,3,N);
%     for i = 1:length(e)-1
%         C       = cross(e(:,i),e(:,i+1));
%         D       = dot(e(:,i),e(:,i+1));
%         NP0   = norm(e(:,i));
%         
%         if ~all(C==0) % check for colinearity    
%             Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
%             R(:,:,i) = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
%         else
%             R(:,:,i) = sign(D) * (norm(e(:,i)) / NP0) ; % orientation and scaling
%         end
%         
%     end
%     R(:,:,end) = 0;
   
    for i = 1:N
        R1 = eye(3) - (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') + (e(:,i)*e(:,i)' + f(:,i)*f(:,i)') * cos((delta_phi(i)));
        R2 = (e(:,i)*f(:,i)' - (f(:,i)*e(:,i)')) * sin((delta_phi(i)));
        R(:,:,i) =  (R1 + R2)';
    end

end

