function [common,shift,orient]=reorder_common_TOT(p,Q)
        
        n=size(Q,1);
        shift = zeros(n,1);
        common = zeros(n,1);
        orient = zeros(n,1);

        %Count overlap
        
        for i = 1:n
            q=Q(i,:);
            c=sum(sum(p'==q));
            common(i)=c;  
        end

        
        %make sure circshift(p,orient) -- circshift(q,shift) is
        %correct t-t sauter structure
        for i = 1:n
            if (common(i)==3)
               shift(i)= ((Q(i,:)==p))*[0;1;2];
               orient(i)=0;
            elseif common(i)==2
                M=(p'==q);
                shift(i)=(sum(M,1)==0)*[2;1;0];
                orient(i)=[2,1,0]*(sum(M,2)==0);
                %q=[q(2) q(1) q(3)];
            else
                M=p'==q;
                shift(i)= sum(M,1)*[0;2;1];
                orient(i)=[0,2,1]*sum(M,2);
            end
        end
        
 end