
%--------------------------------------------------------------------------
% This function segments the 3_D hyperspectral data into superpixels by 
% angle-based SLIC.
% A: M*N*L hyperspectral images
% k: desired number of superpixels
% m: scaling factor

% group: label matrix of superpixel segmentation
% num: the number of superpixels
%--------------------------------------------------------------------------


function [group,num] = hyperspectral_superpixels(A,k,m)

[M,N,L] = size(A);
s = ceil((M*N/k)^0.5);
r = ceil(M/s);
w = ceil(N/s);
k = r*w;
group = zeros(M,N);
center = zeros(k,L+2);
dist = inf*ones(M,N);
filt = 5;

%% initialize center
for i=1:r
    for j=1:w
        if (i<r)
            x = (i-1)*s+fix(s/2);
        else
            x = (i-1)*s+fix(rem(M,s)/2);
        end
        if (j<w)
            y=(j-1)*s+fix(s/2);
        else
            y=(j-1)*s+fix(rem(N,s)/2);
        end
        
        for p = x-filt:x+filt
            if(p >= 1)&&(p <= M)
                for q = y-filt:y+filt
                    flag = inf;
                    derta = 0;
                    if(q >= 1)&&(q <= N)
                        if p-1>0&&q-1>0
                            derta = derta+norm(reshape(A(p-1,q-1,:),1,L) - reshape(A(p,q,:),1,L));
                        end
                         if p-1>0&&q+1<N
                            derta = derta+norm(reshape(A(p-1,q+1,:),1,L) - reshape(A(p,q,:),1,L));
                         end
                         if p+1<M&&q-1>0
                            derta = derta+norm(reshape(A(p+1,q-1,:),1,L) - reshape(A(p,q,:),1,L));
                         end
                         if p+1<M&&q+1<N
                            derta = derta+norm(reshape(A(p+1,q+1,:),1,L) - reshape(A(p,q,:),1,L));
                         end
                                         
                         if flag > derta
                             flag = derta;
                             x_temp = p;y_temp = q;
                         end
                         
                    end
                end
            end
        end
        
        x = x_temp;y = y_temp;
        z = A(x,y,:);
        z = reshape(z,1,L);
        center((i-1)*w+j,:)=[z x y];
    end
end

%% update centers and labels
iter = 1000;
error = k;
for t = 1:iter
    if error < 4
        break;
    end
    error = 0;
    
    for i = 1:k
        for u = center(i,L+1)-s:center(i,L+1)+s
            if(u >= 1)&&(u <= M)
                for v = center(i,L+2)-s:center(i,L+2)+s
                    if(v>=1)&&(v<=N)
                        dc = dis(reshape(A(u,v,:),1,L),center(i,[1:L]));
                        ds = norm([u,v] - center(i,[L+1 L+2]));
                        d = dc + ds*m/s;
                        if d<dist(u,v)
                            dist(u,v) = d;
                            group(u,v) = i; 
                            error = error + 1;
                        end
                    end
                end
            end
        end
    end 
    
    temp_center = zeros(k,L+2);
    num_group = zeros(k,1);
    for p = 1:M
        for q = 1:N
        i = group(p,q);
        num_group(i) = num_group(i) + 1;
        temp_center(i,:) = [reshape(A(p,q,:),1,L) p q] + temp_center(i,:);
        end
    end
    for i = 1:k
    center(i,:) = ceil(temp_center(i,:)/num_group(i));
    end
    
end

%% merge small superpixels
for i=1:k
    bw = zeros(M,N);
    for p =1:M
        for q = 1:N
            if group(p,q) == i
                bw(p,q)=1;
            end
        end
    end
    [L, num] = bwlabel(bw, 4);
    for t = 1:num
        [rr, cc] = find(L==t);
        c1=size(rr,1);
        if c1>0&&c1<10
            for g=1:c1(1)
                if rr(1)-1>=1
                    group(rr(g),cc(g))=group(rr(1)-1,cc(1));
                elseif cc(1)-1>=1
                    group(rr(g),cc(g))=group(rr(1),cc(1)-1);
                elseif cc(1)+1<=N
                    group(rr(g),cc(g))=group(rr(1),cc(1)+1);
                elseif rr(1)+1<=M
                    group(rr(g),cc(g))=group(rr(1)+1,cc(1));
                end
            end
        end
    end
end
k = max(max(group));
num = k;
i = 1;
% while(i<k)
%     temp = size(find(group == i),1);
%     if norm(temp) == 0 &&i == k
%         num = num-1;
%         break
%         
%     elseif  norm(temp) == 0
%        temp_1 = find(group > i);
%        group(temp_1) = group(temp_1) - 1;
%        num = num - 1;
%        i = i-1;
%        k = k-1;
%     end
%     i = i+1;
%     
% end

%% visualize the resultes of superpixel segmentation
temp = A(:,:,[99,67,30]);% the channel of false clolor[30,6,99] [99,67,30][20,100,140]);
for i = 1:3
temp(:,:,i) = temp(:,:,i)/max(max(temp(:,:,i)))*255;
end
temp = uint8(temp);
figure;
imshow(temp,'InitialMagnification','fit');

%surf(group)
figure;
boundary = boundarymask(group);
imshow(imoverlay(temp,boundary,'red'),'InitialMagnification','fit');
end