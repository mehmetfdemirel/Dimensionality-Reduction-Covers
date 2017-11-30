clc
eps = 0.1;

% theta = 0.001*(1/sqrt(113));
% theta = 0.001*(1/sqrt(13));
theta = 0.01;

x1 = 0.2; x2 = 0.3; x3 = 0.4;
qq = sqrt(2);

dim = 3;
npoints=100;

xvect = zeros(dim,1);

for i=1:dim
    xvect(i) = 0.1 * i * qq;
end

x = zeros(100,1); y = zeros(100,1); z = zeros(100,1);
circle = zeros(100,1);


M = zeros(dim,npoints);
% Create the manifold M
for k = 1:npoints
    
%     x(k) = sin(2*pi*(x1 - k*theta)) + sin(4*pi*(x1 - k*theta)) + sin(6*pi*(x1 - k*theta));
% 
%     y(k) = sin(2*pi*(x2 - k*theta)) + sin(4*pi*(x2 - k*theta)) + sin(6*pi*(x2 - k*theta));
% 
%     z(k) = sin(2*pi*(x3 - k*theta)) + sin(4*pi*(x3 - k*theta)) + sin(6*pi*(x3 - k*theta));

    for j=1:dim
        M(j,k) = sin(2*pi*(xvect(j) - k*theta)) + sin(4*pi*(xvect(j) - k*theta)) + sin(6*pi*(xvect(j) - k*theta));
    end

end
figure;
% scatter3(M(1,:),M(2,:),M(3,:),'filled');


clear U
col = 1;

% Find all the normalized secants in M i.e. U(M)

for i=1:npoints
    % vec1 = [x(i) y(i) z(i)]';
      vec1 = M(:,i);
    for j=(i+1):npoints
        % vec2 = [x(j) y(j) z(j)]';
          vec2 = M(:,j);
        y_vect = (vec1 - vec2) / norm(vec1-vec2);
        U(:,col) = y_vect;
        col = col + 1;
    end
end

figure;
scatter3(U(1,:), U(2,:), U(3,:));

% U = (unique(U', 'rows'))'; % Remove duplicates
clear probs
[orgh, orgw] = size(U);
iter_sim = 100;

cover_num = 0;
next = 1;
clear covers
%    row_size = 50;
% amat = randn(row_size,dim)*(1/sqrt(row_size));
% xxx = 0;
% nchoose2 = npoints * (npoints - 1) / 2;
% count = 0;
% 
% [height, width] = size(U);
% 
% 
% 
% for iter=1:1000
%       % amat = randn(300,500)*(1/sqrt(300));
%     rand1 = randi(nchoose2);
%     rand2 = randi(nchoose2);
%     isit = 1;
%     while 1
%         p = U(:,rand1);
%         q = U(:,rand2);
%         if p(1) > 0 && p(2) > 0 && p(3) > 0 && q(1) > 0 && q(2) > 0 && q(3) > 0
%             break;
%         end
%     end
%       %p_q = amat*(p-q);
%       p_minus_q = p-q;
%     if norm(p_minus_q) > 0.1
%         xxx = xxx + 1;
%     end
%     count = count + 1
% end

while 1
    eps_j = (1/100) * (1/power(2,cover_num));
    [height, width] = size(U);
    if cover_num ~= 0
        [height1, w_cov] = size(covers);
    end 
    probs = zeros(1,width);
    
    if cover_num == 0  %constructing C0
        for i=1:iter_sim
            phi = randn(2,dim) * (1/sqrt(2)) * 1/100;
            phiTU = phi * U;
            for j = 1:width
                if norm(phiTU(:,j)) > 1 + eps_j
                    probs(:,j) = probs(:,j) + 1;
                end
            end
        end
       
        [MM,I] = min(probs);
        covers(:,next) = U(:,I(1));
        U(:,I(1)) = [];
        next = next + 1;
    else % Cj where j>0
        probs = zeros(1,width);
        for i=1:width
            p = U(:,i);
            
            for j=1:w_cov
                
                q = covers(:,j);
                diff = p-q;

                for x=1:iter_sim/5
                    phi = randn(2,3) * (1/sqrt(2)) * 1/100;
                    phiTdiff = phi * diff;
                    
                    if norm(phiTdiff) > eps_j
                        probs(:,i) = probs(:,i) + 1;
                    end
                end
            end
        end
        target = power(2,cover_num) - power(2,cover_num-1);
        while target > 0
            cols = size(probs,2);
            P = randperm(cols);
            probs = probs(:,P); 
            [MM,I] = min(probs);
            probs(:,I(1)) = [];
            covers(:,next) = U(:,I(1));
            next = next + 1;
            U(:,I(1)) = [];
            [height, width] = size(U);
            target = target - 1;
        end
    end
    [height1, w_covxx] = size(covers);
    
  
    if width <= 0
        break;
    end
    cover_num
    % covers
    figure;
    scatter3(covers(1,:)', covers(2,:)', covers(3,:)');
    cover_num = cover_num + 1;
end




