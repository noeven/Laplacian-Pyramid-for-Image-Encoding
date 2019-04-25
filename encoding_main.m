% load the image we will experiment with
I = imresize(double(rgb2gray(imread('lena.png'))),[256 256]);

% build the Laplacian pyramid of this image with 6 levels
depth = 6;
L = laplacianpyr(I,depth);

% compute the quantization of the Laplacian pyramid
bins = [16,32,64,128,128,256]; % number of bins for each pyramid level
LC = encoding(L,bins);

function LC = encoding(L, bins)
    depth = numel(bins);
    LC = cell(1,depth);  
    for i = 1:depth
        if i == depth % blurred image in range [0, 256]
            edges =linspace(0,256,bins(i));
        else % difference image in range [-128,128]
            edges =linspace(-128,128,bins(i));
        end
        centers =zeros(1,bins(i)-1);
        for j=2:bins(i)
           centers(j-1)=0.5*(edges(j)+edges(j-1)); 
        end
        LC{i} =discretize(L{i},edges,centers);     
    end
end
function L = laplacianpyr(I,depth)
    L = cell(1,depth);
    G = gausspyr(I,depth);
    for i = 1:depth
        if i < depth
            L{i} = G{i}-expand(G{i+1});
        else
            L{i} = G{i};
        end
    end
    
end
function G = gausspyr(I,depth)

   G = cell(1,depth);
    for i = 1:depth
        if i == 1
            G{i} = I; % original image
        else
            G{i} = reduce(G{i-1});% reduced version of the previous level
        end
    end

end
function g = reduce(I)
    [m,n,x]=size(I);
    G=fspecial('gaussian',5,1);
    I2=imfilter(I,G,'conv');
    g(1:m/2,1:n/2,:)=I2(1:2:m, 1:2:n,:);

end
function g = expand(I)

   [m,n,x]=size(I);
    I2=zeros(2*m,2*n,x);
    I2(1:2:m*2, 1:2:n*2,:) = I(1:m,1:n,:);
    G=fspecial('gaussian',5,1);
    g=4*imfilter(I2,G,'conv');
end