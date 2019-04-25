% load the image we will experiment with
I = imresize(double(rgb2gray(imread('lena.png'))),[256 256]);

% build the Laplacian pyramid of this image with 6 levels
depth = 6;
L = laplacianpyr(I,depth);

% compute the quantization of the Laplacian pyramid
bins = [16,32,64,128,128,256]; % number of bins for each pyramid level
LC = encoding(L,bins);

% compute the entropy for the given quantization of the pyramid
ent = pyramident(LC);

% Use the collapse command of the Lab 3 to recover the image
Ic = collapse(LC);

% compute the snr for the recovered image
snr_c = compute_snr(I,Ic);

% use the code from Lab 2 to compute an approximation image with 
% the same level of compression approximately
[rows,cols] = size(I);
n_0 = rows*cols;
M = n_0/8;
Id = decompress(compress(I,sqrt(M)));
snr_d = compute_snr(I,Id);

% plot the resulting images
subplot(1,3,1); 
imshow(I,[]); title('Original image');
subplot(1,3,2); imshow(Ic,[]); 
title('Laplacian Encoding'); xlabel(['SNR = ' num2str(snr_c)]);
subplot(1,3,3); imshow(Id,[]); 
title('Fourier Approximation'); xlabel(['SNR = ' num2str(snr_d)]);

function [Id] = decompress(Fcomp)

    % Input:
    % F: the compressed version of the image
    % Output:
    % Id: the approximated image

    % Please follow the instructions in the comments to fill in the missing commands.    
    
    % 1) Apply the inverse FFT shift (MATLAB command ifftshift)
    I1=ifftshift(Fcomp);
    % 2) Compute the inverse FFT (MATLAB command ifft2)
    I2=ifft2(I1);
    % 3) Keep the real part of the previous output
    Id=real(I2);
end
function [Fcomp] = compress(I,M_root)

    % Input:
    % I: the input image
    % M_root: square root of the number of coefficients we will keep
    % Output:
    % Fcomp: the compressed version of the image

    % Please follow the instructions in the comments to fill in the missing commands.    
    
    % 1) Perform the FFT transform on the image (MATLAB command fft2).
    I1=fft2(I);
    % 2) Shift zero-frequency component to center of spectrum (MATLAB command fftshift).
    I2=fftshift(I1);
    % We create a mask that is the same size as the image. The mask is zero everywhere, 
    % except for a square with sides of length M_root centered at the center of the image.
    [rows,cols] = size(I);
    idx_rows = abs((1:rows) - ceil(rows/2)) < M_root/2 ; 
    idx_cols = abs((1:cols)- ceil(cols/2)) < M_root/2 ; 
    M = (double(idx_rows')) * (double(idx_cols));
    
    % 3) Multiply in a pointwise manner the image with the mask.
    Fcomp=M.*I2;
end
function ent = pyramident(LC)  
    ent = 0;                % initialization of entropy
    [r, c] = size(LC{1});
    pixI = r*c;             % number of pixels in the original image
    
    for i = 1:numel(LC)
         [r, c] = size(LC{i});
         pixi = r*c;
        e=entropy(LC{i});
        ent=ent+e*pixi/pixI;
    end
    
end
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
function I = collapse(L)

    % Input:
    % L: the Laplacian pyramid of an image
    % Output:
    % I: Recovered image from the Laplacian pyramid

    % Please follow the instructions to fill in the missing commands.
    
    depth = numel(L);
    
    % 1) Recover the image that is encoded in the Laplacian pyramid
    for i = depth:-1:1
        if i == depth
            % Initialization of I with the smallest scale of the pyramid
            I = L{i};
        else
            % The updated image I is the sum of the current level of the
            % pyramid, plus the expanded version of the current image I.
            I = L{i}+expand(I);
        end
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
function snr = compute_snr(I, Id)

    % Input: 
    % I: the original image
    % Id: the approximated (noisy) image
    % Output:
    % snr: signal-to-noise ratio
    
    % Please follow the instructions in the comments to fill in the missing commands.    

    % 1) Compute the noise image (original image minus the approximation)
    I1=I-Id;
    % 2) Compute the Frobenius norm of the noise image
    fapp=norm(I1,'fro')
    % 3) Compute the Frobenius norm of the original image
    f = norm(I,'fro')
    % 4) Compute SNR
    snr=-20*log10((fapp)/f);
    SNR=-20*log(fapp/f)
end