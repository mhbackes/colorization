function colorImageRgb = colorization(grayImageRgb, markedImageRgb)
    
    [n, m, ~] = size(grayImageRgb);
    imageSize = n * m;

    grayImage = rgb2ntsc(grayImageRgb);
    markedImage = rgb2ntsc(markedImageRgb);
    
    diff = abs(grayImageRgb - markedImageRgb) ~= 0;
    markedPixels = diff(:,:,1) | diff(:,:,2) | diff(:,:,3);
    
    colorImage(:, :, 1) = grayImage(:, :, 1);
    
    D = speye(imageSize);
    %W = sparse(imageSize, imageSize);
    wRow = zeros(imageSize,1);
    wCol = zeros(imageSize,1);
    wVal = zeros(imageSize,1);
    wSize = 0;
    
    index = reshape(1:imageSize, n, m);
    
    for y = 1:m
        for x = 1:n
            if(~markedPixels(x,y))
                Y = grayImage(x,y,1);
                [nRows, nCols, nY] = neighbours(x, y, grayImage);           
                weight = w(Y,nY);
                weightSum = sum(weight);
                r = index(x,y);
                for i = 1:size(nY,1)
                    s = index(nCols(i),nRows(i));
                    wSize = wSize + 1;
                    wRow(wSize) = r;
                    wCol(wSize) = s;
                    wVal(wSize) =  weight(i) / weightSum;
                end
            end
        end
    end

    wRow = wRow(1:wSize);
    wCol = wCol(1:wSize);
    wVal = wVal(1:wSize);
    W = sparse(wRow,wCol,wVal,imageSize,imageSize);
    
    A = D - W;
    U = reshape(markedImage(:,:,2) .* markedPixels, imageSize, 1);
    V = reshape(markedImage(:,:,3) .* markedPixels, imageSize, 1);
    newU = A \ U;
    newV = A \ V;
    
    colorImage(:,:,2) = reshape(newU, n, m);
    colorImage(:,:,3) = reshape(newV, n, m);
    
    colorImageRgb = ntsc2rgb(colorImage);
end

function [rows, cols, Y] = neighbours(x, y, image)
    k = 1;
    rows = zeros(8,1);
    cols = zeros(8,1);
    Y = zeros(8,1);
    [n, m, ~] = size(image);
    for i = max(1,x-1):min(n,x+1)
        for j = max(1,y-1):min(m,y+1)
            if(i~=x || j~=y)
                rows(k) = j;
                cols(k) = i;
                Y(k) = image(i, j, 1);
                k = k + 1;
            end
        end
    end
    rows = rows(1:(k-1),:);
    cols = cols(1:(k-1),:);
    Y = Y(1:(k-1),:);
end

function weight = w(Yr, Ys)
    Y = vertcat(Ys, Yr);
    %variance = var(Y, 1);
    variance=mean((Y-mean(Y)).^2)*.6;
    mgv=min((Ys-Yr).^2);
    if (variance<(-mgv/log(0.01))) %bruxarias para evitar divisão por zero
        variance=-mgv/log(0.01);
    end
    if(variance < 0.000002)
        variance = 0.000002;
    end
    weight = exp(-(Yr - Ys).^2/variance);
end
