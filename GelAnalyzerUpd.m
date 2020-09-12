%% Import

if exist('im_import',"var")
    im_import = imread([filedir, filename]);
else
    [filename, filedir] = uigetfile('*.*');
    im_import = imread([filedir, filename]);
end



min_bkgd_im = SubtractBackground(im_import);

band = FindBands(min_bkgd_im);

band_table = struct2table(band);


imshow(imcomplement(min_bkgd_im));
hold on
for b = 1:length(band)
plot(band(b).trace(:,2), band(b).trace(:,1), 'm', 'LineWidth', 2);
end

%imshowpair(im_import, imcomplement(min_bkgd_im), 'montage')



%% Functions



function [img, baseline] = SubtractBackground(UINT8)
    % Function takes input image object and returns uint8 with
    % subtracted background based on UI-selected background
    % Need to make it compatable with uint16

    [rows, columns] = size(UINT8);

    % Conditional to ensure all calculations are done with raw (dark)
    % image, not inverted image
    
    int8_ratio = sum(UINT8,"all")/(256*(rows*columns));
    if int8_ratio > 1 || int8_ratio < 0
        error('Image is not uint8')
    end
    
    % Need to also do uint16
    if int8_ratio < .5
        light = imcomplement(UINT8);
        raw = UINT8;
    else
        light = UINT8;
        raw = imcomplement(UINT8);
    end
    
    imobj = imshow(light);
    
    bkgd_rect = drawrectangle(imobj.Parent);
    
    % This could probably be neater
    xmax = imobj.XData(2);
    ymax = imobj.YData(2);
    
    
    x1 = bkgd_rect.Position(1);
    y1 = bkgd_rect.Position(2);
    mx = bkgd_rect.Position(3);
    my = bkgd_rect.Position(4);
    
    x2 = (x1 + mx);
    y2 = (y1 + my);
    
    % Should double check the math here sometime
    x_row1 = round(rows*x1/xmax);
    x_row2 = round((rows*x2/xmax)-2); % Need to -2 so that index never exceeds row bounds
    y_col1 = round(columns*y1/ymax);
    y_col2 = round((columns*y2/ymax)-2);
    
    bkgd_arr = raw(x_row1:x_row2, y_col1:y_col2);
    avg = mean(bkgd_arr, "all");
    
    avg_arr = ones(rows, columns) * double(avg);
    
    img = raw - uint8(avg_arr);
    baseline = avg;
end

function band = FindBands(GelImage)

    
    band = struct('BandNumber', 0, 'pixels',[], 'trace', [], 'volume', [], 'midpoint', [], 'Rf', 0);

    log_im = bwareaopen(GelImage, 20); % Should conn be hard-coded?
    sz = size(GelImage);
    [traces, ~] = bwboundaries(log_im, 'noholes');


    band_cell = bwconncomp(log_im).PixelIdxList';

    midpts = FindMidpoints(band_cell, sz);
    rfs = FindRfs(midpts, sz);
    vols = FindVolumes(band_cell, GelImage);
    bandnums = AssignBandNum(midpts);
    
    for k = 1:length(band_cell)
        band(k).BandNumber = bandnums(k);
        band(k).pixels = band_cell{k};
        band(k).trace = traces{k};
        band(k).volume = vols(k);
        band(k).midpoint = midpts{k};
        band(k).Rf = rfs(k);
    end


        function volumes = FindVolumes(band_indicies, gel_im)

            volumes{20,1} = [];

            for i = 1:length(band_indicies)
                vol = 0;

                for j = 1:length(band_indicies{i})
                    vol = vol + double(gel_im(band_indicies{i}(j)));
                end

                volumes{i} = vol;
            end
            volumes = volumes(~cellfun('isempty',volumes));
        end

        function midpoints = FindMidpoints(band_indicies, gel_size)
        % Input is a cell of linear indicies from band_cell
        % Output is a cell of corresponding midpoints in [row, column] form

            mp_cell{20,1} = [];

            for i = 1:length(band_indicies)
                i_of_i = roundn(length(band_indicies{i})/2,0);
                mp = band_indicies{i}(i_of_i);
                [row, column] = ind2sub(gel_size, mp);
                mp_cell{i,1} = [row, column];
            end

            midpoints = mp_cell(~cellfun('isempty',mp_cell));

        end

        function Rfs = FindRfs(midpoints, gel_size)
            Rfs = zeros(20,1);
            total_rows = gel_size(1,1);

            for i = 1:length(midpoints)
                mp_row = midpoints{i}(1,1);
                rf = mp_row/total_rows;
                Rfs(i) = rf;
            end
            Rfs = Rfs(arrayfun(@any,Rfs));
        end

        function bandnums = AssignBandNum(midpoints)
            %Function is kind of a mess, could probably be more concise
            
            keys = zeros(length(midpoints), 1);
            bandnums = zeros(length(midpoints));
            
            for i = 1:length(midpoints)
                keys(i) = midpoints{i}(1,1);
            end
            
            keys = sort(keys);
            values = 1:length(keys);
            bandnums_map = containers.Map(keys, values);
            
            for i =1:length(keys)
            bandnums(i) = bandnums_map(midpoints{i}(1,1));
            end
            
        end
    
end

