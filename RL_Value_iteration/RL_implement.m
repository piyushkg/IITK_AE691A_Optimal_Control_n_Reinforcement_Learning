R_left = [  -2,0,0,-2;
            -2,-2,-2,-2;
            -2,0,0,0];  
R_right = [  0,0,0,-2;
            -2,-2,0,-2;
            0,0,0,-2];   
R_up =    [  -2,-2,-2,-2;
            0,-2,0,-2;
            0,-2,0,0];   
R_down =    [ 0,-2,0,-2;
            0,-2,0,-2;
            -2,-2,-2,-2];
R = zeros(3,4,4);
R(:,:,1) = R_left; R(:,:,2) = R_right; R(:,:,3) = R_up; R(:,:,4) = R_down;
Value = zeros(3,4,5);

gamma = 0.8;
prob = 0.8;
reward = [ 0,0,0,1;
           0,0,0,-1;
           0,0,0,0];
iter = 500;
lt = 0; rt = 0; up = 0; dn = 0;

for i = 1:iter
    for j = 1:3 % rows in environment
        for k = 1:4 % cols in environment
            state = [j,k];
            av_act =  possible_actions(state, R);
            rest_prob = (1-prob)/(length(av_act)-1);
            for l = 1:length(av_act) % possible actions
                if av_act(l) == 1 % left move
                    lt = reward(j,k-1) + gamma * Value(j,k-1,5);
                elseif av_act(l) == 2 % right move
                    rt = reward(j,k+1) + gamma * Value(j,k+1,5);
                elseif av_act(l) == 3 % up move
                    up = reward(j-1,k) + gamma * Value(j-1,k,5);
                else % down move
                    dn = reward(j+1,k) + gamma * Value(j+1,k,5);
                end 
%             end
%             for l = 1:length(av_act) % possible actions
                if av_act(l) == 1 % left move
                    Value(j,k,1) = prob * lt + rest_prob * rt + rest_prob * up + rest_prob * dn;
                elseif av_act(l) == 2 % right move
                    Value(j,k,2) = rest_prob * lt + prob * rt + rest_prob * up + rest_prob * dn;
                elseif av_act(l) == 3 % up move
                    Value(j,k,3) = rest_prob * lt + rest_prob * rt + prob * up + rest_prob * dn;
                else % down move
                    Value(j,k,4) = rest_prob * lt + rest_prob * rt + rest_prob * up + prob * dn;
                end
%                 Value(j,k,5) = max(Value(j,k,1:4));
            end
            Value(j,k,5) = max(Value(j,k,1:4));
            lt = 0; rt = 0; up = 0; dn = 0;
        end
    end
end
dp = Value(:,:,5)+reward
imagesc(dp);
colorbar;

% imagesc(dp);            % Create a colored plot of the matrix values
% colormap(flipud(gray));  % Change the colormap to gray (so higher values are
%                          %   black and lower values are white)
% 
% textStrings = num2str(dp, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
% [x, y] = meshgrid(1:3,1:4);  % Create x and y coordinates for the strings
% hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
%                 'HorizontalAlignment', 'center');
% midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
% textColors = repmat(mat(:) > midValue, 1, 3);  % Choose white or black for the
%                                                %   text color of the strings so
%                                                %   they can be easily seen over
%                                                %   the background color
% set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
% 
% set(gca, 'XTick', 1:4, ...                             % Change the axes tick marks
%          'XTickLabel', {'A', 'B', 'C', 'D'}, ...  %   and tick labels
%          'YTick', 1:3, ...
%          'YTickLabel', {'A', 'B', 'C'}, ...
%          'TickLength', [0 0]);
     
