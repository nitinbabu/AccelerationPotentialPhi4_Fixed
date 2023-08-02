function [nearestNum, position] = findNearestNumber(array, targetNum)
    % Initialize variables
    minDifference = inf;
    position = -1;

    % Iterate over each element in the array
    for i = 1:length(array)
        % Calculate the difference between the target number and the current element
        difference = abs(targetNum - array(i));

        % Update the nearest number and its position if the current difference is smaller
        if difference < minDifference
            minDifference = difference;
            nearestNum = array(i);
            position = i;
        end
    end
end
