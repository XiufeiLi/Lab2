function next = getNextd(rw,d)
% Get next possible point fulfilling the self-avoiding requirement given a random walk
    lastPoint = rw(end,:);
    if isequal(lastPoint,zeros(1,d)) && (size(rw,1) > 1)
        next = [];
        return
    end
    possibleNext = zeros(2*d,d);
    for i = 1:d
        point = lastPoint;
        point(i) = point(i) - 1;
        possibleNext(2*i-1,:) = point;
        
        point = lastPoint;
        point(i) = point(i) + 1;
        possibleNext(2*i,:) = point;
    end
    next = [];
    for i = 1:size(possibleNext,1)
        if ~ismember(possibleNext(i,:),rw,'rows') % ~isequal(possibleNext(i,:), rw(end-1,:)) && 
            next(end+1,:) = possibleNext(i,:);
        end
    end
end