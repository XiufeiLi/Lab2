function next = getNext(rw)
% Get next possible point fulfilling the self-avoiding requirement given a random walk
    lastPoint = rw(end,:);
    if isequal(lastPoint,[0,0]) && (size(rw,1) > 1)
        next = [];
        return
    end
    x = lastPoint(1);
    y = lastPoint(2);
    possibleNext = [[x+1,y];[x-1,y];[x,y+1];[x,y-1]];
    next = [];
    for i = 1:size(possibleNext,1)
        if ~ismember(possibleNext(i,:),rw,'rows') % ~isequal(possibleNext(i,:), rw(end-1,:)) && 
            next(end+1,:) = possibleNext(i,:);
        end
    end
end