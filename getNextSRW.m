function next = getNextSRW(lastPoint)
    x = lastPoint(1);
    y = lastPoint(2);
    possibleNext = [[x+1,y];[x-1,y];[x,y+1];[x,y-1]];
    next = possibleNext;
end