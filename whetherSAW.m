function result = whetherSAW(rw)
    result = ~ismember(rw(end,:),rw(1:end-1,:),'rows');
end