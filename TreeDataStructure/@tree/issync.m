function flag = issync(obj, anotherobj)

    if numel(obj.Parent) ~= numel(anotherobj.Parent)
        flag = false;
    else
        flag = all (  obj.Parent == anotherobj.Parent );
    end

end