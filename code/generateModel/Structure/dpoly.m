function d=dpoly(p,pv)
    np=size(p,1);
    nvs=size(pv,1)-1;
    ds=dsegment(p,pv);
    d=min(ds,[],2);
    d=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d;
end