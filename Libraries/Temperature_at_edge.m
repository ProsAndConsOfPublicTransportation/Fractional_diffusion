function u_selected=Temperature_at_edge(u,y,x,pb)
    u_selected=[];
    x_min=x;
    if x~=-1
        for tt=1:size(pb.Mesh.Nodes,2)
            if (pb.Mesh.Nodes(2,tt)==y) && (x_min>abs(pb.Mesh.Nodes(1,tt)-x))
                u_selected=u(tt,:);
                x_min=abs(pb.Mesh.Nodes(1,tt)-x);
            end
        end
    else
        for tt=1:size(pb.Mesh.Nodes,2)
            if (pb.Mesh.Nodes(2,tt)==y)
                u_selected=[u_selected;u(tt,:)];
            end
        end
        u_selected=mean(u_selected);
    end
end