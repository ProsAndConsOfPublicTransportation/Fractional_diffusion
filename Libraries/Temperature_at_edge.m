function u_selected=Temperature_at_edge(u,y,pb)
    u_selected=[];
    for tt=1:size(pb.Mesh.Nodes,2)
        if pb.Mesh.Nodes(2,tt)==y
            u_selected=[u_selected;u(tt,:)];
        end
    end
    u_selected=mean(u_selected);
end