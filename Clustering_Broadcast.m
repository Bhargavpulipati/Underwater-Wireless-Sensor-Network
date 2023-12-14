function [z, out] = Clustering_Broadcast(m, X)

    % Calculate Distance Matrix
    d = pdist2(X, m);
    
    % Assign Clusters and Find Closest Distances
    [dmin, ind] = min(d, [], 2);
    
    % Transmission cost
    Tc = sum(dmin);
    
    z=Tc;

    out.d=d;
    out.dmin=dmin;
    out.ind=ind;
    out.SelThresh=Tc;
    
end