using FLoops, Distances

export threaded_pairwise, non_threaded_pairwise, threaded_inbounds_pairwise


function non_threaded_pairwise(d::Metric,x,y)
    C = zeros(length(x),length(y))
    for i in 1:size(C,1)
        for j in 1:size(C,2)
            C[i,j] = d(x[i],y[j])
        end 
    end 
    return C 
end 

function threaded_pairwise(d::Metric, x, y)
    C = zeros(length(x),length(y))
    @floop for (i,j) in Iterators.product(eachindex(x), eachindex(y))
            C[i,j] = d(x[i],y[j])
    end
    
    return C 
end 

function threaded_inbounds_pairwise(d::Metric, x, y)
    C = zeros(length(x),length(y))
    @inbounds begin
        @floop for (i,j) in Iterators.product(eachindex(x), eachindex(y))
            C[i,j] = d(x[i],y[j])
        end
    end  
    return C 
end 
