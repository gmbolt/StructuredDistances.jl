@testset "matching_distances" begin

    V = 1:10
    x = [rand(V, rand(4:7)) for i in 1:5]
    y = [rand(V, rand(4:7)) for i in 1:10]

    dists = [
        MatchingDistance(LCS()),
        FastMatchDist(LCS(), 100),
        FixPenMatchDist(LCS(), 1.0)

    ]
    for d in dists 
        @test d(x,y) == d(y,x)
    end 
end
