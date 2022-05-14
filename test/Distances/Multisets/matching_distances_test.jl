@testset "matching_distances" begin

    V = 1:10
    n = 30
    x = [rand(V, rand(4:7)) for i in 1:n]
    y = [rand(V, rand(4:7)) for i in 1:2n]

    dists = [
        MatchingDistance(LCS()),
        FastMatchDist(LCS(), 2n),
        FixPenMatchDist(LCS(), 1.0)

    ]
    # Match and FastMatch should be equal
    @test dists[1](x,y) == dists[2](x,y)  

    # Test properties of distances 
    for d in dists 
        @test d(x,x) == 0.0
        @test d(x,y) == d(y,x)
        @test d(nothing, nothing) == 0.0
    end 
end
