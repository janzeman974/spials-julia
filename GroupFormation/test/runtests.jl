using GroupFormation
using Random
using Statistics
using Graphs
using Test
using Random
using Logging

#include("groupformation.jl")
#import .GroupFormation would import only module name


#=
setup() consists of:
Centroid(),Agent()
patchcolor() R
getconformity() R
gethkboundary() R
-getsdiro() R
!getsigmoids()
!getcolor()
!getplace()
-facexy!
updateldistancesweights()
-opiniondistance2 R
computeidentitythresholds()
computepolarisationrepeatedly()
-computecentroidspositions()
-opiniondistance3 R
-opiniondistance R
-ashpolarisation R
--updatingcentroidsopinion
--updateagentsopiniongroup
setgroupidentities()

go consists of:
x prepareeverythingforthestep
preparingmyself!
changeopinionhk
-rollopiniondice
--computesigmoid
-rollidentitydice
updatingpatchesandglobals
x avoidingruntimeerrors
x recording_situation_and_computing_polarisation
=#
@testset "setuptests" begin
    GroupFormation.whocounter = 0
    empty!(agents)
    empty!(centroids)
    empty!(distancematrices)
    empty!(ids_and_ns_of_id_groups)
    t1::Agent=Agent()
    @test t1.who == 1
    t2::Agent=Agent()
    @test t2.who == 2
    push!(agents, t1)
    push!(agents, t2)
    t3::Centroid=Centroid()
    
    @testset "turtlestest" begin
        @test t3.who == t1.who+2
    end
    
    @testset "patchcolortest" begin
        p::Patch = Patch()
        p.pxcor = 2
        p.pycor = 2
        t1.xcor = 2.5
        t1.ycor = 2.5
        col1 = patchcolor(p)
        t2.xcor = 2.5
        t2.ycor = 2.5
        col2 = patchcolor(p)
        @test col1 != col2
    end
    
    @testset "getconformity" for i in (0.1:0.1:0.9), j in (0:2)
        data.conformitydistribution = DistributionE(j)
        data.conformitymean = i
        data.conformitystd = 0.222
        a::Float64 = getconformity()
        @test a>=0 && a<=1
    end
    
    @testset "gethkboundary" for i in (0.1:0.01:0.9), j in (0:2)
        data.boundarydistribution = DistributionE(j)
        data.boundarymean = i
        data.boundarystd = 0.222
        a::Float64 = gethkboundary()
        @test a>=0 && a<=1
    end
    
    @testset "getsdiro" for i in (0.1:0.1:0.9), j in (0:2)
        data.sdirodistribution = DistributionE(j)
        data.sdiromean = i
        data.sdirostd = 0.222
        a::Float64 = getsdiro()
        @test a>=0 && a<=1
    end
    
    @testset "opiniondistancetest" begin
        t1.opinion = rand(data.numberofopiniondimensions) .* 2 .- 1 # We set opinions...
        t1.opinion = round.(t1.opinion; digits=3)     
        t1.previousopinion = copy(t1.opinion)
        t2.opinion = copy(t1.opinion)     
        t2.previousopinion = copy(t2.opinion)
        @test opiniondistance(t1, t2) == 0
        t1.opinion = [1.0, 0.0]
        t2.opinion = [0.0, -1.0]
        data.normalizedistances = false
        @test opiniondistance(t1, t2) == sqrt(2)
        data.normalizedistances = true
        @test opiniondistance(t1, t2) == 0.5
        data.normalizedistances = true
        @test opiniondistance2(t1.opinion, t2.opinion) == 0.5
        @test opiniondistance3(t1.opinion, t2.opinion) == 0.5
        data.normalizedistances = false
        @test opiniondistance3(t1.opinion, t2.opinion) == floor(sqrt(2); digits=10)
        #@test_logs (:info,"x") opiniondistance3(t1.opinion, t2.opinion)
    end
    
    @testset "getsigmoidtest" for i=1:10
        getsigmoids!(t1)
        @test t1.opinionsigmoidxoffset >= 0.0 && t1.opinionsigmoidxoffset <= 1.0
        @test t1.opinionsigmoidsteepness >= 0.0 && t1.opinionsigmoidsteepness <= 1.0
        @test t1.identitysigmoidxoffset >= 0.0 && t1.identitysigmoidxoffset <= 1.0
        @test t1.identitysigmoidsteepness >= 0.0 && t1.identitysigmoidsteepness <= 1.0
    end
    
    @testset "getcolortest" begin
        getcolor!(t1)
        @test t1.size != 0.0
        @test t1.color >= 9 && t1.color <= 21
    end
    
    @testset "facexytest" begin
        t1.xcor = 1.0
        t1.ycor = -2.0
        facexy!(t1, 3.0, 2.0) # new position
        eps::Float64 = 0.00001
        @test t1.heading >= acos(1/sqrt(5)) - eps && t1.heading <= acos(1/sqrt(5)) + eps        
    end
    
    @testset "getplacetest" for i=1:10
        seeda = i #abs(rand(Int64))
        Random.seed!(seeda) 
        t1.opinion = rand(data.numberofopiniondimensions) .* 2 .- 1 # We set opinions...
        t1.opinion = round.(t1.opinion; digits=3)
        t1.heading = 0.0
        
        beforex = t1.xcor
        beforey = t1.ycor
        getplace!(t1)
        @test t1.xcor != beforex || t1.ycor != beforey
        if isequal(t1.ycor, beforey) == false
            @test t1.heading != 0.0
        end   
    end
    
    @testset "updateldistancetest" begin
        data.numberofopiniondimensions = 2
        t1.opinion = [1.0,0.0] 
        t2.opinion = [0.0,-1.0]
        
        ldistance = Ldistance(t1, t2)
        get!(ldistances, Pair(t1, t2), ldistance)
        updateldistancesweights()
        @test getindex(ldistances, Pair(t1, t2)).lweight == 0.5 #1-sqrt(4+4)/2sqrt(2)
        delete!(ldistances, Pair(t1, t2))
    end
    
    @testset "computeidentitythresholdstest" begin
        data.maximumsdiro = 0.89
        data.minimumsdiro = 0.29
        data.identitylevels = 7
        t1.sdiro = getsdiro()
        beforet1 = t1.sdiro 
        t2.sdiro = getsdiro()
        beforet2 = t2.sdiro
        computeidentitythresholds()
        sdirosettst =[0.29, 0.39, 0.49, 0.59, 0.69, 0.79, 0.89] 
        @test all(in.(sdiroset, Ref(sdirosettst)))
        @test t1.sdiro in sdiroset
        if in(beforet1, sdirosettst) == false 
            @test beforet1 != t1.sdiro
        end     
        @test beforet1 <= 0.19 || beforet1 + 0.1 > t1.sdiro
        @test t2.sdiro in sdiroset
        if in(beforet2, sdirosettst) == false 
            @test beforet2 != t2.sdiro
        end     
        @test beforet2 + 0.1 > t2.sdiro
    end
    
    @testset "updatingcentroidsopiniontest" begin
        t4 = Centroid()
        @test t4.who == 4
        t4.opinion = round.(ones(data.numberofopiniondimensions) - rand(data.numberofopiniondimensions) .* 2; digits=8)
        t4.who = 20
        get!(centroids, 20, t4)
        t5 = Centroid()
        @test t5.who == 5
        t5.who = 21
        t5.opinion = round.(ones(data.numberofopiniondimensions) - rand(data.numberofopiniondimensions) .* 2; digits=8)
        get!(centroids, 21, t5)
        data.numberofopiniondimensions = 2
        t6::Agent = Agent()
        push!(agents, t6)
        t1.opinion = [1, 2]
        t2.opinion = [3, 4]
        t6.opinion = [5, 6]
        t1.groupnumber = 20
        t2.groupnumber = 20
        t6.groupnumber = 21 
        updatingcentroidsopinion(20, 21)
        @test t4.opinion == [2.0, 3.0]
        @test t5.opinion == [5.0, 6.0]
        t1.groupnumber = 20
        t2.groupnumber = 21 
        t6.groupnumber = 21
        updatingcentroidsopinion(20, 21)
        @test t4.opinion == [1.0, 2.0]
        @test t5.opinion == [4.0, 5.0]
        t1.groupnumber = 20
        t2.groupnumber = 21 
        t6.groupnumber = 20
        updatingcentroidsopinion(20, 21)
        @test t4.opinion == [3.0, 4.0]
        @test t5.opinion == [3.0, 4.0]
    end
    
    @testset "updateagentsopiniongrouptest" begin
        t6::Agent = agents[3]
        t1.opinion = [1,2]
        t2.opinion = [5,6]
        t6.opinion = [5,7]
        
        t4::Centroid = getindex(centroids, 20)  
        t5::Centroid = getindex(centroids, 21)
        
        t1.groupnumber = 20
        t2.groupnumber = 20
        t6.groupnumber = 21
        
        t4.opinion = [1,2]
        t5.opinion = [5,6]

        updateagentsopiniongroup(20, 21)
        @test t1.groupnumber == 20
        @test t2.groupnumber == 20 # stayed
        @test t6.groupnumber == 21
        
        #all 3 wrong - changed will be
        t1.groupnumber = 21 # wrongly at grp0
        t2.groupnumber = 20 # wrongly at grp1
        t6.groupnumber = 20 # wrongly at grp1 -- more distant from cent0
        
        #i.e.
        updateagentsopiniongroup(20, 21)
        @test t1.groupnumber == 20
        @test t2.groupnumber == 20 # stayed
        @test t6.groupnumber == 21
        
        #4 - same number of wrongs means change all
        t1.groupnumber = 21 
        t2.groupnumber = 20 
        t6.groupnumber = 20
        t7::Agent = Agent()
        t7.opinion = [1,3]
        t7.groupnumber = 21
        push!(agents, t7)
        updateagentsopiniongroup(20, 21)
        @test t1.groupnumber == 20
        @test t2.groupnumber == 21
        @test t6.groupnumber == 21
        @test t7.groupnumber == 20
        popat!(agents, 4)
        empty!(centroids)   
    end
    
    @testset "ashpolarisationtest" begin
        t6::Agent = agents[3]
        t1.opinion = [1,1]
        t2.opinion = [1,2]
        t6.opinion = [1,5]
        t7::Agent = Agent()
        t7.opinion = [1,5]
        push!(agents, t7)
        
        data.esbgfurthestout = 1
        data.normalizedistances = true
        data.numberofagents = 4
        ap = ashpolarisation()
        eps = 0.0000001
        @test ap >= 3/(2*sqrt(2)) - eps && ap <= 3/(2*sqrt(2)) + eps
        data.normalizedistances = false
        ap = ashpolarisation()
        eps = 0.0000001
        @test ap >= 3/4 * sqrt(2) - eps && ap <= 3/4 * sqrt(2) + eps
        data.esbgfurthestout = 0
        ap = ashpolarisation()
        eps = 0.0000001
        diva = 0.5 / sqrt(8)
        @test ap >= 3.5/4 * sqrt(2) / (1 + diva) - eps && ap <= 3.5/4 * sqrt(2) / (1 + diva) + eps
    end
    
    @testset "computepolarisationrepeatedlytest" begin
        data.polarrepeats = 5
        computepolarisationrepeatedly()
        eps = 0.0000001
        diva = 0.5 / sqrt(8)
        tst = round(3.5/4 * sqrt(2) / (1 + diva); digits=3)
        @test esbgpolarisation == tst       
        popat!(agents, 4)
    end
    
    @testset "computecentroidspositionstest" begin
        #copied from updatingcentroidsopinion
        t4 = Centroid()
        t4.who = 20
        t5 = Centroid()
        t5.who = 21
        t4.opinion = [1.0, 1.0]
        t5.opinion = [2.0, 2.0]
        get!(centroids, 20, t4)
        get!(centroids, 21, t5)        
        
        t6::Agent = agents[3]
        t1.opinion = [1, 2]
        t2.opinion = [3, 4]
        t6.opinion = [5, 6]
        
        t1.groupnumber = 21
        t2.groupnumber = 21 
        t6.groupnumber = 21
        computecentroidspositions(agents)
        @test t4.opinion == [1.0, 1.0]
        @test t5.opinion == [3.0, 4.0]
        t1.groupnumber = 20
        t2.groupnumber = 20
        t6.groupnumber = 21 
        computecentroidspositions(agents)
        @test t4.opinion == [2.0, 3.0]
        @test t5.opinion == [5.0, 6.0]
        t1.groupnumber = 20
        t2.groupnumber = 21 
        t6.groupnumber = 21
        computecentroidspositions(agents)
        @test t4.opinion == [1.0, 2.0]
        @test t5.opinion == [4.0, 5.0]
        t1.groupnumber = 20
        t2.groupnumber = 21 
        t6.groupnumber = 20
        computecentroidspositions(agents)
        @test t4.opinion == [3.0, 4.0]
        @test t5.opinion == [3.0, 4.0]
    end
    
    @testset "setgroupidentitiestest" begin
        #1) identitytype global
        data.sdiromean = 0.45
        empty!(sdiroset)
        push!(sdiroset,data.sdiromean)
        # Secondly, we set 'own-SDIRO' of agents to the constant value:
        for agent in agents
            agent.sdiro = data.sdiromean
        end
        
        t6::Agent = agents[3]
        t1.opinion = [3.1, 4.1]
        t2.opinion = [2.9, 3.9]
        t6.opinion = [3, 4]
        t6.who  = 3
        t7::Agent = Agent()
        t8::Agent = Agent()
        t9::Agent = Agent()
        t7.opinion = [8,9]
        t8.opinion = [7.9,8.9]
        t9.opinion = [8.1,9.1]        
        
        t7.who  = 4
        t8.who  = 5
        t9.who  = 6
        push!(agents, t7)
        push!(agents, t8)
        push!(agents, t9)
        #=t10::Agent = Agent()
        t11::Agent = Agent()
        push!(agents, t10)
        push!(agents, t11)
        t10.opinion = [3.1,4.1]
        t11.opinion = [4.1,5.1]=#
        for i = 1 : (length(agents) - 1)
            ag0 = agents[i]
            for j = (i + 1) : length(agents)
                ag1 = agents[j]
                ldistance = Ldistance(ag0, ag1)
                get!(ldistances, ldistance.ends, ldistance)
                push!(ag0.ldistances, ldistance)
                push!(ag1.ldistances, ldistance) #JZ always the same obj.
                ldistance.lweight = -10000.0 #TOdO
            end
        end
        updateldistancesweights() # Setting distances links' weights
        for ldistance in values(ldistances) # Hiding links for saving comp. resources
            ldistance.hidden = true
        end
        
        @test isempty(distancematrices) == true
        @test isempty(ids_and_ns_of_id_groups) == true
        for agent in agents
            agent.groupnumber = 0
        end
        setgroupidentities()
        for agent in agents
            @test agent.groupnumber != 0
        end
        @test isempty(distancematrices) == false
        @test isempty(ids_and_ns_of_id_groups) == false
        @test getindex(ids_and_ns_of_id_groups, 0.45)[1] == 2
        mincentwho = minimum(agent.groupnumber for agent in agents)
        @test getindex(agents[1].which_group_has_each_sdiro_sorted_me_in, 0.45) == mincentwho
        @test getindex(agents[4].which_group_has_each_sdiro_sorted_me_in, 0.45) == mincentwho+1
        @test getindex(distancematrices, 0.45)[1,2] == round(5*sqrt(2); digits=3)
    end
    
    @testset "setupmethodtest" begin
        #data.numberofagents = 129 že stačí
        data.useidentity = false
        data.killingcentroids = false
        #cannot test @test esbgpolarisation == 0.0
        empty!(distancematrices)
        empty!(centroids)
        setup()
        @test isempty(centroids) == true
        @test isempty(distancematrices) == true
        data.useidentity = true
        data.identitytype = globala
        setup()
        @test esbgpolarisation != 0.0
        @test isempty(centroids) == false
        @test isempty(distancematrices) == false
        @test sdiroset[1] == data.sdiromean
        data.setseed = true
        data.useidentity = true
        data.identitytype = individual 
        setup()
        @test length(sdiroset) != 1
        @test esbgpolarisation != 0.0
        @test isempty(distancematrices) == false
        #Todo test vzdal. centroidu? nevim kolik jich bude / leda jestl
        # to v nule nebude haprovat
    end
    
    @testset "gotests" begin  
        @testset "prepareeverythingforthesteptest" begin
            prepareeverythingforthestep()
            @test length(agents) != 0
        end
    
        @testset "preparingmyselftest" begin
            agent::Agent = Agent()
            beforex = agent.xcor
            beforey = agent.ycor
            agent.opinion = [0.5, 0.4]
            preparingmyself!(agent)
            @test agent.size != 0.0
            @test agent.color >= 9 && t1.color <= 21
            @test t1.xcor != beforex || t1.ycor != beforey
            if isequal(t1.ycor, beforey) == false
                @test t1.heading != 0.0
            end   
        end
        
        @testset "updatingpatchesandglobalstest" begin
            empty!(agents)
            data.recordlength = 2
            push!(agents, Agent())
            push!(agents, Agent())
            push!(agents[1].record, 0.0)
            push!(agents[1].record, 0.0)
            push!(agents[2].record, 0.0)
            push!(agents[2].record, 0.0)
            agents[1].opinion = [0.5, 0.5] 
            #not the same as previousopinion prepends record 0
            updatingpatchesandglobals()
            @test mainrecord[1] == 0.25
            agents[1].record[1] = 1.0
            agents[1].record[2] = 0.0
            agents[2].record[1] = 1.0
            agents[2].record[2] = 0.0
            #not the same as previousopinion prepends record 0
            updatingpatchesandglobals()
            @test mainrecord[1] == 0.75
            @test mainrecord[2] == 0.25
            agents[1].record[1] = 1.0
            agents[1].record[2] = 0.0
            agents[2].record[1] = 1.0
            agents[2].record[2] = 0.0
            agents[1].previousopinion = [0.5, 0.5] 
            #same as previousopinion prepends record 1
            updatingpatchesandglobals()
            @test mainrecord[1] == 1.0
            @test mainrecord[2] == 0.75
        end
        
        @testset "computesigmoidtest" begin
            data.maxsteepness = 2
            val::Float64 = computesigmoid(2,1,1) # x, xoffset, steepness
            e::Float64 = 0.001
            @test val - e <= 1/(1+exp(2)) && val + e >= 1/(1+exp(2)) 
        end
        
        @testset "rolldicestest" begin
            data.opinionsigmoid == true 
            data.maxsteepness = 700 
            agents[1].opinionsigmoidxoffset = 1
            agents[1].opinionsigmoidsteepness = 1
            rollopiniondice(agents[1], agents[2], 20000000.0)
            @test agents[2].opiniondice == false
            agents[1].identitysigmoidxoffset = 1
            agents[1].identitysigmoidsteepness = 1
            rollidentitydice(agents[1], agents[2], 20000000.0)
            @test agents[2].identitydice == false
        end
        
        @testset "changeopinionhktest" begin
            agents[1].sdiro = data.sdiromean
            get!(agents[1].which_group_has_each_sdiro_sorted_me_in, data.sdiromean, 1)
            agents[2].sdiro = data.sdiromean
            get!(agents[2].which_group_has_each_sdiro_sorted_me_in, data.sdiromean , 1)
            get!(distancematrices, data.sdiromean, [-20000000.0;;])
            get!(ids_and_ns_of_id_groups, data.sdiromean, [1, 1])
            agents[1].opinion = [0.0, 0.0]
            befopinion = [0.0, 0.0]
            agents[2].opinion = [20000000.0 , 0.0]
            #agents[2] will not be an influencer means the same opinion
            changeopinionhk!(agents[1])
            @test agents[1].opinion == befopinion
            agents[2].conformity = data.conformitymean
            agents[2].opinion = [0.0 , 0.0]
            changeopinionhk!(agents[2])
            @test agents[2].opinion == befopinion
            data.updating = 2 # both opinion dimensions updating
            agents[2].identitysigmoidsteepness = 1.0
            agents[2].opinionsigmoidsteepness = 1.0 #all others will be influentials
            agents[2].opinionsigmoidxoffset = 2000000.0 #all others will be influentials
            agents[1].opinion = [1.0 , 1.0]
            changeopinionhk!(agents[2])
            @test agents[2].opinion != befopinion
        end
    end
end