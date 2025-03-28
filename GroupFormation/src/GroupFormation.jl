module GroupFormation

using Random
using Statistics
using Graphs
using SimpleWeightedGraphs
using Logging

include("CommunityDetection.jl")
using .CommunityDetection

macro exportinstances(enum)
    eval = GlobalRef(Core, :eval)
    return :($eval($__module__, Expr(:export, map(Symbol, instances($enum))...)))
end
@enum DistributionE constant uniform normal
@enum IdentityE individual globala nonidentity
@exportinstances IdentityE
@enum ShapeE circle square
@enum ModelE hk

import Random:
    rand, seed! 

import Base:
    Vector, Set, Pair, ==, copy, length, zeros, collect, ones, in, round, rand
    
import Statistics:
    mean

export
    Agent, Turtle, Centroid, Ldistance, Patch, DistributionE, IdentityE, agents, ldistances,
    centroids, data, sdiroset, esbgpolarisation, mainrecord, whocounter, opiniondistance, 
    opiniondistance2, opiniondistance3, patchcolor, getconformity, getsdiro, 
    gethkboundary, getsigmoids!, getcolor!, getplace!, facexy!, 



minpxcor::Int64 = -16
minpycor::Int64 = -16
maxpxcor::Int64 = 16 
maxpycor::Int64 = 16 
ticks::Int64 = 1

mutable struct Patch
    pxcor::Int64
    pycor::Int64 
    pcolor::Float64
    
    Patch(pxcor::Int64 = 0, pycor::Int64 = 0, 
        pcolor::Float64 = 0.0) = new(pxcor, pycor, pcolor)
end
Patch((a,b)::Tuple{Int64, Int64}) = Patch(a,b)
patches::Matrix{Patch} = map(Patch, (i,j) 
    for i=minpxcor:maxpxcor, j=minpycor:maxpycor)





updateldistancesweights, 
    computeidentitythresholds, updateagentsopiniongroup, updatingcentroidsopinion,
    ashpolarisation, computepolarisationrepeatedly, setgroupidentities, 
    computecentroidspositions, ids_and_ns_of_id_groups, distancematrices,
    setup, preparingmyself!, prepareeverythingforthestep, updatingpatchesandglobals,
    computesigmoid, rollopiniondice, rollidentitydice, changeopinionhk!


abstract type Turtle end

mutable struct Data
    rs::Int64
	xopinion::Int64
	yopinion::Int64
	updating::Int64
	numberofopiniondimensions::Int64
	identitysigmoid::Bool
	
	centroidcolor::Bool
	showdicerolls::Bool
	opinionsigmoid::Bool
	maxsteepness::Float64
	numberofagents::Int64
	ncentroids::Int64
	avoidseedcontrol::Bool
	setseed::Bool
	useidentity::Bool
	normalizedistances::Bool
	centroidschange::Float64
	recordlength::Int64
	identitylevels::Int64
	maximumsdiro::Float64
	minimumsdiro::Float64
	boundarydistribution::DistributionE
	boundarymean::Float64
	boundarystd::Float64
	conformitydistribution::DistributionE
	sdirodistribution::DistributionE
	sdirostd::Float64
	identitytype::IdentityE
	sdiromean::Float64
	conformitystd::Float64
	conformitymean::Float64
	esbgfurthestout::Int64
	killingcentroids::Bool
	mean_opinion_sigmoid_steepness::Float64
	std_opinion_sigmoid_steepness::Float64
	mean_identity_sigmoid_xoffset::Float64
	std_identity_sigmoid_xoffset::Float64
	mean_identity_sigmoid_steepness::Float64
	std_identity_sigmoid_steepness::Float64
	polarrepeats::Int64
	polarisation_each_n_steps::Int64
	
	maxticks::Int64
	
    #Data(rs::Int64=1235) = new(rs)
	#Data() = (x=new();x.rs=123;return x)
    #=test
        data.rs = 1.648985926E9 #inputbox
        data.setseed = true
        data.numberofagents = 129 # 10-1000
        data.numberofopiniondimensions = 2 # 1-50
        data.boundarymean = 0.2 # 0-1
        data.boundarydistribution = uniform
        data.xopinion = 1 # 1-50
        data.yopinion = 1 # 1-50
        data.updating = 1 # 1-50
        data.recordlength = 15 # 10-100
        data.maxticks = 730 # 100-10000
        data.conformitymean = 0.56 # 0-1
        data.ncentroids = 2
        data.centroidschange = 0.00001 # 0.000000001-0.001
        data.centroidcolor = true
        data.killingcentroids = true
        data.sdiromean = 0.44 # 0.01-1
        data.polarisation_each_n_steps = 400 # 0-10000
        data.polarrepeats = 50 # 1-100
        data.useidentity = true
        data.esbgfurthestout = 5 # 0-100
        data.sdirodistribution = uniform
        data.identitytype = individual
        data.identitylevels = 7 #1-10
        data.sdirostd = 0.075 #0-1
        data.conformitystd = 0.222 #0-1
        data.boundarystd = 0.05 #0-1
        data.mean_opinion_sigmoid_steepness = 0.5 #0-1
        data.maxsteepness = 700 #0-700
        data.std_opinion_sigmoid_steepness = 0.1 #0-1
        data.mean_identity_sigmoid_steepness = 0.7 #0-1
        data.std_identity_sigmoid_steepness = 0.125 #0-1
        data.minimumsdiro = 0.29 #0-0.5
        data.maximumsdiro = 0.89 #0.55-1
        data.normalizedistances = true
        data.showdicerolls = true
        data.mean_identity_sigmoid_xoffset = 0.15 #0-1
        data.std_identity_sigmoid_xoffset = 0.1 #0-1
        data.opinionsigmoid = true 
        data.identitysigmoid = true 
        data.avoidseedcontrol = true
        =#
	Data(rs::Int64=Int(1.648985926E9), xopinion::Int64=1, yopinion::Int64=1, updating::Int64=1,numberofopiniondimensions::Int64=2, identitysigmoid::Bool=true,
        centroidcolor::Bool=true, 
        showdicerolls::Bool=true, opinionsigmoid::Bool=true, maxsteepness::Float64=30.0,
        numberofagents::Int64=129, ncentroids::Int64=2, avoidseedcontrol::Bool=false,
        setseed::Bool=false, useidentity::Bool=true, normalizedistances::Bool=true,
        centroidschange::Float64=0.00001, recordlength::Int64=15, identitylevels::Int64=7,
        maximumsdiro::Float64=0.89, minimumsdiro::Float64=0.29, boundarydistribution::DistributionE=uniform,
        boundarymean::Float64=0.2, boundarystd::Float64=0.05, conformitydistribution::DistributionE=uniform,
        sdirodistribution::DistributionE=uniform, sdirostd::Float64=0.075, identitytype::IdentityE=individual,
        sdiromean::Float64=0.44, conformitystd::Float64=0.222, conformitymean::Float64=0.56,
        esbgfurthestout::Int64=5, killingcentroids::Bool=true, mean_opinion_sigmoid_steepness::Float64=0.5,
        std_opinion_sigmoid_steepness::Float64=0.1, mean_identity_sigmoid_xoffset::Float64=0.15, 
        std_identity_sigmoid_xoffset::Float64=0.1,
    mean_identity_sigmoid_steepness::Float64=0.7, std_identity_sigmoid_steepness::Float64=0.125, polarrepeats::Int64=50,
        polarisation_each_n_steps::Int64=400, maxticks::Int64=730) = new(rs, xopinion, yopinion,      
            updating, numberofopiniondimensions,
            identitysigmoid, centroidcolor, 
            showdicerolls, opinionsigmoid, maxsteepness, numberofagents,
            ncentroids, avoidseedcontrol, setseed, useidentity, normalizedistances,
            centroidschange, recordlength, identitylevels, maximumsdiro,
            minimumsdiro, boundarydistribution, boundarymean, boundarystd,
            conformitydistribution, sdirodistribution, sdirostd, identitytype, sdiromean,
            
            conformitystd, conformitymean, esbgfurthestout, killingcentroids,
            mean_opinion_sigmoid_steepness, std_opinion_sigmoid_steepness, mean_identity_sigmoid_xoffset,
            std_identity_sigmoid_xoffset, mean_identity_sigmoid_steepness, std_identity_sigmoid_steepness,
            polarrepeats, polarisation_each_n_steps, maxticks)
end
data::Data = Data()

model::ModelE = hk

mutable struct Ldistance
    ends::Pair{Turtle, Turtle}
    lweight::Float64
    hidden::Bool
    
    Ldistance(ends::Pair{Turtle, Turtle}, 
        lweight::Float64, hidden::Bool) = new(ends, lweight, hidden)
end
Ldistance(end1::Turtle, end2::Turtle) = Ldistance(Pair{Turtle,Turtle}(end1, end2), 0.0, false)
ldistances::Dict{Pair{Turtle, Turtle}, Ldistance} = Dict()

whocounter::Int64 = 0 
mutable struct Agent <: Turtle
    who::Int64
    color::Float64
    size::Float64
    xcor::Float64
    ycor::Float64
    
    opinion::Vector{Float64} # the number of dimensions will be determined by numberofopiniondimensions 

	which_group_has_each_sdiro_sorted_me_in::Dict{Float64, Int64}
	ldistances::Vector{Ldistance}
	sdiro::Float64
	identitydice::Bool
	opiniondice::Bool
	distancetocentroid::Float64
	opinionsigmoidxoffset::Float64
	identitysigmoidxoffset::Float64
	identitysigmoidsteepness::Float64
	opinionsigmoidsteepness::Float64
	boundary::Float64
	previousopinion::Vector{Float64}
	record::Vector{Float64}
	conformity::Float64
	groupnumber::Int64
	
	heading::Float64
	
	Agent(who::Int64=(global whocounter = whocounter + 1), color::Float64=0.5, size::Float64=1.0, xcor::Float64=25.0, ycor::Float64=25.0,
       opinion::Vector{Float64}=Vector{Float64}(), which_group_has_each_sdiro_sorted_me_in::Dict{Float64, Int64}=Dict{Float64, Int64}(),
	   ldistances::Vector{Ldistance}=Vector{Ldistance}(), sdiro::Float64=0.0, identitydice::Bool=true, opiniondice::Bool=true,
	   distancetocentroid::Float64=0.0, opinionsigmoidxoffset::Float64=0.0, identitysigmoidxoffset::Float64=0.0,
	   identitysigmoidsteepness::Float64=0.0, opinionsigmoidsteepness::Float64=0.0, boundary::Float64=0.0,
	   previousopinion::Vector{Float64}=Vector{Float64}(), record::Vector{Float64}=Vector{Float64}(), conformity::Float64=0.0, groupnumber::Int64=0,
	   heading::Float64=0.0) = new(who, color, size, xcor, ycor, opinion, which_group_has_each_sdiro_sorted_me_in,
	   ldistances, sdiro, identitydice, opiniondice, distancetocentroid, opinionsigmoidxoffset, 
       identitysigmoidxoffset, identitysigmoidsteepness, opinionsigmoidsteepness, boundary,
	   previousopinion, record, conformity, groupnumber, heading) 
end
agents::Vector{Agent} = []

mutable struct Centroid <: Turtle
    who::Int64
    color::Float64
    size::Float64
    xcor::Float64
    ycor::Float64
    
    opinion::Vector{Float64} # the number of dimensions will be determined by numberofopiniondimensions
    previousopinion::Vector{Float64}
    shape::ShapeE
    heading::Float64
    
    Centroid(who::Int64=(global whocounter = whocounter + 1),
        color::Float64=0.0, size::Float64=0.0, xcor::Float64=25.0,
        ycor::Float64=25.0, opinion::Vector{Float64}=Float64[],
        previousopinion::Vector{Float64}=Float64[],shape::ShapeE=circle,heading::Float64=0.0) = new(who, color, size,
        xcor, ycor, opinion, previousopinion, shape, heading)  
end
centroids::Dict{Int64, Centroid} = Dict()

mainrecord::Vector{Float64} = Float64[]
esbgpolarisation::Float64 = 0.0
sdiroset::Vector{Float64} = Float64[] 
ids_and_ns_of_id_groups::Dict{Float64, Vector{Int64}} = Dict()
distancematrices::Dict{Float64, Matrix{Float64}} = Dict()
minpxcor::Int64 = -16
minpycor::Int64 = -16
maxpxcor::Int64 = 16 
maxpycor::Int64 = 16 
ticks::Int64 = 1

mutable struct Patch
    pxcor::Int64
    pycor::Int64 
    pcolor::Float64
    
    Patch(pxcor::Int64 = 0, pycor::Int64 = 0, 
        pcolor::Float64 = 0.0) = new(pxcor, pycor, pcolor)
end
Patch((a,b)::Tuple{Int64, Int64}) = Patch(a,b)
patches::Matrix{Patch} = map(Patch, (i,j) 
    for i=minpxcor:maxpxcor, j=minpycor:maxpycor)


function facexy!(turtle::Turtle, facex::Float64, facey::Float64)
    turtle.heading = acos(
       (facex - turtle.xcor) 
       / sqrt((facex - turtle.xcor)^2 + (facey - turtle.ycor)^2))
    if (facey - turtle.ycor) < 0
        turtle.heading = 2*pi - turtle.heading
    end
end

# We compute polarisation several times and then set it for the average
function computepolarisationrepeatedly()
    # Initialization of temporal variables
    r::Int64 = 0
    ap::Vector{Float64} = []
    updateldistancesweights()
    
    # Repeating cycle
    while r < data.polarrepeats
        push!(ap, ashpolarisation())
        r += 1
    end
    
    # Setting variables back
    global esbgpolarisation = round(mean(ap); digits=3)
end

# Sub-routine for updating Distances links' weights,
# according opinion distance of both their ends
function updateldistancesweights()
    # We use function 'opinion-distance2', which needs two opinion positions as input and
    # receives their distance as output, but this distance is converted to weight:
    # weight = 1 means that both positions are same, weight = 0 means that their distance is maximal,
    # i.e. both positions are in oposit corners of respective N-dimensional space.
    for lens in keys(ldistances)
        ldistance::Ldistance = getindex(ldistances, lens)
        larr::Pair{Turtle,Turtle} = ldistance.ends
        ldistance.lweight = opiniondistance2(larr[1].opinion, larr[2].opinion)
        ldistance.hidden = true
    end
end

# Initialization and setup
function setup()
    # We erase the world and clean patches
    empty!(agents)
    global whocounter = 0
    empty!(centroids)
    #global centroidwhocounter = 0
    empty!(ldistances)
    
    for patch in patches
        patch.pcolor = patchcolor(patch)
    end
    
    # To avoid some random artificialities we have to set random seed.
    # But for BehaviorSearch we have to avoid seed control, since BeahaviorSearch does it.
    seeda::Int64 = 0
    if data.avoidseedcontrol == false # Note: if switch 'avoid_seed_control?' is TRUE, then we avoid setting the seed.
        if data.setseed == true
            seeda = data.rs
        else
            seeda = abs(rand(Int64))
        end
        seed!(seeda) 
        if data.setseed == false
            data.rs = seeda
        end
    end
    
    # Then we migh initialize agents/turtles
    for i = 1 : data.numberofagents
        a::Agent = Agent()
        push!(agents, a)
    end
    
    for agent in agents
        agent.opinion = rand(data.numberofopiniondimensions) .* 2 .- 1 # We set opinions...
        agent.opinion = round.(agent.opinion; digits=3)     
        agent.previousopinion = copy(agent.opinion) #  ...set last opinion as present opinion...
        # ... we prepare indicator of turtle's stability, at all 
        # positions we set 0 as non-stability...
        agent.record = zeros(data.recordlength)
        agent.conformity = getconformity()  # setting individual conformity level, and ...
        agent.boundary = gethkboundary()    # ... setting value of HK boundary.
        agent.sdiro = getsdiro() # Individual sensitivity for group tightness/threshold.
        getsigmoids!(agent) #;;getting sigmoid parameters for opinion and identity influence probabilities
        getcolor!(agent) # Coloring the agents according their opinion.
        getplace!(agent) # Moving agents to the opinion space according their opinions.
    end
    # We also have to set up some globals and create links for measuring opinion distances:
    # JZ not necessary in julia
    # agents = copy(turtles) # Note: If we just write 'agents = turtles', then variable 'agents' 
    # is a synonym for 'turtles', so it will contain in the future created centroids!
    ag0::Agent=agents[1]
    ag1::Agent=agents[1]
    # In case we don't use identity, we don't need to create links
    if data.useidentity == true
        #  Creating full network for computing groups and polarisation
        for i = 1 : (length(agents) - 1)
            ag0 = agents[i]
            for j = (i + 1) : length(agents)
                ag1 = agents[j]
                ldistance = Ldistance(ag0, ag1)
                get!(ldistances, ldistance.ends, ldistance)
                push!(ag0.ldistances, ldistance)
                push!(ag1.ldistances, ldistance) #JZ always the same obj.
            end
        end
        
        updateldistancesweights() # Setting distances links' weights
        for ldistance in values(ldistances) # Hiding links for saving comp. resources
            ldistance.hidden = true
        end
    end 
    
    # Setting identity levels according identity scenario:
    if data.useidentity == true && data.identitytype == individual
        # If we use 'individual' perceptions of identity groups, we have to set up levels of identity sensitivity:
        if data.identitylevels == 1
            data.maximumsdiro = data.sdiromean
        end
        computeidentitythresholds()
    else
        # If we use 'global' perception, then there is only one 
        # level and each agent has same value of 'own-SDIRO':
        # Firstly, we set 'SDIRO_set' as list of one constant 
        # value: 'SDIRO_Mean':
        global sdiroset = [data.sdiromean]
        
        # Secondly, we set 'own-SDIRO' of agents to the constant value:
        for agent in agents
            agent.sdiro = data.sdiromean
        end
    end 
    
    # Setting agents' identity groups
    # Everything is prepared for equal processing in all three 
    # cases of 'non-identity', 'global' and 'individual' perception 
    # of identity groups,
    # that's why we use only one procedure here for all three 
    # scenarios. But we handle them equally:
    # we process it for all identity levels -- in case of 
    # non-identity and global, there is just one.
    if data.useidentity == true
        setgroupidentities()
    end
    # Coloring patches according the number of agents/turtles on them.
    for patch in patches
        patch.pcolor = patchcolor(patch)
    end
    # Setting the indicator of change for the whole simulation, again as non-stable.
    global mainrecord = zeros(data.recordlength)

    # Compute polarisation
    computepolarisationrepeatedly()

    # resetticks()
    global ticks = 1
end


# for computing thresholds and assigning levels to each agent as per drawing parameters
function computeidentitythresholds()
    # Computing values:
    # We know that SDIRO_Mean lower than 0.4 produces one identity group, 
    # that's why the most tolerant group/level will have threshold 0.4. 
    # We also know that thresholds beyond 0.80 fracture the public to many groups and it makes
    # no sense use higher threshold than 0.8, so that's why the most intolerant group/level will have threshold 0.8.
    # If we will use just two identity levels, they will be 0.4 and 0.8. If we will use more levels, we will smoothly
    # distribute them between 0.4 and 0.8.
    global sdiroset = fill(data.maximumsdiro, data.identitylevels)
    for i = 1 : data.identitylevels 
        sdiroset[i] = round(data.minimumsdiro + ((i-1) * ((data.maximumsdiro - data.minimumsdiro) / (data.identitylevels - 1))); digits=3) 
    end
    
    # We "ceiling" values 'own-SDIRO' of agents to the values of 'SDIRO_set', i.e.
    # we find the closest higher value of 'SDIRO_set' to the 'own-SDIRO' of agent and
    # change it for the closest value of 'SDIRO_set'.
    for agent in agents
        # Firstly, we have to process the 'SDIRO_set': substract it from 'own-SDIRO' of agent.
        diff::Vector{Float64} = Float64[] 
        for sdiro in sdiroset
            push!(diff, abs(sdiro - agent.sdiro))
        end
        
        # Secondly, we have to find the position of minimal 'diff':
        pos::Int64 = 0  
        for i = 1 : length(diff)        
            if diff[i] == minimum(diff)
                pos = i
                break
            end
        end
        
        # Thirdly, we have to check, whether the closest value of 'SDIRO_set' is higher,
        # if not, we have to point agent to the higher value, i.e. 'set pos pos  + 1",
        # and also check whether the 'pos' points on correct items on 'SDIRO_set', i.e.
        # whether we do/not jump out of range. In case we point 'outside the list',
        # we have to set 'pos' to maximal value.
        if !(pos == length(sdiroset) || agent.sdiro < sdiroset[pos])
            pos = pos + 1 
        end
        
        # Finally, we set 'own-SDIRO' as the closest value of 'SDIRO_set':
        agent.sdiro = sdiroset[pos]
    end
    
    # Very lastly, we have to check whether all values of 'SDIRO_set' are represented in the agents,
    # it is possible, that for some low SD of random normal distributions of 'SDIRO_Mean' some values of
    # 'SDIRO_set' would not be represented by any agent.
    # Then it is obsolete to perform Louvain for non-represented values of 'SDIRO_set'. So...
    # Now we go through 'SDIRO_set' list value by value and only represented remain.
    #for i = length(sdiroset) - 1 : 1
    #   if sdiroset[i] not in agents.sdiro
    #       remove(sdiroset, sdiroset[i]) 
    #   end
    #end
    
    empty!(sdiroset)
    for agent in agents
        if in(agent.sdiro, sdiroset) == false
            push!(sdiroset, agent.sdiro)
       end
    end
    sort!(sdiroset)
end

# Setting identity groups via threshold levels to account for differences in sensitivity to group relationships
function setgroupidentities()
    # Firstly, we have to erase lists 'which_group_has_each_sdiro_sorted_me_in' of agents,
    # since every step we have to fill it by new IDs of identity centroids:
    for agent in agents
        empty!(agent.which_group_has_each_sdiro_sorted_me_in) 
    end
    # We erase/create these global tables, as well.
    empty!(ids_and_ns_of_id_groups) # NOTE: the key is the level of group sensitivity, 
        # stored is the list: the first value in the list is number of groups, second/last 
        # is the lowest ID of group centroid.
    empty!(distancematrices)

    # Secondly, we go one effective ID threshold level after another, perform Louvain and 
    # k-means clusters for every level and store memberships.
    for idtl in sdiroset
        # Cleaning environment
        empty!(centroids)

        #;; Detection of clusters via Louvain: Detection itself
        cntgivenldistances::Int64 = 0
        selectedagents::Vector{Agent} = []
        for agent in agents
            cntgivenldistances = 0
            for ldistance in agent.ldistances
                if ldistance.lweight >= idtl
                    cntgivenldistances += 1
                end
            end
            if 2 <= cntgivenldistances
                #println("pushing $(agent.who)")
                push!(selectedagents, agent)
            end
        end
        # Note: We take into account only not loosely connected agents
        selectedldistances::Vector{Ldistance} = []
        for ldistance in values(ldistances)
            if ldistance.lweight >= idtl  && 
                ldistance.ends[1] in selectedagents && 
                ldistance.ends[2] in selectedagents # necessary / otherwise empty second end possible
                    push!(selectedldistances, ldistance)        
            end
        end

        # JZ there could be graph of nodes with no edges
        # these arrays would not be accessed, because the particular agents
        # would have different sdiro
        if isempty(selectedagents) || isempty(selectedldistances)
            for agent in agents
                get!(agent.which_group_has_each_sdiro_sorted_me_in, idtl, 0) 
            end
            data.ncentroids = 0
            get!(ids_and_ns_of_id_groups, idtl, [0, 0])
            get!(distancematrices, idtl, Matrix{Float64}([;;]))
            continue
        end

        adj = zeros(Int64, length(selectedagents), length(selectedagents))
        ag1::Agent=agents[1]
        ag2::Agent=agents[1]
        for ldistance in selectedldistances
           for i = 1 : length(selectedagents)
               ag1 = selectedagents[i]
               if ag1 in ldistance.ends
                   for j = (i + 1) : length(selectedagents)
                       ag2 = selectedagents[j]
                       if ag2 in ldistance.ends
                           adj[i, j] = 1
                           adj[j, i] = 1
                           break
                       end
                   end
               end
           end
        end
                       
        g = SimpleGraph(adj)
        #;; For starting centroids we take into account only 
        #not loosely connected agents, but later we set groups 
        #for all.  
        
        #JZ whole louvain in julia
        c = community_detection_louvain(g)
        data.ncentroids = maximum(c)
        
        # Computing clusters' mean 'opinion'
        positionsclusters::Vector{Vector{Float64}} = [] # List with all positions of all clusters
        agentmembers::Vector{Agent} = []
      
        for i = 1 : maximum(c)
            for j = 1 : length(c) # todo more efficient
                if c[j] == i
                    push!(agentmembers, selectedagents[j])
                end    
            end
            one::Vector{Float64} = []  # List for one position of one cluster
            mean = 0.0
            for o=1:data.numberofopiniondimensions
                for agent in agentmembers
                    mean += agent.opinion[o]
                end
                if length(agentmembers) != 0
                    mean /= length(agentmembers)
                end
                push!(one, round(mean; digits=3))
            end
            push!(positionsclusters, one)
        end

        # Preparation of centroids -- feedeing them with communities
        #is zero at beginning
        totalminwho::Int64 = -1
        if isempty(centroids) == false
            totalminwho = minimum(keys(centroids)) - 1        
        end
        for i = 1: data.ncentroids
            centroid::Centroid = Centroid()
            if totalminwho == -1
                totalminwho = centroid.who - 1
            end
            centroid.heading = centroid.who - totalminwho
            centroid.opinion = positionsclusters[Int(centroid.heading)] # We set opinions, we try to do it smoothly...
            centroid.shape = circle
            centroid.size = 1.5
            centroid.color = 5 + (centroid.who - totalminwho) * 10 
            getplace!(centroid)
            get!(centroids, centroid.who, centroid)
        end
        
        mindist::Float64 = 1000000000000.0 
        minwho::Int64 = 0
        dist::Float64 = 0.0
    
        # Assignment of agents to groups
        for agent in agents
            mindist = 1000000000000.0 
            minwho = 0
            for centroid in values(centroids)
                dist = opiniondistance(centroid, agent)
                if dist < mindist
                    mindist = dist
                    minwho = centroid.who
                end
            end 
            agent.groupnumber = minwho
            # Sic! Here we intentionally use all agents, including loosely connected.
        end
        
        
        # Computation of centroids possitions
        computecentroidspositions(agents)
        
        # Iterating cycle -- looking for good match of centroids
        sum::Float64 = 0.0
        for centroid in values(centroids)
            sum += opiniondistance3(centroid.previousopinion, centroid.opinion)
        end
        while sum > data.centroidschange
            for agent in agents
                mindist = 1000000.0
                minwho = 0
                for centroid in values(centroids)
                    dist = opiniondistance(centroid, agent) 
                    if dist < mindist
                        mindist = dist
                        minwho = centroid.who
                    end
                end
                agent.groupnumber = minwho
            end
            
            # Computation of centroids possitions
            computecentroidspositions(agents)
            
            sum = 0.0
            for centroid in values(centroids)
                sum += opiniondistance3(centroid.previousopinion, centroid.opinion)
            end
        end
        
        # Storing ID of in-groups for the present identity level
        for agent in agents
            get!(agent.which_group_has_each_sdiro_sorted_me_in, idtl, agent.groupnumber) 
        end
        
        
        # Killing centroids without connected agents
        #JZ problematic, see JZ below
        #=keysa::Vector{Int64} = collect(keys(centroids))
        for i = length(keysa) : -1 : 1 #in keys(centroids)#
            wom::Int64 = keysa[i]
#            wom::Int64 = centroid.who
            centroid::Centroid = getindex(centroids, wom)
            grpfound::Bool = false
            for agent in agents
                if agent.groupnumber == wom
                    grpfound = true
                    break
                end
            end
            if grpfound == false
                println("gggggggggggggggggggggggg")
                delete!(centroids, wom) 
            end
        end
        data.ncentroids = length(centroids)=#
        
        # Storing numbers and IDs of groups and centroid distances in tables
        # Firstly, number and lowest ID, it's the easiest
        totalminwho = 100000
        for who in keys(centroids)
            if who < totalminwho
                totalminwho = who
            end
        end
        totalmaxwho = 0
        for who in keys(centroids)
            if who > totalmaxwho
                totalmaxwho = who
            end
        end
        get!(ids_and_ns_of_id_groups, idtl, [data.ncentroids, totalminwho]) 
        # NOTE: the first value in the list is number of groups, 
        # second/last is the lowest ID of group centroid.
        
        # Secondly, we fill in respective item in 'distance-matrices'
        # We create distance matrix ... 
        # JZ would be wrong if killed centroid without agents would be in the middle   
        # ... so do we have to know the last centroid's ID.
        #m::Matrix{Float64} = zeros(Float64, 1, 1)  # future distance matrix -- initialized as empty row list
        m::Vector{Vector{Float64}} = []
        i::Int64 = totalminwho
        while i <= totalmaxwho
            j = totalminwho
            row::Vector{Float64} = []  # future row of distance matrix -- initialized as empty list
            while j <= totalmaxwho
                if j == i
                    push!(row, 0)
                else
                    push!(row, round(opiniondistance3(getindex(centroids, i).opinion, getindex(centroids, j).opinion); digits=3))
                end
                j += 1
            end
            push!(m,row)
            i += 1
        end
        
        dm::Matrix{Float64} = reduce(vcat, transpose.(m))
        # And then save it as the table entry...
        get!(distancematrices, idtl, dm) # matrix:from-row-list 
        
        # Final coloring and killing of centroids
        if data.centroidcolor == true
            for agent in agents
                agent.color = (5 + 10 * (agent.groupnumber - totalminwho))
            end
        end
        if data.killingcentroids == true 
            empty!(centroids)
        end
  end
  #print IDs-and-ns-of-id-groups
  #print distance-matrices
end

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;   G O !   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function recording_situation_and_computing_polarisation()
    # Recording condition:
    # 1) We reached end, e.g. number of steps specified in MAX-TICKS
    # 2) Recording and computing polarisation on the fly, e.g. we reached POLARISATION-EACH-N-STEPS
    if (ticks == data.maxticks) || (ticks / data.polarisation_each_n_steps) == floor(ticks / data.polarisation_each_n_steps) 
        computepolarisationrepeatedly()
    end
end

# Main routine
function go()
    # ;;;; Preparation part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # All preparations of agents, globals, avoiding errors etc. in one sub-routine
    prepareeverythingforthestep()
    
    # ;;;; Main part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    for agent in agents
        if model == hk 
            changeopinionhk!(agent)
        end
        # Note: Now here is only Hegselmann-Krause algorithm, but in the future we might easily employ other algorithms here!
    end
      
    # ;;;; Final part ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    # ;; Recoloring patches, agents, computing how model settled down
    updatingpatchesandglobals()
    
    # Todo repaint leda
    
    global ticks += 1
    # Finishing condition:
    # 1) We reached number of steps specified in MAX-TICKS
    recording_situation_and_computing_polarisation()
    if ticks == data.maxticks 
        return
    end
end

function test()
    for i = 1 : 10
        go()
    end
end

#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;   SUB/PROCEURES !   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function prepareeverythingforthestep()
    # Just checking and avoiding runtime errors part of code
    avoidingruntimeerrors() 
    if data.useidentity == true
        updateldistancesweights()
        setgroupidentities()
    end
    for agent in agents
        preparingmyself!(agent)
    end
end

function preparingmyself!(agent::Agent)
    # Updating color and place
    getcolor!(agent)
    getplace!(agent)
    
    # storing previous opinion position as 'own-previous-opinion'
    agent.previousopinion = copy(agent.opinion) 
end

#/////////////////////////////////////////////////////
#// sub-routine for updating opinion position of turtle 
#//according the Hegselmann-Krause (2002) model
function changeopinionhk!(agent::Agent)
    #println("pre-zmena3")
    # We define all other agents as NEIGHBS, i.e. potential 
    # INFLUENTIALS
    neighs::Vector{Agent} = Agent[]
    neighs = copy(agents)
    i::Int64 = 1
    meana = 0.0
    for i = 1 : length(agents)
        if agents[i] === agent
            break
        end 
    end
    popat!(neighs, i) 
    
    # In the first block of code we have to determine 
    # INFLUENTIALS -- the agents to whom the  updating agents 
    # listens to
    # Since rolling identity dice is computationally less 
    # demanding, we start with id dice:
    # If we use identity, we check distance of group 
    # centroids of NEIGHBS and filter them out.
    # In case we dont't use identity, them NEIGHBS are 
    # still all other agents.
    
    if data.useidentity == true
        #println(sdiroset)
        #println(bitstring(agent.sdiro))
        # Identity check -- preparation of myself from global variables:
        mygroupid = getindex(agent.which_group_has_each_sdiro_sorted_me_in, agent.sdiro)
        #println(distancematrices)
        distancematrix = getindex(distancematrices, agent.sdiro)
        groupinfol = getindex(ids_and_ns_of_id_groups, agent.sdiro)
        # NOTE: the first value in the list is number of groups,
        # second/last is the lowest ID of group centroid.
        
        myi = mygroupid - groupinfol[2]
        # 'my-i' is distance matrix column we want to use, so we 
        # subtract group ID of myself from minimal group ID in 
        # respective level
        
        # Sigmoid code, ask other agents
        for otheragent in agents 
            if otheragent === agent
                continue
            end
            # Firstly, each neighbor has to find 'her-j', i.e. 
            # distance matrix row for identity check:
            hergroupid = getindex(otheragent.which_group_has_each_sdiro_sorted_me_in, agent.sdiro)           
            herj = hergroupid - groupinfol[2]
            # 'her-j' is distance matrix row we want to use, 
            # so we subtract group ID of self from minimal group 
            # ID in respective level
            # After all the computations we finally get distance of centroids from distance matrix...
            # Secondly, we roll the identity dice...
            ourdistance = distancematrix[myi + 1, herj + 1] 
                        
            # NOTE: The commented out part below can be uncommented 
            # if one wants to speed up simulation by hard-setting 
            # probabilities of
            # extreme distance values ( 0 or 1 ) to extreme values 
            # (p = 1 or p = 0). This would lead to differences in 
            # outcome especially
            # when the sigmoid steepnesses are low.

            # ifelse our-distance > 0 and our-distance < 1 [
                rollidentitydice(agent, otheragent, ourdistance)
            # ;][
            # ifelse our-distance = 0[
            #     set identity-dice? true
            # ][
            #      set identity-dice? false
            # ]
        end
        # Now we set successful AGENTS as NEIGHBS:
        empty!(neighs)
        for otheragent in agents 
            if otheragent === agent
                continue
            end
            if otheragent.identitydice == true 
                push!(neighs, otheragent)
            end
        end
        if data.showdicerolls
            println(length(neighs))
        end
    end
    
    # NOW we roll probabilistic dice based on opinion distance from influencer - if an agent is too 
    # far they are less likely to be heard this time.
    # first find out which agent is going to be heard this time.
    for other in neighs 
        opdist = opiniondistance3(other.opinion, agent.opinion)
        # Note: We compute opinion distance ('op-dist') of 
        # NEIGHBS member and updating agent and ...
        flooredopdist = floor(opdist; digits=3)
        # ... pass it rounded to 3 digits as the 
        # argument of ROLL-OPINION-DICE 
        rollopiniondice(agent, other, flooredopdist)
    end
    
    influentials::Vector{Agent} = Vector{Agent}()
    # now only follow them, the successful rollers and set them as INFLUENTIALS...
    for other in neighs 
        if other.opiniondice == true
            push!(influentials, other)
        end
    end
    if data.showdicerolls == true
        println(length(influentials))
    end
    
    # 3) we also add the updating agent into 'influentials'
    push!(influentials, agent)
    #  we check whether there is someone else then calling/updating agent in the agent set 'influentials'
    if length(influentials) > 1
        # here we draw a list of dimensions which we will update:
        # by 'range opinions' we generate list of integers from '0' to 'opinions - 1',
        # by 'n-of updating' we randomly take 'updating' number of integers from this list
        # by 'shuffle' we randomize order of the resulting list
        alloplist::Vector{Int64} = [i for i in 1:data.numberofopiniondimensions]
        oplist::Vector{Int64} = shuffle(alloplist[rand(1:end, data.updating)])   
    
        # we initialize counter 'step'
        step::Int64 = 1
    
        j::Int64 = 0
        # we go through the while-loop 'updating' times:
        while step <= data.updating
            # we initialize/set index of updated opinion dimension according the items on the 'op-list',
            # note: since we use while-loop, we go through each item of the 'op-list', step by step, iteration by iteration.
            j = oplist[step]
        
            # ad 1: averge position computation
            meana = 0.0
            for v in influentials
                meana += v.opinion[j]
            end
            if length(influentials) != 0
                meana /= length(influentials)
            end
            val::Float64 = floor(meana; digits=3)
            # NOTE: H-K model really assumes that agent 
            # adopts immediatelly the 'consesual' position
            
            # ad 2: updating/weighting 'val' by 'Conformity' 
            # and own opinion
            my::Float64 = agent.opinion[j]
            val = my + ((val - my) * agent.conformity)
            
            # ad 3: assigning the value 'val'
            agent.opinion[j] = val
            
            # advancement of counter 'step'
            step += 1
        end
    end
end

function rollidentitydice(agent::Agent, other::Agent, ourdistance::Float64) 
    if data.identitysigmoid == true 
        probabilityofinteraction = computesigmoid(ourdistance, agent.identitysigmoidxoffset, agent.identitysigmoidsteepness)
        
        dice::Float64 = rand()
        
        if dice < probabilityofinteraction
            other.identitydice = true
        else 
            other.identitydice = false
        end
        if data.showdicerolls == true && probabilityofinteraction > 0 && probabilityofinteraction < 1
            println("Identity distance: $(ourdistance), Identity distance threshold: "  
                * "$(agent.identitysigmoidxoffset),\n Sigmoidial probability: " 
                * "$(probabilityofinteraction), Rolled dice: $(dice)" 
                * "; Result: $(other.identitydice)")
        end
    else 
        # If we don't use identity sigmoid, then we set 'identitydice' according sharp difference
        if ourdistance <= agent.identitysigmoidxoffset
            other.identitydice = true
        else 
            other.identitydice = false
        end
    end
end

function rollopiniondice(agent::Agent, other::Agent, opdist::Float64) 
    if data.opinionsigmoid == true 
        probabilityofinteraction = computesigmoid(opdist, agent.opinionsigmoidxoffset, agent.opinionsigmoidsteepness)
        dice = rand()
        if dice < probabilityofinteraction
            other.opiniondice = true
        else 
            other.opiniondice = false
        end
        
        if data.showdicerolls == true && probabilityofinteraction > 0 && probabilityofinteraction < 1
            println("Opinion distance: $(opdist), Boundary: "  
                * "$(agent.boundary),\n Sigmoidial probability: " 
                * "$(probabilityofinteraction), Rolled dice: $(dice)" 
                * "; Result: $(other.opiniondice)")
        end
    else
        #println("opsigmoid false") 
        #println("opdist $(opdist), sigmoidxoffset $(agent.opinionsigmoidxoffset)")
        # If we don't use identity sigmoid, then we set 'identitydice' according sharp difference
        if opdist <= agent.opinionsigmoidxoffset
            other.opiniondice = true
        else 
            other.opiniondice = false
        end
    end
end

# subroutine for computing inverted (y decreases as x increases) sigmoid of input
function computesigmoid(x, xoffset, steepness)::Float64
  # Slider 'maxSteepness' determines maximum possible steepness of the sigmoids.
  # Sigmoid outputs from 1 to 0 as x increases from 0 to 1
  return round(1.0 / (1.0 + exp(steepness * data.maxsteepness * (x - xoffset))); digits=3)

end

# sub-routine for computing opinion distance of two comparing agents
function opiniondistance(above::Turtle, aboveabove::Turtle)::Float64
    # we store in temporary variable the opinion of the called and compared agent
    my::Vector{Float64} = above.opinion
    # we store in temporary variable the opinion of the calling and comparing agent
    her::Vector{Float64} = aboveabove.opinion
    
    # we initialize counter of step of comparison -- we will compare as many times 
    # as we have dimensions
    step::Int64 = 1
    # we initialize container where we will store squared distance in each dimensio
    dist::Float64 = 0.0
    # while loop going through each dimension, computiong distance in each 
    # dimension, squarring it and adding in the container
    while step <= data.numberofopiniondimensions
        # computiong distance in each dimension, squarring it and adding in the container
        dist = dist + (my[step] - her[step])^2
        # advancing 'step' counter by 1
        step += 1
    end
    # computing square-root of the container 'dist' -- computing 
    # Euclidean distance -- and setting it as 'dist'
    dist = sqrt(dist)
    # Computing normalization constant according switch 'normalize_sigmoid_distances?':
    normalizationa::Float64 = 0.0
    # reporting Euclidean distance
    if data.normalizedistances == true
        normalizationa = 1.0 / sqrt(4 * data.numberofopiniondimensions)
    else
        normalizationa = 1.0
    end
    return dist * normalizationa
end

# sub-routine for computing opinion distance of two comparing 
# opinion positions  -- relative distance weighted as 1 for minimal 
# distance and 0 for the maximal one
function opiniondistance2(my::Vector{Float64}, her::Vector{Float64})::Float64
    # we initialize counter of step of comparison -- we will compare 
    # as many times as we have dimensions
    step::Int64 = 1 
    # we initialize container where we will store squared distance in each dimension
    dist::Float64 = 0.0
    # while loop going through each dimension, computiong distance in each 
    # dimension, squarring it and adding in the container
    while step <= data.numberofopiniondimensions
        
        # computiong distance in each dimension, squarring it and adding in the container
        dist = dist + (my[step] - her[step])^2
        # advancing 'step' counter by 1
        step += 1
    end
    # computing square-root of the container 'dist' -- computing Euclidean distance -- and setting it as 'dist'
    dist = sqrt(dist)
    
    # Turning 'dist' into 'weight'
    weight::Float64 = (sqrt(4*data.numberofopiniondimensions) - dist) / sqrt(4*data.numberofopiniondimensions)
    
    # reporting weight of distance
    return floor(weight; digits=10)
end

# sub-routine for computing opinion distance of two comparing opinion positions  -- absolute distance without weighting
function opiniondistance3(my::Vector{Float64}, her::Vector{Float64})::Float64
    # we initialize counter of step of comparison -- we will compare as many times as we have dimensions
    step::Int64 = 1
    # we initialize container where we will store squared distance in each dimension
    dist::Float64 = 0.0
    # while loop going through each dimension, computiong distance in each 
    # dimension, squarring it and adding in the container
    while step <= data.numberofopiniondimensions
        # computiong distance in each dimension, squarring it and adding in the container
        dist = dist + (my[step] - her[step])^2
        # advancing 'step' counter by 1
        step += 1
    end
        
    # computing square-root of the container 'dist' -- computing Euclidean distance -- and 
    # setting it as 'dist'
    dist = sqrt(dist)
    # Computing normalization constant according switch 
    # 'normalize_sigmoid_distances?':
    normalizationa::Float64 = 0.0
    # reporting weight of distance
    
    if data.normalizedistances == true
        normalizationa = 1.0/sqrt(4*data.numberofopiniondimensions)
    else 
        normalizationa = 1.0
    end
    return floor(dist * normalizationa; digits=10)
end

# Updating patches and global variables
function updatingpatchesandglobals()
    # Patches update color acco;rding the number of turtles on it.
    for patch in patches 
        patch.pcolor = patchcolor(patch)
    end
    
    # We have to check here the change of opinions, resp. how many agents changed,
    # and record it for each agent and also for the whole simulation
    # Turtles update their record of changes:
    for a in agents
        # we take 1 if opinion is same, we take 0 if opinion changes, then
        # we put 1/0 on the start of the list Record, but we omit the last item from Record
        
        if a.previousopinion == a.opinion
            prepend!(a.record, 1.0) 
        else
            prepend!(a.record, 0.0)
        end
        pop!(a.record)
    end
    
    # Then we might update it for the whole:    
    meana = 0.0
    glmean::Float64 = 0.0
    for a in agents
        a.color = 0.0
        meana = mean(a.record)
        glmean += meana
    end
    if length(agents) != 0
        glmean /= length(agents)
    end
    glmean = floor(glmean; digits=3)
    prepend!(mainrecord, glmean)
    pop!(mainrecord)
    
    # Coloring agents according identity group
    if data.centroidcolor == true 
        mingroupnumber::Int64 = 100000000
        for a in agents 
            if a.groupnumber < mingroupnumber 
                mingroupnumber = a.groupnumber 
            end
        end
        for a in agents
            a.color = 5 + 10*(a.groupnumber - mingroupnumber) 
        end
    end
end

# Procedure reporting ESBG/Ashwin's polarisation
function ashpolarisation()::Float64
    cent::Centroid = Centroid()
    cent.shape = square
    cent.opinion = round.(ones(data.numberofopiniondimensions) - rand(data.numberofopiniondimensions) .* 2; digits=8)
    get!(centroids, cent.who, cent)
    
    cent = Centroid()
    cent.shape = square
    cent.opinion = round.(ones(data.numberofopiniondimensions) - rand(data.numberofopiniondimensions) .* 2; digits=8)
    get!(centroids, cent.who, cent)
    
    #;; Storing 'who' of twO new centroids
    cent1 = maximum(keys(centroids))
    cent0 = cent1 - 1
    
    # Random assignment of agents to the groups -- we create random list of 'cent0s' and 'cent1s' 
    # of same length as agents number and then assign them based on agent's WHO
    secondval::Int64 = 0
    firstval::Int64 = 0
    #setrounding not working in Julia
    
    if 0 == mod(data.numberofagents, 2)
        firstval = data.numberofagents / 2
        secondval = data.numberofagents / 2        
    else
        firstval = (data.numberofagents + 1) / 2
        secondval = (data.numberofagents - 1) / 2
    end
    
    membership::Vector{Int64} = shuffle(vcat(fill(cent0, firstval),fill(cent1, secondval)))
    for i = 1: length(agents)
        agent::Agent = agents[i]
        agent.groupnumber = membership[i]
    end
    updatingcentroidsopinion(cent0, cent1)  # Initial update
    
    # Iterating until centroids are stable
    sum::Float64 = 0.0
    for centroid in values(centroids)
        if centroid.who >= cent0
            sum += opiniondistance3(centroid.opinion, centroid.previousopinion)
        end
    end
    while data.centroidschange < sum
        updateagentsopiniongroup(cent0, cent1)
        updatingcentroidsopinion(cent0, cent1)
        sum = 0.0
        for centroid in values(centroids)
            if centroid.who >= cent0
                 sum += opiniondistance3(centroid.opinion, centroid.previousopinion)
            end
        end
    end
    
    # Computing polarisation -- cutting-out agents too distant from centroids
    a0::Vector{Agent} = Agent[]
    a1::Vector{Agent} = Agent[] 
    for agent in agents
        agent.distancetocentroid = opiniondistance(agent, getindex(centroids, agent.groupnumber))
        if agent.groupnumber == cent0
            push!(a0, agent)
        elseif agent.groupnumber == cent1
            push!(a1, agent)
        end
    end
    
    a0cutted::Vector{Agent} = Agent[]
    a0distlist::Vector{Float64} = []
    indvaluemap::Dict{Float64, Vector{Agent}} = Dict() 
    a::Agent = agents[1]
    aglist::Vector{Agent} = []
    opdista::Float64 = 0.0
    for a in a0
        opdista = opiniondistance3(a.opinion, getindex(centroids, cent0).opinion)
        push!(a0distlist, opdista)
        aglist = get(indvaluemap, opdista, [])
        push!(aglist, a)
        get!(indvaluemap, opdista, aglist)
    end
    sort!(a0distlist) 
    for i = 1 : length(a0) - data.esbgfurthestout
        doubleval::Float64 = a0distlist[i]
        aglist = getindex(indvaluemap, doubleval)
        agenttoadd::Agent = pop!(aglist)
        push!(a0cutted, agenttoadd) 
    end
    
    a1cutted::Vector{Agent} = Agent[]
    a1distlist::Vector{Float64} = []
    empty!(indvaluemap)
    for a in a1
        opdista = opiniondistance3(a.opinion, getindex(centroids, cent1).opinion)
        push!(a1distlist, opdista)
        aglist = get(indvaluemap, opdista, [])
        push!(aglist, a)
        get!(indvaluemap, opdista, aglist)
    end
    sorted = collect(keys(indvaluemap))
    sort!(a1distlist)
    for i = 1 : length(a1) - data.esbgfurthestout
        doubleval::Float64 = a1distlist[i]
        aglist = getindex(indvaluemap, doubleval)
        agenttoadd::Agent = pop!(aglist)
        push!(a1cutted, agenttoadd) 
    end
    
    # Updating centroids and agents opinion position (without furthest agents)
    meana = 0.0
    for i = 1 : data.numberofopiniondimensions
        meana = 0.0 
        for a in a0cutted
            meana += a.opinion[i]
        end
        if length(a0cutted) != 0
            meana /= length(a0cutted)
        end
        getindex(centroids,cent0).opinion[i] = floor(meana; digits=8)
    end
    getplace!(getindex(centroids, cent0))
    for a in a0cutted
        a.distancetocentroid = opiniondistance(a, getindex(centroids,cent0))
    end
    
    for i = 1 : data.numberofopiniondimensions
        meana = 0.0 
        for a in a1cutted
            meana += a.opinion[i]
        end
        if length(a1cutted) != 0
            meana /= length(a1cutted)
        end
        getindex(centroids, cent1).opinion[i] = floor(meana; digits=8)
    end
    getplace!(getindex(centroids, cent1))
    for a in a1cutted
        a.distancetocentroid = opiniondistance(a, getindex(centroids, cent1))
    end
    
    # Preparing final distances and diversity
    normalization::Float64 = 0.0 
    if data.normalizedistances == true 
        normalization = 1.0
    else 
        normalization = 1.0 / sqrt(4 * data.numberofopiniondimensions)
    end
    # NOTE: If the switch 'data.normalizedistances?' is true, then functions 
    #   'opiniondistance' and 'opiniondistance3' compute 
    # normalized distances,
    # then we must avoid normalization here, but if the switch 
    # is false, then the function computes plain distance and we 
    # must normalize here --
    # so we use same concept and almost same code as in the 
    # functions 'opiniondistance' and 'opiniondistance3' here, 
    # but reversed:
    # the true switch means no normalization here, the false switch means 
    # normalization here, so in the end we normalize values in 
    # ESBG polarization exactly once.
    
    centdist = normalization * opiniondistance3( 
        getindex(centroids,cent0).opinion, getindex(centroids,cent1).opinion)
    meana = 0.0
    for i = 1 : length(a0cutted)
        meana += a0cutted[i].distancetocentroid
    end
    if length(a0cutted) != 0
        meana /= length(a0cutted)
    end
    div0::Float64 = normalization * meana
    meana = 0.0
    for i = 1 : length(a1cutted)
        meana += a1cutted[i].distancetocentroid
    end
    if length(a1cutted) != 0
        meana /= length(a1cutted)
    end
    div1::Float64 = normalization * meana
    
    toremove = Centroid[]
    for i in keys(centroids) 
        if i >= cent0
            #push!(toremove, getindex(centroids, i))
            delete!(centroids, i)
        end
    end
    
    return centdist / (1 + div0 + div1)
end

function updateagentsopiniongroup(cent0::Int64, cent1::Int64)
    #;; Checking thoe assignment -- is the assigned centroid the nearest? If not, reassign!
    minopdist::Float64 = 10000000.0
    curopdist::Float64 = 0.0
    who::Int64 = 0
    for agent in agents
        minopdist = 10000000.0
        curopdist = 0.0
        who = 0
        for centroid in values(centroids)
            if centroid.who >= cent0
                curopdist = opiniondistance(centroid, agent)
                if curopdist < minopdist
                    minopdist = curopdist
                    who = centroid.who
                end
            end
        end
        agent.groupnumber = agent.groupnumber - who
    end
    # ;commented set color 15 + group * 10]
    wronglyatgrp0::Vector{Agent} = []
    for agent in agents
        if agent.groupnumber == -1 
            push!(wronglyatgrp0, agent)
        end
    end # they are in 0, but should be in 1: 0 - 1 = -1
    wronglyatgrp1::Vector{Agent} = []
    for agent in agents
        if agent.groupnumber == 1
            push!(wronglyatgrp1, agent)
        end
    end # they are in 1, but should be in 0: 1 - 0 = 1
    if length(wronglyatgrp0) == length(wronglyatgrp1)
        for agent in agents
            minopdist = 10000000.0
            curopdist = 0.0
            who = 0
            for centroid in values(centroids)
                if centroid.who >= cent0
                    curopdist = opiniondistance(centroid, agent)
                    if curopdist < minopdist
                        minopdist = curopdist
                        who = centroid.who
                    end
                end
            end
            agent.groupnumber = who
        end
    else
        peleton::Vector{Agent} = []
        for agent in agents
            if agent.groupnumber == 0
                push!(peleton, agent)
            end
        end
        pridano::Int64 = 0
        opdist::Float64 = 0.0
        optdistagent::Dict{Float64, Agent} = Dict()
        sorted::Vector{Float64} = []
        stayed::Vector{Agent} = []
        
        if length(wronglyatgrp0) < length(wronglyatgrp1)
            for k in wronglyatgrp0
                push!(peleton, k)
            end
            if length(wronglyatgrp0) != 0
                pridano = 0
                opdist = 0.0
                empty!(optdistagent)
                for a in wronglyatgrp1
                    opdist = opiniondistance3(a.opinion, getindex(centroids, cent0).opinion)
                    get!(optdistagent, opdist, a)
                end
                empty!(sorted)
                for k in keys(optdistagent)
                    push!(sorted, k)
                end
                sort!(sorted)
                for i::Int64 = reverse(1:length(sorted))
                    if pridano >= length(wronglyatgrp0)
                        break
                    end 
                    push!(peleton, getindex(optdistagent, sorted[i]))
                    pridano += 1
                end
            end
            # all agents assigned correctly + smaller group of wrong + from bigger group 'n of size of smaller group'
        else
            ## all agents assigned correctly + smaller group of wrong + 
                ## from bigger group 'n of size of smaller group'
            for k in wronglyatgrp1
                push!(peleton, k)
            end
            if length(wronglyatgrp1) != 0
                pridano = 0
                opdist = 0.0
                empty!(optdistagent)
                for a in wronglyatgrp0
                    opdist = opiniondistance3(a.opinion, getindex(centroids, cent1).opinion) 
                    get!(optdistagent, opdist, a)
                end
                empty!(sorted) 
                for j in keys(optdistagent)  
                    push!(sorted, j)
                end
                sort!(sorted)
                for i = length(sorted) : -1 : 1 
                    if pridano >= length(wronglyatgrp1)
                        break
                    end
                    push!(peleton, getindex(optdistagent, sorted[i]))
                    pridano += 1
                end
            end
        end
        empty!(stayed)
        for agent in agents
            if in(agent, peleton) == false
                push!(stayed, agent)
            end 
        end
                      
        for a::Agent in peleton
            minwhoopdist = 10000000
            minopdist = 1000000000.0
            opdist = 0.0
            for b::Int64 in keys(centroids) 
                if b >= cent0 
                    opdist = opiniondistance(a, getindex(centroids, b))
                    if opdist < minopdist
                        minopdist = opdist
                        minwhoopdist = b
                    end
                end
            end
            a.groupnumber = minwhoopdist
        end
        # set color 15 + group * 10
        if length(wronglyatgrp0) < length(wronglyatgrp1)
            for b::Agent in stayed
                b.groupnumber = cent1
            end
        else
            for b::Agent in stayed
                b.groupnumber = cent0
            end
        end
        #;set color 15 + group * 10
    end
end

function updatingcentroidsopinion(cent0::Int64, cent1::Int64)
    # Storing opinion as own-previous-opinion
    for centroid in values(centroids)
        if centroid.who >= cent0
            centroid.previousopinion = centroid.opinion
        end
    end
    # Computing groups mean 'own-opinion'
    positionscluster::Vector{Vector{Float64}} = []
    # List with all positions of both 2 groups
    oneposition::Vector{Float64} = []
    meana::Float64 = 0.0
    cnt::Int64 = 0
    for i::Int64 = 0 : 1 
        oneposition = []
        for j::Int64 = 1 : data.numberofopiniondimensions 
            meana = 0.0
            cnt = 0
            for a::Agent in agents 
                if a.groupnumber == cent0 + i 
                    meana = meana + a.opinion[j]
                    cnt += 1
                end
            end
            if cnt != 0
                meana /= cnt
            end
            meana = floor(meana; digits=8)
            push!(oneposition, meana)
        end
        push!(positionscluster, oneposition)
    end
    getindex(centroids, cent0).opinion = positionscluster[1]
    getplace!(getindex(centroids, cent0))
    getindex(centroids, cent1).opinion = positionscluster[2]
    getplace!(getindex(centroids, cent1))
end

# Sub-routine of polarization routine
function computecentroidspositions(selagents::Vector{Agent})
    # Preparation
    for centroid in values(centroids)
        centroid.previousopinion = copy(centroid.opinion)
    end

    # Computation of centoids positions
    grp::Int64 = 100000
    maxwho::Int64 = 0
    for who in keys(centroids)
        if who <= grp
            grp = who
        end 
        if who > maxwho
            maxwho = who
        end
    end
    dim::Int64 = 0
    meana::Float64 = 0.0
    curselagents::Vector{Agent} = []
    while grp <= maxwho
        cent::Centroid = getindex(centroids, grp)
        isagentwithgrp::Bool = false
        for agent in agents
            if agent.groupnumber == grp
                isagentwithgrp = true
                break
            end
        end
        if isagentwithgrp == false
            cent.opinion = copy(cent.previousopinion)
        else
            dim = 1
            meana = 0.0
            while dim <= data.numberofopiniondimensions
                empty!(curselagents)
                for agent in selagents
                    if agent.groupnumber == grp
                        push!(curselagents, agent)
                    end
                end
                meana = 0.0
                for agent in curselagents
                    meana += agent.opinion[dim]
                end
                if length(curselagents) != 0
                    meana /= length(curselagents)
                end
                cent.opinion[dim] = meana
                dim += 1
            end
        end
        grp += 1
    end
    for centroid in values(centroids)
        getplace!(centroid)
    end
end

# ;; Sub-routine for assigning value of conformity
function getconformity()::Float64
    # We have to initialize empty temporary variable
    cvalue::Float64 = 0.0
    # Then we draw the value according the chosen method
    if data.conformitydistribution == constant
        # ;; NOTE! 'random-float 0' is here for consuming one 
        # pseudorandom number to consume same number of 
        # pseudorandom numbers as uniform
        cvalue = data.conformitymean + rand() * 0
    elseif data.conformitydistribution == uniform
        if data.conformitymean <= 0.5
            cvalue = round(rand() * (2 * data.conformitymean); digits=3) 
        else
            cvalue = 1.0 - round(rand() * (2 * (1.0 - data.conformitymean)); digits=3)
        end
    elseif data.conformitydistribution == normal
        cvalue = round(randn() * data.conformitystd + data.conformitymean; digits=3)
        while cvalue > 1 || cvalue <= 0 
            cvalue = round(randn() * data.conformitystd + data.conformitymean; digits=3)
        end
    end
    return cvalue
end

# Sub-routine for assigning value of
function getsdiro()::Float64
    # We have to initialize empty temporary variable
    gtvalue::Float64 = 0.0
        
    # Then we draw the value according the chosen method
    if data.sdirodistribution == constant 
    # ;; NOTE! 'random-float 0' is here for consuming one 
    # pseudorandom number to cunsume same number of pseudorandom 
    # numbers as "uniform
        gtvalue = data.sdiromean + rand() * 0
    elseif data.sdirodistribution == uniform
        if data.sdiromean <= 0.5
            gtvalue = round(rand() * (2 * data.sdiromean); digits=3) 
        else
            gtvalue = 1.0 - round(rand() * (2 * (1.0 - data.sdiromean)); digits=3) 
        end
    elseif data.sdirodistribution == normal
        gtvalue = round(randn() * data.sdirostd + data.sdiromean; digits=3)
        while gtvalue > 1 || gtvalue <= 0 
            gtvalue = round(randn() * data.sdirostd + data.sdiromean; digits=3)
        end
    end
    return gtvalue
end

# sub-routine for assigning value of boundary to the agent
function gethkboundary()::Float64
    # We have to initialize empty temporary variable
    uvalue::Float64 = 0.0
    # Then we draw the value according the chosen method
    if data.boundarydistribution == constant 
        # NOTE! 'random-float 0' is here for consuming one 
        # pseudorandom number to cunsume same number of 
        # pseudorandom numbers as "uniform"
        uvalue = data.boundarymean + rand() * 0
    elseif data.boundarydistribution == uniform
        #uvalue = round(rand() * (2 * data.boundarymean); digits=3)
        if data.boundarymean <= 0.5
            uvalue = round(rand() * (2 * data.boundarymean); digits=3) 
        else
            uvalue = 1.0 - round(rand() * (2 * (1.0 - data.boundarymean)); digits=3) 
        end 
    elseif data.boundarydistribution == normal
        uvalue = round(randn() * data.boundarystd + data.boundarymean; digits=3)
        while uvalue > 1 || uvalue <= 0 
            uvalue = round(randn() * data.boundarystd + data.boundarymean; digits=3)
        end
    end
    return uvalue
end


# Sub procedure for generating needed parameters for sigmoids
function getsigmoids!(agent::Agent) # JZ meaning setsigmoids
    # ;; Opinion ones
    agent.opinionsigmoidxoffset = agent.boundary
    agent.opinionsigmoidsteepness = round( 
        randn() * data.std_opinion_sigmoid_steepness + data.mean_opinion_sigmoid_steepness; digits=3)
    while agent.opinionsigmoidsteepness > 1 ||  agent.opinionsigmoidsteepness < 0 
        agent.opinionsigmoidsteepness = round( 
            randn() * data.std_opinion_sigmoid_steepness + data.mean_opinion_sigmoid_steepness; digits=3)
    end
    # Identity ones
    agent.identitysigmoidxoffset = round(randn() * data.std_identity_sigmoid_xoffset + data.mean_identity_sigmoid_xoffset; digits=3)
    while agent.identitysigmoidxoffset > 1 || agent.identitysigmoidxoffset < 0 
        agent.identitysigmoidxoffset = round(randn() * data.std_identity_sigmoid_xoffset + data.mean_identity_sigmoid_xoffset; digits=3) 
    end
    agent.identitysigmoidsteepness = round(randn() * data.std_identity_sigmoid_steepness + data.mean_identity_sigmoid_steepness; digits=3)
    while agent.identitysigmoidsteepness > 1    ||  agent.identitysigmoidsteepness < 0 
        agent.identitysigmoidsteepness = round(randn() * data.std_identity_sigmoid_steepness + data.mean_identity_sigmoid_steepness; digits=3)
    end
end

# sub-routine for graphical representation -- 
# it takes two opinion dimension and gives 
# the agent on XY coordinates accordingly
function getplace!(turtle::Turtle)
    # check whether our cosen dimension is not bigger than maximum of dimensions in the simulation
    if data.xopinion > data.numberofopiniondimensions
        data.xopinion = 1
    end
    if data.yopinion > data.numberofopiniondimensions
        data.yopinion = 1
    end
    
    # then we rotate the agent towards the future place
    if typeof(turtle) == Agent
        facexy!(turtle, turtle.opinion[data.xopinion] * maxpxcor, turtle.opinion[data.yopinion] * maxpycor)
    end
    
    # lastly we move agent on the place given by opinion dimensions chosen for X and Y coordinates
    turtle.xcor = turtle.opinion[data.xopinion] * maxpxcor
    turtle.ycor = turtle.opinion[data.yopinion] * maxpycor
end

# sub routine for coloring agents according their average opinion across all dimensions --
# useful for distinguishing agents with same displayed coordinates, but differing in other opinion dimensions,
# then we see at one place agents with different colors.
function getcolor!(agent::Agent) 
    # speaking agents are colored from very dark red (average -1) 
    # through red (average 0) to very light red (average +1)
    agent.color = 15 + 4 * mean(agent.opinion) 
    agent.size = maxpxcor / 10
end

# sub-routine for visual purposes -- colors empty patches white, 
# patches with some agents light green, with many agents dark green, 
# with all agents black
function patchcolor(p::Patch)::Float64
    agentsherecnt::Int64 = 0
    for t in agents
        if t.xcor >= p.pxcor && t.xcor <= p.pxcor + 1.0 && 
            t.ycor >= p.pycor && t.ycor <= p.pycor + 1.0
            agentsherecnt += 1
        end
    end
    return 59.9 - (9.8 * (log(1 + agentsherecnt) / log(data.numberofagents)))
end

# Sub routine just for catching run-time errors
function avoidingruntimeerrors()
    # Check whether we set properly parameter 'updating' --
    # if we want update more dimensions than exists in simulation, then we set 'updating' to max of dimensions, i.e. 'opinions'
    if data.updating > data.numberofopiniondimensions 
        data.updating = data.numberofopiniondimensions
    end
end

function testhkbenchmark()
    data.maxticks = 5000
    data.rs = 1 # stepped 1:10
    data.numberofagents = 129 # 257 513
    data.boundarymean = 0.1 # 0.2 0.3
    data.xopinion = 1 
    data.yopinion = 1 
    data.setseed = true
    data.recordlength = 100
    data.numberofopiniondimensions = 1
end

function testhkbenchmarksos()
    data.maxticks = 5000
    data.rs = 1 # stepped 1:10
    data.numberofagents = 129 # 257 513
    data.boundarymean = 0.1 # 0.2 0.3
    data.xopinion = 1 
    data.yopinion = 1 
    data.setseed = true
    data.recordlength = 100
    data.numberofopiniondimensions = 1
end

function testhkbenchmarkv03()
    data.maxticks = 5000
    data.rs = 1 # stepped 1:73
    data.numberofagents = 129 # 257 513
    data.boundarymean = 0.1 # 0.2 0.3
    data.xopinion = 1 
    data.yopinion = 1 
    data.setseed = true
    data.recordlength = 100
    data.numberofopiniondimensions = 1
end

function testnewexperiment()
    data.maxticks = 1000
    data.rs = 101 # stepped 101:110
    data.numberofagents = 101 # 1001
    data.boundarymean = 0.1 # 0.2 0.3
    data.xopinion = 1 
    data.yopinion = 2
    data.setseed = true
    data.recordlength = 100
    data.numberofopiniondimensions = 1 # 8
end

# complement_ofHillClimbingSearchInClassicHK
function testRS001_200_complement()
    data.rs = 1 # 200
    data.setseed = true
    data.avoidseedcontrol = false
    data.boundarymean = 0.1 # -0.5
    data.conformitymean = 0.1 # -1
    data.numberofagents = 129
    data.numberofopiniondimensions = 1
    data.polarrepeats = 50
    data.yopinion = 1
    data.maxticks = 365
    data.centroidschange = 1.0E-5
    data.mean_identity_sigmoid_steepness = 0.85
    data.std_identity_sigmoid_steepness = 0.35
    data.opinionsigmoid = false
    data.polarisation_each_n_steps = 400
    data.showdicerolls = false
    data.esbgfurthestout = 5
    data.identitysigmoid = false
    data.maxsteepness = 700
    data.identitylevels = 7
    data.recordlength = 15
    data.centroidcolor = true
    data.useidentity = false
    data.mean_opinion_sigmoid_steepness = 0.45
    data.ncentroids = 3
    data.conformitystd = 0.15
    data.std_identity_sigmoid_xoffset = 0.1
    data.killingcentroids = true
    data.boundarystd = 0.05
    data.mean_identity_sigmoid_xoffset = 0.35
    data.normalizedistances = true
    data.std_opinion_sigmoid_steepness = 0.15
    data.identitytype = individual
    data.xopinion = 1
    data.minimumsdiro = 0.25
    data.maximumsdiro = 0.85
    data.sdiromean = 0.25
    data.sdirostd = 0.25
end

end # end of the module 