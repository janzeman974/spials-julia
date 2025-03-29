module simpleHK
using Random

@enum ShapeE house

import Random:
	rand, seed! 

import Statistics:
	mean, median

export
    Turtle, whocounter, HumanAb, Human, MediumAb, Medium, Patch, humans, patches, media, data, setup, go, ESBG, colorPatches!, changeAttitudes!, side # TODO posledni nevim


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


mutable struct Data
	humanSize::Float64 #/dole je slider min 0.5, max 5,  stred 2.0 a posouvam po 0.5
	n::Int64 #/dole je slider min 50, max 2000,  stred 1000 a posouvam po 50
	m::Int64 #/dole je slider min 0, max 5,  stred 5 a posouvam po 1 
	epsilon::Float64 #/dole je slider min 0.01, max 1,  stred 0.3 a posouvam po 0.01
	SD_epsilon::Float64 #/dole je slider min 0, max 0.5,  stred 0.3 a posouvam po 0.01
	
	Data(humanSize::Float64=2.0, n::Int64=1000,m::Int64=5,epsilon::Float64=0.3,SD_epsilon::Float64=0.3) = 
		new(humanSize, n, m, epsilon, SD_epsilon)
end
data::Data = Data()

abstract type Turtle end
#=struct Turtle{T} <: AbstractTurtle
	x::T,
	Turtle{T}(args...) where T = new{T}(T(args...))
end=#

whocounter::Int64 = 0  

abstract type HumanAb <: Turtle end  
mutable struct Human <: HumanAb
 	who::Int64
	color::Float64
	size::Float64
	xcor::Float64
	ycor::Float64

	heading::Float64
 
	attitude::Float64
	tolerance::Float64 

	Human(who::Int64=(global whocounter = whocounter + 1), color::Float64=0.5, size::Float64=1.0, xcor::Float64=25.0, ycor::Float64=25.0, heading::Float64=90.0, attitude::Float64=0.0, tolerance::Float64=0.0) = new(who, color, size, xcor, ycor, heading, attitude, tolerance)
end
humans::Vector{Human} = []
 
abstract type MediumAb <: Turtle end 
mutable struct Medium <: MediumAb
	who::Int64
	color::Float64
	size::Float64
	xcor::Float64
	ycor::Float64

	hidden::Bool
	shape::ShapeE
	attitude::Float64 
	tolerance::Float64

	Medium(who::Int64=(global whocounter = whocounter + 1), color::Float64=0.5, size::Float64=1.0, xcor::Float64=25.0, ycor::Float64=25.0, hidden::Bool=false, shape::ShapeE=house, attitude::Float64=0.0, tolerance::Float64=0.0) = new(who, color, size, xcor, ycor, hidden, shape, attitude, tolerance)
end
media::Vector{Medium} = []
 
function colorPatches!()
	# ;; Recoloring patches, agents, computing how model settled down
	# Patches update color according the number of turtles on it.	
	for patch in patches  
		agentsherecnt::Int64 = 0
    		for turtle in humans
        		if turtle.xcor >= patch.pxcor && turtle.xcor <= patch.pxcor + 1.0 && 
            		turtle.ycor >= patch.pycor && turtle.ycor <= patch.pycor + 1.0
            			agentsherecnt += 1
				 
 			end
           	end
         	patch.pcolor = Float64(agentsherecnt) / data.n * 10.0
	  	if patch.pcolor != 0.0
	 		  # print("colorpatches:[$(patch.pxcor),$(patch.pycor)]:$(agentsherecnt)/$(data.n)*10.0=$(patch.pcolor)")
		end
		# JZ: original commentary; if pcolor > 0.1 [ask neighbors [set pcolor [pcolor] of myself]]	 
	end
end

function side(human)
	if human.attitude > 0.0
		return 90.0
	else
		return 270.0
	end
end

function setup()
	# We erase the world and clean patches
	empty!(humans)
	global whocounter = 0
	empty!(media)

	# clear the canvas
	# turtlesCnt = 0

	# create n humans
	# TODO nenastavuju semeno
	local human::Human
	#print("before hum")
	for i::Int64 = 1 : data.n # JZ TODO check ma bejt vcetne 
		human = Human()
		human.size = data.humanSize
		human.attitude = 1.0 - rand() * 2
		human.xcor = human.attitude * 100
		human.heading = side(human)
		human.tolerance = 0.0
		while (human.tolerance <= 0.01 || human.tolerance > 2.0)
			human.tolerance = (randn() * data.SD_epsilon + data.epsilon) * 2.0 
		end
		push!(humans, human)
	end

	local medium::Medium
	# create m media
	for i::Int64 = 1 : data.m 
		medium = Medium()
		medium.size = 3 # JZ solid number
		medium.attitude = 1.0 - rand() * 2
		medium.xcor = medium.attitude * 100
		medium.color = 0.5 # JZ TODO enum nebo hexa orange
		medium.shape = house
		push!(media, medium)
	end

    colorPatches!()
	
	# JZ reset-ticks
	global ticks = 1
end

function go()
	for human in humans
		changeAttitudes!(human)
	end
	colorPatches!()
	global ticks = ticks + 1 
end

function ESBG()
	# computing centroids of two equal halves and diversity of each half
	# attitudes::Vector{Float64} = []
	# for human in humans
	#	push!(attitudes, human.attitude)
	# end
	# cutOff::Float64 = median(attitudes)

	g1::Vector{Human} = []
	g2::Vector{Human} = []
	minAtt::Float64 = -2.0 #pocitame za minimum
	minHum = nothing
	for i = 1 : data.n / 2
		# todo asi blbe, tady zakladam novej objekt
		minHum = nothing
		# get je na slovnik, ok
        for human in humans
         	if (human.attitude > minAtt) && (minHum == nothing || human.attitude < minHum.attitude)
                minHum = human
			end
			# JZ: a bit slow because it searches through all humans even those processed
                end
		if minHum != nothing
			push!(g1, minHum)
            minAtt = minHum.attitude
		end
	end

    for human in humans
		if in(human, g1) == false
            push!(g2, human)
		end
	end
	
	# JZ: mean of stg. inside an object cannot be done by Statistics.mean(
	meana::Float64 = 0.0
	for human in g1
		meana += human.attitude
	end
	meana /= length(g1)
        cent1::Float64 = meana / 2.0
	
    
	meana = 0.0
	for human in g2
		meana += human.attitude
	end
	meana /= length(g2)
        cent2::Float64 = meana / 2.0
        
    meana = 0.0
    	for human in g1
		meana += abs(human.attitude - cent1)
	end
	meana /= length(g1)
        dif1::Float64 = meana / 2.0
	meana = 0.0
    	for human in g2
		meana += abs(human.attitude - cent2)
	end
	meana /= length(g2)
        dif2::Float64 = meana / 2.0

	# returning ESBG
	return abs(cent1 - cent2) / (1 + dif1 + dif2)
end

function changeAttitudes!(human::Human)
	# computing perceived consensus of tolerated group of humans
	lowerBound::Float64 = human.attitude - human.tolerance
	upperBound::Float64 = human.attitude + human.tolerance
	# 7. 3. zakomentovano, jen kdyby to nefungovalo 
    # print("$(human.who)-lB:$(lowerBound),uB:$(upperBound)")
	tolerated::Vector{Turtle} = Turtle[]
	# JZ: turtles are humans and media
	for human2 in humans
		if (human2.attitude > lowerBound && human2.attitude < upperBound)
			push!(tolerated, human2)
		end
	end	
	for medium in media
		if (medium.attitude > lowerBound && medium.attitude < upperBound)
			push!(tolerated, medium)
		end
	end	
    groupConsensus::Float64 = 0.0
    if length(tolerated) != 0
        meana = 0.0
        for turt in tolerated
    		meana += turt.attitude
    	end
    	#print("meanasum $(meana)")
    	meana /= length(tolerated)
    	#print("lento $(length(tolerated))")
    	groupConsensus = meana
    end 
	#print("meana-cosensus $(meana)")

	# setting attitude
	human.attitude = groupConsensus
	human.xcor = human.attitude * 100.0
	#ycor also not in the original netlogo
	human.heading = side(human)
end
 
end	# module end 