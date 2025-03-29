using simpleHK
using Statistics
using Test

#include("groupformation.jl")
#import .GroupFormation would import only module name
 
@testset "setuptests" begin
	whocounter = 0
	empty!(humans)
	empty!(media) 
	t1::Human=Human()
 	@test t1.who == 1
	t2::Human=Human()
	@test t2.who == 2
	push!(humans, t1)
	push!(humans, t2) 
	t3::Medium=Medium() 	
    
	@testset "turtlestest" begin
         	@test t3.who == t1.who+2
	end

	#patches are tiles of the current field
	#therefore colorpatches! counts again the whole field, each tile 

 	@testset "colorpatchestest" begin
		#data.n = 10
		t1.xcor = 2.5
		t1.ycor = 2.5
		colorPatches!()
		# @test patches == [0, 0, 0; 0, 0.01, 0; 0, 0, 0]
		@test patches[-simpleHK.minpxcor + 1 + 2, -simpleHK.minpycor + 1 + 2].pcolor == 0.01
		
		t2.xcor = 2.5
		t2.ycor = 2.5
		colorPatches!()
		patches[-simpleHK.minpxcor + 1 + 2, -simpleHK.minpycor + 1 + 2].pcolor == 0.02
	end
 	
	@testset "setupmethodtest" begin
		empty!(humans)
		empty!(media)		
		setup()
         	@test isempty(humans) == false		
		@test isempty(media) == false

		@test humans[1].attitude > -1.0 
		@test humans[1].attitude <  1.0
		@test humans[1].tolerance > 0.0
		@test humans[1].tolerance < 2.0
		@test media[1].attitude > -1.0        
		@test media[1].attitude <  1.0 

		#TODO check
        	@test ESBG() != 0.0
		#OT udela ESBG neco s media? ne, pouze navrati diference
		# heading nelze otestovat, attitude muze byt < i > 0	
		#@test humans[1].heading == 270.0
	end

	@testset "esbgtest" begin
		#data.n = 50
		#setup()	
		#FOR DEBUGGING UNCOMMENT THIS AND PRINTLNS IN SIMPLEHK.ESBG() 
		@info "$(humans[1].attitude)"
		@info "$(humans[2].attitude)"
		@info "$(humans[3].attitude)"
		ESBG() 
	end

	@testset "changeattitudestest" begin
		t1.attitude = 1.0 - rand() * 2
		t1.xcor = 2.5
		t1.heading = side(t1)
		t1.tolerance = 0.0
		while (t1.tolerance <= 0.01 || t1.tolerance > 2.0)
			t1.tolerance = (randn() * data.SD_epsilon + data.epsilon) * 2.0 
		end 
		changeAttitudes!(t1) # changes t1.attitude, t1.xcor, t1.heading
		pomAtt = t1.attitude
		#print("calling t1 for the 2nd time")
 		changeAttitudes!(t1)
		@test t1.attitude != pomAtt
		pomAtt = t1.attitude
		# the same coordinates of t2=>the same attitude
		#print("calling t2")
		t2.attitude = 1.0 - rand() * 2
		t2.xcor = 2.5
		t2.heading = side(t2)
		t2.tolerance = 0.0
		while (t2.tolerance <= 0.01 || t2.tolerance > 2.0)
			t2.tolerance = (randn() * data.SD_epsilon + data.epsilon) * 2.0 
		end
		pomAtt = t2.attitude
		changeAttitudes!(t2)
		#with introduction of determinity, this will  work 
		@test t2.attitude != pomAtt
	end 
end

@testset "gotest" begin
	empty!(humans)
	empty!(media)		
  	maxticks::Int64 = 100 # 
	data.humanSize::Float64 = 2.0
	# JZ number of humans
	data.n::Int64 = 50 # 150, 1000 
	# JZ number of media
	data.m::Int64 = 0 # 0, 1, 2, 3, 4, 5
	data.epsilon::Float64 = 0.1 # 0.1, 0.2, 0.3
	data.SD_epsilon::Float64 = 0 # 0, 0.1, 0.2, 0.3
	setup()	
	for i::Int64 = 1 : maxticks
		go()
		ESBG()
	end
	# co testovat?
end
  