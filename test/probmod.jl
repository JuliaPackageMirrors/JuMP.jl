#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/probmod.jl
# Testing problem modification
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, FactCheck, BaseTestNext

# If solvers not loaded, load them (i.e running just these tests)
!isdefined(:lp_solvers) && include("solvers.jl")

const TOL = 1e-4
approx_eq(x::Variable, y) = isapprox(getValue(x), y, atol=TOL)
approx_eq(m::Model, y) = isapprox(getObjectiveValue(m), y, atol=TOL)


@testset "Problem modification" begin
print_with_color(:magenta, "Problem modification\n")

@testset "LP Model (with $typeof(solver))" for solver in lp_solvers
    # max 1.1x + 1.0y
    # st     x +    y <= 3
    #     0 <= x <= 3
    #     1 <= y <= 3
    # x* = 2, y* = 1
    m = Model(solver=solver)
    @defVar(m, 0 <= x <= 3)
    @defVar(m, 1 <= y <= 3)
    @setObjective(m, Max, 1.1x + 1.0y)
    maincon = @addConstraint(m, x + y <= 3)
    @test solve(m) == :Optimal
    @test m.internalModelLoaded
    @test approx_eq(x, 2.0)
    @test approx_eq(y, 1.0)

    # Test adding a variable
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     1 <= y <= 3
    #     0 <= z <= 5
    # x* = 0, y* = 1, z* = 2
    @defVar(m, 0 <= z <= 5,
            objective=100.0,
            inconstraints=[maincon],
            coefficients=[1.0])
    @test solve(m) == :Optimal
    @test approx_eq(x, 0.0)
    @test approx_eq(y, 1.0)
    @test approx_eq(z, 2.0)
    @test approx_eq(m, 201.0)

    # Test changing bounds
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 2
    # x* = 1, y* = 0, z* = 2
    setLower(y, 0.0)
    setUpper(z, 2.0)
    @test solve(m) == :Optimal
    @test approx_eq(x, 1.0)
    @test approx_eq(y, 0.0)
    @test approx_eq(z, 2.0)
    m.internalModelLoaded = false
end  # "LP Model (with $solver)"


@testset "IP Model (with $typeof(solver))" for solver in ip_solvers
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     1 <= y <= 3
    #     0 <= z <= 0
    m = Model(solver=solver)
    @defVar(m, 0 <= x <= 3)
    @defVar(m, 0 <= y <= 3)
    @defVar(m, 0 <= z <= 0)
    @setObjective(m, Max, 1.1x + 1.0y + 100.0z)
    maincon = @addConstraint(m, x + y + z <= 3)
    @test solve(m) == :Optimal

    # Test changing problem type
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 3
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 2, y* = 0, z* = 1
    setUpper(z, 1.5)
    setCategory(z, :Int)
    @test solve(m) == :Optimal
    @test approx_eq(x, 2.0)
    @test approx_eq(y, 0.0)
    @test approx_eq(z, 1.0)

    # Test changing constraint bound (<= constraint)
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 1, y* = 0, z* = 1
    chgConstrRHS(maincon, 2.0)
    @test solve(m) == :Optimal
    @test approx_eq(x, 1.0)
    @test approx_eq(y, 0.0)
    @test approx_eq(z, 1.0)

    # Test adding a constraint
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #        x        +      z <= 0
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 0, y* = 2, z* = 0
    xz0ref = @addConstraint(m, x + z <=0)
    @test solve(m) == :Optimal
    @test approx_eq(x, 0.0)
    @test approx_eq(y, 2.0)
    @test approx_eq(z, 0.0)

    # Test changing constraint bound (>= constraint)
    # max 1.1x + 1.0y + 100.0z
    # st     x +    y +      z <= 2
    #        x        +      z <= 1
    #        x +    y          >= 1
    #     0 <= x <= 3
    #     0 <= y <= 3
    #     0 <= z <= 1.5, Integer
    # x* = 1, y* = 0, z* = 1
    chgConstrRHS(xz0ref, 2.0)
    xyg0ref = @addConstraint(m, x + y >= 0)
    @test solve(m) == :Optimal
    chgConstrRHS(xyg0ref, 1)
    @test solve(m) == :Optimal
    @test approx_eq(x, 1.0)
    @test approx_eq(y, 0.0)
    @test approx_eq(z, 1.0)
end  # "IP Model (with $solver)"

end  # "Problem modification"


facts("[probmod] Test adding a range constraint and modifying it") do
    m = Model()
    @defVar(m, x)
    rangeref = @addConstraint(m, -10 <= x <= 10)
    @fact_throws chgConstrRHS(rangeref, 11)
end


facts("[probmod] Test adding a 'decoupled' variable (#205)") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, x >= 0)
    @setObjective(m, Min, x)
    solve(m)
    @defVar(m, y >= 0)
    @addConstraint(m, x + y == 1)
    @setObjective(m, Min, 2x+y)
    solve(m)
    @fact getValue(x) --> roughly(0.0, TOL)
    @fact getValue(y) --> roughly(1.0, TOL)
end
end
end


facts("[probmod] Test buildInternalModel") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, x >= 0)
    @defVar(m, y >= 0)
    @addConstraint(m, x + y == 1)
    @setObjective(m, Max, y)
    buildInternalModel(m)
    @fact getInternalModel(m) --> not(nothing)
    @fact m.internalModelLoaded --> true
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getValue(x) --> roughly( 0.0, TOL)
    @fact getValue(y) --> roughly( 1.0, TOL)
    @fact getObjectiveValue(m) --> roughly(1.0, TOL)
    @fact getDual(x)  --> roughly(-1.0, TOL)
    @fact getDual(y)  --> roughly( 0.0, TOL)
end
end
end


facts("[probmod] Test buildInternalModel with MIP") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, x >= 0, Int)
    @defVar(m, y, Bin)
    @addConstraint(m, x + y == 1)
    @setObjective(m, Max, y)
    buildInternalModel(m)
    @fact getInternalModel(m) --> not(nothing)
    @fact m.internalModelLoaded --> true
    stat = solve(m)
    @fact stat --> :Optimal
    @fact getValue(x) --> roughly( 0.0, TOL)
    @fact getValue(y) --> roughly( 1.0, TOL)
    @fact getObjectiveValue(m) --> roughly(1.0, TOL)
end
end
end


facts("[probmod] Test adding empty constraints") do
for solver in lp_solvers
context("With solver $(typeof(solver))") do
    m = Model(solver=solver)
    @defVar(m, 0 <= x <= 9)
    @setObjective(m, Max, x)
    @addConstraint(m, x <= 5)
    @addConstraint(m, 0 <= 1)
    solve(m)
    @addConstraint(m, 0 <= 1)
    @fact solve(m) --> :Optimal
    @fact getValue(x) --> roughly(5.0, TOL)
end
end
end


facts("[probmod] Test bound modification on binaries") do
for solver in ip_solvers
context("With solver $(typeof(solver))") do
    # Test semantics for modifying bounds on binary variables:
    # Variables should be restricted to the intersection of
    # {0,1} and their bounds.
    mod = Model(solver=solver)
    @defVar(mod, x, Bin)
    @setObjective(mod, Max, x)
    solve(mod)
    @fact getValue(x) --> roughly(1.0)
    setUpper(x, 2.0)
    solve(mod)
    @fact getValue(x) --> roughly(1.0)
    setUpper(x, 0.0)
    solve(mod)
    @fact getValue(x) --> roughly(0.0)
    # same thing, other direction
    mod = Model(solver=solver)
    @defVar(mod, x, Bin)
    @setObjective(mod, Min, x)
    solve(mod)
    @fact getValue(x) --> roughly(0.0)
    setLower(x, -1.0)
    solve(mod)
    @fact getValue(x) --> roughly(0.0)
    setLower(x, 1.0)
    solve(mod)
    @fact getValue(x) --> roughly(1.0)
end
end
end

facts("[probmod] Applicable regressions") do

function methods_test(solvername, solverobj, supp)
    mod = Model(solver=solverobj)
    @defVar(mod, x >= 0)
    @addConstraint(mod, 2x == 2)
    solve(mod, suppress_warnings=true)
    internal_mod = getInternalModel(mod)
    context(solvername) do
        for (it,(meth, args)) in enumerate(mpb_methods)
            if supp[it]
                @fact applicable(meth, internal_mod, args...) --> true
                @fact method_exists(meth, map(typeof, tuple(internal_mod, args...))) --> true
            end
        end
    end
end

const mpb_methods =
    [(MathProgBase.addquadconstr!, (Cint[1],Float64[1.0],Cint[1],Cint[1],Float64[1],'>',1.0)),
     (MathProgBase.setquadobjterms!, (Cint[1], Cint[1], Float64[1.0])),
     (MathProgBase.addconstr!,   ([1],[1.0],1.0,1.0)),
     (MathProgBase.addsos1!,     ([1],[1.0])),
     (MathProgBase.addsos2!,     ([1],[1.0])),
     (MathProgBase.addvar!,      ([1],[1.0],1.0,1.0,1.0)),
     (MathProgBase.setvarLB!,    ([1.0],)),
     (MathProgBase.setvarUB!,    ([1.0],)),
     (MathProgBase.setconstrLB!, ([1.0],)),
     (MathProgBase.setconstrUB!, ([1.0],)),
     (MathProgBase.setobj!,      ([1.0],)),
     (MathProgBase.setsense!,    (:Min,)),
     (MathProgBase.setvartype!,  ([:Cont],)),
     (MathProgBase.getinfeasibilityray, ()),
     (MathProgBase.getunboundedray, ()),
     (MathProgBase.getreducedcosts, ()),
     (MathProgBase.getconstrduals, ()),
     (MathProgBase.setwarmstart!, ([1.0]))]

if grb
    supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
    methods_test("Gurobi", Gurobi.GurobiSolver(OutputFlag=0), supp)
end
if cpx
    supp = (true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true)
    methods_test("CPLEX", CPLEX.CplexSolver(), supp)
end
if cbc
    supp = (false,false,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
    methods_test("Clp", Clp.ClpSolver(), supp)
end
if cbc
    supp = (false,false,false,false,false,false,false,false,false,false,false,true,true,false,false,false,false,false)
    methods_test("Cbc", Cbc.CbcSolver(), supp)
end
if glp
    supp = (false,false,true,false,false,true,true,true,true,true,true,true,false,true,true,true,true,false)
    methods_test("GLPK LP", GLPKMathProgInterface.GLPKSolverLP(), supp)
    supp = (false,false,true,false,false, true,true,true,true,true,true,true,true,false,false,false,false,false)
    methods_test("GLPK MIP", GLPKMathProgInterface.GLPKSolverMIP(), supp)
end
if mos
    supp = (true,true,true,false,false,true,true,true,true,true,true,true,true,true,true,true,true,false)
    methods_test("Mosek", Mosek.MosekSolver(), supp)
end

end
