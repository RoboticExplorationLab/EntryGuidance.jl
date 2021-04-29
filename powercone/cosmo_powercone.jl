using LinearAlgebra, COSMO, JuMP, MathOptInterface
const MOI = MathOptInterface

# model = COSMO.Model()
#
# # variable x = [t v1 v2 v3]
#
# # first we stup the SOC constraint with t
# cs1 = COSMO.Constraint(Array(float(I(4))), zeros(4), COSMO.SecondOrderCone)
#
# # now we do the power stuff
c = 342.0
α = 1/4
# cs2 = COSMO.Constraint([0 0 0 0 ;0 0 0 0  ; 1.0 0 0 0], [c;1;0],COSMO.PowerCone(α))
# P = randn(4,4);P = P'*P + 1e-3*I;
# P = float(zeros(4,4))
#
# custom_settings = COSMO.Settings()
# custom_settings.verbose = true
# COSMO.assemble!(model, P, [0;randn(3)], [cs1;cs2],settings = custom_settings)
# results = COSMO.optimize!(model);

m = JuMP.Model(COSMO.Optimizer)

@variable(m,t)
@variable(m,v[1:3])

@constraint(m, [t,v[1],v[2],v[3]] in SecondOrderCone())
@constraint(m, [c,1,t] in MOI.PowerCone(α) )

@objective(m, Min, dot(randn(3),v))
optimize!(m)
