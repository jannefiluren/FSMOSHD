# Idealized energy balance tests for the soil thermal routine

using FlexibleSnowModelOSHD
using Test

"""
    create_minimal_fsm(Tf, Ti)

Create a minimal FSM structure with only the fields needed for soil!()
"""
function create_minimal_fsm(Tf::Type, Ti::Type)

    # Initialize model with minimal configuration
    fsm = FSM{Tf, Ti}(
        Nx = 1,
        Ny = 1,
    )

    # Set tile fraction
    fsm.surface.tilefrac[1, 1] = Tf(1.0)

    # Set default thermal properties
    # Typical soil volumetric heat capacity: ~2.0e6 J/m³/K
    # csoil in the model is areal heat capacity = volumetric × thickness
    for k in 1:fsm.grid.Nsoil
        fsm.diag.csoil[k, 1, 1] = Tf(2.0e6) * fsm.grid.Dzsoil[k]  # J/m²/K
        fsm.diag.ksoil[k, 1, 1] = Tf(1.5)                    # W/m/K (typical soil)
    end

    # Initialize temperatures to a reasonable value
    fsm.state.Tsoil[:, 1, 1] .= Tf(285.0)

    # Initialize Gsoil (will be set in each test)
    fsm.diag.Gsoil[1, 1] = Tf(0.0)

    return fsm
end

"""
    compute_soil_energy(fsm)

Compute total thermal energy stored in the soil column.

Energy is computed relative to absolute zero:
    E = Σ_k [csoil[k] × Tsoil[k]]

where csoil[k] is the areal heat capacity (J/m²/K) of layer k.
"""
function compute_soil_energy(fsm)
    E = 0.0
    for k in 1:fsm.grid.Nsoil
        E += fsm.diag.csoil[k, 1, 1] * fsm.state.Tsoil[k, 1, 1]
    end
    return E
end

@testset "Soil Energy Balance Tests" begin

    # Common parameters
    Tf = Float32
    Ti = Int32

    # Test 1: Uniform temperature, positive heat flux (warming from above)
    @testset "Warming from above" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set uniform initial temperature
        T_init = 280.0  # K
        fsm.state.Tsoil[:, 1, 1] .= T_init

        # Set positive heat flux into soil (warming)
        Gsoil_in = 50.0  # W/m²
        fsm.diag.Gsoil[1, 1] = Gsoil_in

        # Calculate initial energy
        E_initial = compute_soil_energy(fsm)

        # Run soil routine
        soil!(fsm)

        # Calculate final energy
        E_final = compute_soil_energy(fsm)

        # Expected energy change from top boundary flux
        E_expected_from_top = Gsoil_in * fsm.params.dt
        E_actual_change = E_final - E_initial

        # The energy change should be very close to the top flux
        # (small deviations due to bottom boundary treatment)
        @test isapprox(E_actual_change, E_expected_from_top, rtol = 1.0e-4)

        # Check that top layer warmed (physics sanity check)
        @test fsm.state.Tsoil[1, 1, 1] > T_init
    end

    # Test 2: Uniform temperature, negative heat flux (cooling from above)
    @testset "Cooling from above" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set uniform initial temperature
        T_init = 290.0  # K
        fsm.state.Tsoil[:, 1, 1] .= T_init

        # Set negative heat flux (heat leaving soil)
        Gsoil_in = -30.0  # W/m²
        fsm.diag.Gsoil[1, 1] = Gsoil_in

        # Calculate initial energy
        E_initial = compute_soil_energy(fsm)

        # Run soil routine
        soil!(fsm)

        # Calculate final energy
        E_final = compute_soil_energy(fsm)

        # Expected energy change from top boundary flux
        E_expected_from_top = Gsoil_in * fsm.params.dt
        E_actual_change = E_final - E_initial

        # The energy change should be very close to the top flux
        @test isapprox(E_actual_change, E_expected_from_top, rtol = 1.0e-4)

        # Check that top layer cooled (physics sanity check)
        @test fsm.state.Tsoil[1, 1, 1] < T_init
    end

    # Test 3: Zero heat flux at top, uniform temperature (should be stable)
    @testset "Isothermal stability" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set uniform temperature at absolute reference (Tm = 273.15 K)
        # This minimizes the bottom boundary flux effect
        T_init = 273.15  # K (melting point - used as reference)
        fsm.state.Tsoil[:, 1, 1] .= T_init

        # Set zero heat flux at top
        Gsoil_in = 0.0  # W/m²
        fsm.diag.Gsoil[1, 1] = Gsoil_in

        # Run soil routine
        soil!(fsm)

        # Check that temperatures didn't change much (isothermal should stay isothermal)
        for k in 1:fsm.grid.Nsoil
            @test isapprox(fsm.state.Tsoil[k, 1, 1], T_init, atol = 1.0e-4)
        end
    end

    # Test 4: Non-uniform thermal properties
    @testset "Non-uniform thermal properties" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set varying thermal properties by layer (areal heat capacity)
        fsm.diag.csoil[1, 1, 1] = 1.5e5  # J/m²/K
        fsm.diag.csoil[2, 1, 1] = 4.0e5  # J/m²/K
        fsm.diag.csoil[3, 1, 1] = 1.0e6  # J/m²/K
        fsm.diag.csoil[4, 1, 1] = 1.6e6  # J/m²/K

        fsm.diag.ksoil[1, 1, 1] = 1.0  # W/m/K
        fsm.diag.ksoil[2, 1, 1] = 1.5  # W/m/K
        fsm.diag.ksoil[3, 1, 1] = 2.0  # W/m/K
        fsm.diag.ksoil[4, 1, 1] = 1.5  # W/m/K

        # Set temperature gradient
        fsm.state.Tsoil[1, 1, 1] = 285.0
        fsm.state.Tsoil[2, 1, 1] = 283.0
        fsm.state.Tsoil[3, 1, 1] = 281.0
        fsm.state.Tsoil[4, 1, 1] = 280.0

        # Set heat flux
        Gsoil_in = 25.0  # W/m²
        fsm.diag.Gsoil[1, 1] = Gsoil_in

        # Calculate initial energy
        E_initial = compute_soil_energy(fsm)

        # Run soil routine
        soil!(fsm)

        # Calculate final energy
        E_final = compute_soil_energy(fsm)

        # Expected energy change from top boundary flux
        E_expected_from_top = Gsoil_in * fsm.params.dt
        E_actual_change = E_final - E_initial

        # Energy change should be approximately equal to top flux
        # (allow slightly larger tolerance due to temperature gradient)
        @test isapprox(E_actual_change, E_expected_from_top, rtol = 1.0e-3)
    end

    # Test 5: Multiple timesteps with cumulative energy tracking
    @testset "Multiple timesteps" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set initial conditions
        fsm.state.Tsoil[:, 1, 1] .= 275.0
        Gsoil_in = 20.0  # W/m²
        fsm.diag.Gsoil[1, 1] = Gsoil_in

        # Calculate initial energy
        E_initial = compute_soil_energy(fsm)

        # Run multiple timesteps
        n_steps = 10
        for _ in 1:n_steps
            soil!(fsm)
        end

        # Calculate final energy
        E_final = compute_soil_energy(fsm)

        # Expected total energy change from top boundary
        E_expected_from_top = Gsoil_in * fsm.params.dt * n_steps
        E_actual_change = E_final - E_initial

        # Check energy conservation over multiple steps
        @test isapprox(E_actual_change, E_expected_from_top, rtol = 1.0e-4)
    end

    # Test 6: Physics sanity checks
    @testset "Physics sanity checks" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set cold top, warm bottom (inverted gradient)
        fsm.state.Tsoil[1, 1, 1] = 270.0
        fsm.state.Tsoil[2, 1, 1] = 275.0
        fsm.state.Tsoil[3, 1, 1] = 280.0
        fsm.state.Tsoil[4, 1, 1] = 285.0

        # Zero flux at top - heat should flow upward internally
        fsm.diag.Gsoil[1, 1] = 0.0

        T_top_init = fsm.state.Tsoil[1, 1, 1]

        # Run soil routine
        soil!(fsm)

        # Top layer should warm (heat flows from warm bottom to cold top)
        @test fsm.state.Tsoil[1, 1, 1] > T_top_init

        # Temperature gradient should decrease (smoothing)
        gradient_after = fsm.state.Tsoil[4, 1, 1] - fsm.state.Tsoil[1, 1, 1]
        @test gradient_after < 15.0  # Was 15 K initially
    end

    # Test 7: Set zero conductivity and only warm the top soil layer
    @testset "Warming of top layer" begin
        fsm = create_minimal_fsm(Tf, Ti)

        # Set uniform initial temperature
        T_init = 280.0  # K
        fsm.state.Tsoil[:, 1, 1] .= T_init

        # Set zero conductivity
        fsm.diag.ksoil[:, 1, 1] .= 0.0  # W/m/K

        # Add a heat flux to increase the top layer soil temperature by one degree
        fsm.diag.Gsoil[1, 1] = 2.0e5 / fsm.params.dt

        # Run soil routine
        soil!(fsm)

        # Top soil layer temperature should have increased by one degree
        @test isapprox(fsm.state.Tsoil[1, 1, 1], T_init + 1, atol = 1.0e-4)
    end
end
