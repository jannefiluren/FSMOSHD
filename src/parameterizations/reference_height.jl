# Reference height
struct AboveGround{Tf} <: AbstractReferenceHeight{Tf} end
struct AboveCanopy{Tf} <: AbstractReferenceHeight{Tf} end
AboveGround{Tf}(grid::Grid; kwargs...) where {Tf} = AboveGround{Tf}()
AboveCanopy{Tf}(grid::Grid; kwargs...) where {Tf} = AboveCanopy{Tf}()

@inline reference_heights(::AboveGround, zU, zT, hcan) = (zU, zT)
@inline reference_heights(::AboveCanopy, zU, zT, hcan) = (zU + hcan, zT + hcan)
