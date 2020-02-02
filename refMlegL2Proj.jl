#=
refMlegL2Proj.jl

ref = reference element
M = Modal
leg = Legendre

A program that returns the coordinates of the L2 projection of a function onto the space spanned by three legendre polynomials on the reference interval [-1, 1].

This amounts too finding uᶜ such that:

               ∑ (ϕᵢ, ϕⱼ) uᶜᵢ = (u , ϕᵢ)     ∀i


since these are legendre (aka orthanormal) M = ∑ (ϕᵢ, ϕⱼ) ∀i = I. Or:

Uᶜᵢ = ∫ uϕᵢ dx ≈ ∑ u(xⱼ)ϕᵢ(xⱼ)wⱼ

Which means that this really just amounts to finding a good (gaussian) quadtrature for uϕᵢ as seen above.
We don't know anything about u (if u is the solution to unsolved diff eq), but we know that ϕᵢ is an order i polynomial. 


Since we have 3 legendre polynomials right now our system looks like:

  (u,ϕ₁) ≈ u(x₁)ϕ₁(x₁)w₁ +  u(x₂)ϕ₁(x₂)w₂ + u(x₃)ϕ₁(x₃)w₃
  (u,ϕ₂) ≈ u(x₁)ϕ₂(x₁)w₁ +  u(x₂)ϕ₂(x₂)w₂ + u(x₃)ϕ₂(x₃)w₃
  (u,ϕ₃) ≈ u(x₁)ϕ₃(x₁)w₁ +  u(x₂)ϕ₃(x₂)w₂ + u(x₃)ϕ₃(x₃)w₃
  
=#


u(x) = x^3

ϕ₁(x) = 1/sqrt(2)
ϕ₂(x) = sqrt(3/2)*x
ϕ₃(x) = 3*sqrt(5)/2*sqrt(2)*x^2 - sqrt(5)/2*sqrt(2)
