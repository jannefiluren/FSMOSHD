function tridiag!(x, Nvec, gamma, Nmax, a, b, c, r)

  fill!(gamma, zero(Float64))

  beta = b[1]
  x[1] = r[1] / beta

  for n = 2:Nvec
    gamma[n] = c[n-1] / beta
    beta = b[n] - a[n] * gamma[n]
    x[n] = (r[n] - a[n] * x[n-1]) / beta
  end

  for n = Nvec-1:-1:1
    x[n] = x[n] - gamma[n+1] * x[n+1]
  end

end